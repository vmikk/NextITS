/*
============================================================================
  NextITS: Pipeline to process eukaryotic ITS amplicons
============================================================================
  License: Apache-2.0
  Github : https://github.com/vmikk/NextITS
  Website: https://Next-ITS.github.io/
----------------------------------------------------------------------------
*/

// Subworkflow for primer trimming and extraction of rRNA regions
//
// The workflow is as follows:
//   1. Trim primers and dereplicate at sample level
//   2. Split the dereplicated primer-trimmed sequences (at sample level) into chunks, while preserving metadata
//   3. Run the ITS extractor (ITSx 1.x or ITSx2) on each chunk, in parallel
//   4. Group results back by sample ID, concatenate + convert to Parquet
//   5. Pool the extracted regions across all samples
//
// If `params.its_region == "none"`, steps 2-5 are skipped and the dereplicated primer-trimmed sequences are used for the downstream analysis
//
// Extractor selection (`params.itsx_tool`):
//   "ITSx"  - ITSx v1.x  (HMMER-based, slow, supports taxonomic profiles and partial regions)
//   "ITSx2" - ITSx2      (Infernal covariance models, much faster; no partial regions, no taxonomic profiles, no `problematic`/`extraction.results` output)

// Trim primers and dereplicate at sample level
process primer_trim {

    label "main_container"

    publishDir(
      params.its_region == "none" ? "${params.outdir}/03_PrimerTrim" : "${params.outdir}/03_ITSx",
      mode:   "${params.storagemode}",
      saveAs: { fn -> params.its_region == "none"
                        ? fn
                        : fn.replaceAll(/\.fa\.gz$/, '_derep.fasta.gz') }
    )
    // cpus 2

    // Add sample ID to the log file
    tag "${meta.id}"

    input:
      tuple val(meta), path(fastq)

    output:
      tuple val(meta), path("${meta.id}.fa.gz"),             emit: derep,  optional: true
      tuple val(meta), path("${meta.id}_hash_table.txt.gz"), emit: hashes, optional: true
      tuple val(meta), path("${meta.id}_uc.uc.gz"),          emit: uc,     optional: true
      tuple val(meta), path("${meta.id}_primertrimmed_sorted.fq.gz"), emit: trimmed_seqs, optional: true
      tuple val("${task.process}"), val('cutadapt'), eval('cutadapt --version'), topic: versions
      tuple val("${task.process}"), val('vsearch'), eval('vsearch --version 2>&1 | head -n 1 | sed "s/vsearch //g" | sed "s/,.*//g" | sed "s/^v//" | sed "s/_.*//"'), topic: versions
      tuple val("${task.process}"), val('seqkit'), eval('seqkit version | sed "s/seqkit v//"'), topic: versions
      tuple val("${task.process}"), val('phredsort'), eval('phredsort -v | sed "s/phredsort //"'), topic: versions
      tuple val("${task.process}"), val('seqhasher'), eval('seqhasher -v | sed "s/SeqHasher //"'), topic: versions
      tuple val("${task.process}"), val('parallel'), eval('parallel --version | head -n 1 | sed "s/GNU parallel //"'), topic: versions

    script:
    def sampID = "${meta.id}"
    """
    echo -e "Primer trimming and dereplication at sample level\\n"
    echo -e "Input sample: "   ${sampID}
    echo -e "Forward primer: " ${params.primer_forward}
    echo -e "Reverse primer: " ${params.primer_reverse}

    ## Reverse-complement rev primer
    RR=\$(rc.sh ${params.primer_reverse})
    echo -e "Reverse primer RC: " "\$RR"

    ## Trim primers
    echo -e "\\nTrimming primers"
    cutadapt \
      -a ${params.primer_forward}";required;min_overlap=${params.primer_foverlap}"..."\$RR"";required;min_overlap=${params.primer_roverlap}" \
      --errors ${params.primer_mismatches} \
      --revcomp --rename "{id}" \
      --discard-untrimmed \
      --minimum-length ${params.trim_minlen} \
      --cores ${task.cpus} \
      --action trim \
      --output ${sampID}_primertrimmed.fq.gz \
      ${fastq}

    echo -e "..Done\\n"

    ## Check if there are sequences in the output
    NUMSEQS=\$( seqkit stat --tabular --quiet ${sampID}_primertrimmed.fq.gz | awk -F'\t' 'NR==2 {print \$4}' )
    echo -e "Number of sequences after primer trimming: " \$NUMSEQS
    if [ \$NUMSEQS -lt 1 ]; then
      echo -e "\\nIt looks like no reads remained after trimming the primers\\n"
      rm -f ${sampID}_primertrimmed.fq.gz
      exit 0
    fi

    ## Estimate sequence quality and sort sequences by quality
    echo -e "\\nSorting by sequence quality"
    seqkit replace -p "\\s.+" ${sampID}_primertrimmed.fq.gz \
      | phredsort -i - -o - --metric meep --header avgphred,maxee,meep \
      | gzip -${params.gzip_compression} > ${sampID}_primertrimmed_sorted.fq.gz
    echo -e "..Done"

    ## Remove the intermediate file as early as possible (it can be large)
    rm -f ${sampID}_primertrimmed.fq.gz

    ## Hash sequences, add sample ID to the header
    ## columns: Sample ID - Hash - PacBioID - AvgPhredScore - MaxEE - MEEP - Sequence - Quality - Length
    echo -e "\\nCreating hash table"
    seqhasher --hash sha1 --name ${sampID} ${sampID}_primertrimmed_sorted.fq.gz - \
      | seqkit fx2tab --length \
      | sed 's/;/\t/ ; s/;/\t/ ; s/ avgphred=/\t/ ; s/ maxee=/\t/ ; s/ meep=/\t/' \
      > ${sampID}_hash_table.txt
    echo -e "..Done"

    ## Dereplicate at sample level
    ## (use quality-sorted sequences, so that the representative sequence is the one with the highest quality)
    echo -e "\\nDereplicating at sample level"
    seqkit fq2fa -w 0 ${sampID}_primertrimmed_sorted.fq.gz \
      | vsearch \
        --derep_fulllength - \
        --output - \
        --strand both \
        --fasta_width 0 \
        --threads 1 \
        --relabel_sha1 \
        --sizein --sizeout \
        --minseqlength ${params.trim_minlen} \
        --uc ${sampID}_uc.uc \
        --quiet \
      > ${sampID}.fa

    echo -e "..Done"

    ## Compress results
    echo -e "\\nCompressing results"
    parallel -j${task.cpus} "gzip -${params.gzip_compression} {}" ::: \
      ${sampID}_hash_table.txt \
      ${sampID}_uc.uc \
      ${sampID}.fa

    echo -e "..Done"
    """
}


// Extract rRNA regions with ITSx (v1.x) from a single chunk of dereplicated sequences
// NB. In input data, sequence header should not contain spaces!
process itsx {

    label "main_container"

    // Chunk-level results are not published - they are concatenated per sample first
    // cpus 3

    tag "${meta.id}__chunk${meta.chunk_id}"

    input:
      tuple val(meta), path(input)   // FASTA file with dereplicated sequences

    output:
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}.full.fasta.gz"), emit: itsx_full, optional: true
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}.SSU.fasta.gz"),  emit: itsx_ssu,  optional: true
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}.ITS1.fasta.gz"), emit: itsx_its1, optional: true
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}.5_8S.fasta.gz"), emit: itsx_58s,  optional: true
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}.ITS2.fasta.gz"), emit: itsx_its2, optional: true
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}.LSU.fasta.gz"),  emit: itsx_lsu,  optional: true
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}.positions.txt"),   emit: itsx_positions,   optional: true
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}.problematic.txt"), emit: itsx_problematic, optional: true
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}_no_detections.fasta.gz"), emit: itsx_nondetects,     optional: true
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}_no_detections.txt"),      emit: itsx_nondetects_txt, optional: true
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}.summary.txt"),           emit: itsx_summary, optional: true
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}.extraction.results.gz"), emit: itsx_details, optional: true
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}.SSU.full_and_partial.fasta.gz"),  emit: itsx_ssu_part,  optional: true
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}.ITS1.full_and_partial.fasta.gz"), emit: itsx_its1_part, optional: true
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}.5_8S.full_and_partial.fasta.gz"), emit: itsx_58s_part,  optional: true
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}.ITS2.full_and_partial.fasta.gz"), emit: itsx_its2_part, optional: true
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}.LSU.full_and_partial.fasta.gz"),  emit: itsx_lsu_part,  optional: true
      tuple val("${task.process}"), val('ITSx'), eval('ITSx --help 2>&1 | head -n 3 | tail -n 1 | sed "s/Version: //"'), topic: versions
      tuple val("${task.process}"), val('seqkit'), eval('seqkit version | sed "s/seqkit v//"'), topic: versions
      tuple val("${task.process}"), val('parallel'), eval('parallel --version | head -n 1 | sed "s/GNU parallel //"'), topic: versions
      tuple val("${task.process}"), val('brename'), eval('brename --help | head -n 4 | tail -1 | sed "s/Version: //"'), topic: versions

    script:
    def sampID      = "${meta.id}"
    def chunkPrefix = "${meta.id}_chunk${meta.chunk_id}"
    def itsx_heuristics = params.ITSx_heuristics
        ? '--heuristics T'
        : '--heuristics F'

    // Allow inclusion of sequences that only find a single domain, given that they meet the given E-value and score thresholds, on with parameters 1e-9,0 by default
    // singledomain = params.ITSx_singledomain ? "--allow_single_domain 1e-9,0" : ""

    """
    echo -e "Extraction of rRNA regions using ITSx\\n"
    echo -e "Input sample: " ${sampID}
    echo -e "Chunk ID: "     ${meta.chunk_id}

    ## ITSx cannot read gz-compressed input
    ## Check if the input file is gz-compressed (by magic bytes `1f 8b`)
    tmpfile=""
    if head -c 2 -- "${input}" | LC_ALL=C od -An -tx1 | tr -d ' \n' | grep -qi '^1f8b'; then
      echo -e "Input file is gz-compressed, decompressing..."
      tmpfile="\$(mktemp "tmp.decompressed.input.XXXXXX")"
      gunzip -c -- "${input}" > "\$tmpfile"
      itsxinput="\$tmpfile"
    else
      itsxinput="${input}"
    fi

    ## ITSx extraction
    echo -e "\\nITSx extraction"
    ITSx \
      -i "\$itsxinput" \
      --complement ${params.ITSx_complement} \
      --save_regions all \
      --graphical F \
      --detailed_results T \
      --positions T \
      --not_found T \
      -E ${params.ITSx_evalue} \
      -t ${params.ITSx_tax} \
      ${itsx_heuristics} \
      --partial ${params.ITSx_partial} \
      --cpu ${task.cpus} \
      --preserve T \
      -o "${chunkPrefix}"

    echo -e "..Done"

    ## Remove the decompressed copy of the input as soon as ITSx is finished
    if [ -n "\$tmpfile" ]; then rm -f -- "\$tmpfile"; fi

      # ITSx.full.fasta
      # ITSx.SSU.fasta
      # ITSx.ITS1.fasta
      # ITSx.5_8S.fasta
      # ITSx.ITS2.fasta
      # ITSx.LSU.fasta
      # ITSx.positions.txt
      # ITSx.problematic.txt
      # ITSx_no_detections.fasta
      # ITSx_no_detections.txt
      # ITSx.summary.txt
      # ITSx.extraction.results
      # ITSx.SSU.full_and_partial.fasta
      # ITSx.ITS1.full_and_partial.fasta
      # ITSx.5_8S.full_and_partial.fasta
      # ITSx.ITS2.full_and_partial.fasta
      # ITSx.LSU.full_and_partial.fasta

    ## If partial sequences were required, remove empty sequences
    if [ \$(find . -type f -name "*.full_and_partial.fasta" | wc -l) -gt 0 ]; then
      echo -e "\\nPartial files found, removing empty sequences"

      find . -name "*.full_and_partial.fasta" \
        | parallel -j${task.cpus} "seqkit seq -m 1 -w 0 {} > {.}_tmp.fasta"

      rm *.full_and_partial.fasta
      brename -p "_tmp" -r "" -f "_tmp.fasta\$"

    fi

    ## Remove empty files (no sequences)
    echo -e "\\nRemoving empty files"
    find . -type f -name "*.fasta" -empty -print -delete
    echo -e "..Done"

    ## Compress results
    echo -e "\\nCompressing files"

    find . -type f -name "${chunkPrefix}*.fasta" \
      | parallel -j${task.cpus} "gzip -${params.gzip_compression} {}"

    gzip -${params.gzip_compression} "${chunkPrefix}".extraction.results

    echo -e "..Done"
    """
}



// ITSx processing workflow
workflow ITSx {


// Extract rRNA regions with ITSx2 from a single chunk of dereplicated sequences
// NB. ITSx2 delimits the cistron pan-eukaryotically with covariance models, therefore
//     the ITSx v1.x options `-t`, `-E`, `--partial`, `--complement` and `--heuristics`
//     have no counterpart here and are deliberately not passed
process itsx2 {

    label "main_container"

    // Chunk-level results are not published - they are concatenated per sample first
    // cpus 3

    tag "${meta.id}__chunk${meta.chunk_id}"

    input:
      tuple val(meta), path(input)   // FASTA file with dereplicated sequences (may be gz-compressed)

    output:
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}.full.fasta.gz"), emit: itsx_full, optional: true
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}.SSU.fasta.gz"),  emit: itsx_ssu,  optional: true
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}.ITS1.fasta.gz"), emit: itsx_its1, optional: true
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}.5_8S.fasta.gz"), emit: itsx_58s,  optional: true
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}.ITS2.fasta.gz"), emit: itsx_its2, optional: true
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}.LSU.fasta.gz"),  emit: itsx_lsu,  optional: true
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}.positions.txt"),       emit: itsx_positions,  optional: true
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}_no_detections.txt"),   emit: itsx_nondetects_txt, optional: true
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}.summary.txt"),         emit: itsx_summary,    optional: true
      tuple val(meta), path("${meta.id}_chunk${meta.chunk_id}.jsonl.gz"),            emit: itsx_jsonl,      optional: true
      tuple val("${task.process}"), val('ITSx2'), eval('itsx2 itsx --version | sed "s/itsx2 itsx //"'), topic: versions
      tuple val("${task.process}"), val('Infernal'), eval('cmsearch -h 2>&1 | sed -n "s/^# INFERNAL \\([0-9][0-9a-z.]*\\).*/\\1/p" | head -n 1'), topic: versions
      tuple val("${task.process}"), val('parallel'), eval('parallel --version | head -n 1 | sed "s/GNU parallel //"'), topic: versions

    script:
    def sampID      = "${meta.id}"
    def chunkPrefix = "${meta.id}_chunk${meta.chunk_id}"
    """
    echo -e "Extraction of rRNA regions using ITSx2\\n"
    echo -e "Input sample: " ${sampID}
    echo -e "Chunk ID: "     ${meta.chunk_id}

    ## ITSx2 reads gz-compressed FASTA/FASTQ natively - no decompression needed
    echo -e "\\nITSx2 extraction"
    itsx2 itsx \
      -i "${input}" \
      -o "${chunkPrefix}" \
      --save_regions all \
      --fasta_out \
      --cpu ${task.cpus}

    echo -e "..Done"

      # ITSx2.full.fasta
      # ITSx2.SSU.fasta
      # ITSx2.ITS1.fasta
      # ITSx2.5_8S.fasta
      # ITSx2.ITS2.fasta
      # ITSx2.LSU.fasta
      # ITSx2.positions.txt
      # ITSx2_no_detections.txt
      # ITSx2.summary.txt
      # ITSx2.jsonl

    ## Remove empty files (no sequences)
    echo -e "\\nRemoving empty files"
    find . -type f -name "*.fasta" -empty -print -delete
    echo -e "..Done"

    ## Compress results
    ## NB. `find` is used to avoid gzipping the symlinked input chunk
    echo -e "\\nCompressing files"

    find . -type f -name "${chunkPrefix}*.fasta" \
      | parallel -j${task.cpus} "gzip -${params.gzip_compression} {}"

    if [ -f "${chunkPrefix}".jsonl ]; then
      gzip -${params.gzip_compression} "${chunkPrefix}".jsonl
    fi

    echo -e "..Done"
    """
}


// Concatenate the extractor output from all chunks (per sample)
// + convert the extracted rRNA regions to Parquet
process itsx_concatenate {

    label "main_container"

    publishDir "${params.outdir}/03_ITSx", mode: "${params.storagemode}"
    // cpus 1

    tag "${meta.id}"

    input:
      tuple val(meta), path(chunk_files, stageAs: "chunks/")  // all files produced by the extractor, for all chunks of a sample

    output:
      tuple val(meta), path("${meta.id}.full.fasta.gz"), emit: itsx_full, optional: true
      tuple val(meta), path("${meta.id}.SSU.fasta.gz"),  emit: itsx_ssu,  optional: true
      tuple val(meta), path("${meta.id}.ITS1.fasta.gz"), emit: itsx_its1, optional: true
      tuple val(meta), path("${meta.id}.5_8S.fasta.gz"), emit: itsx_58s,  optional: true
      tuple val(meta), path("${meta.id}.ITS2.fasta.gz"), emit: itsx_its2, optional: true
      tuple val(meta), path("${meta.id}.LSU.fasta.gz"),  emit: itsx_lsu,  optional: true
      tuple val(meta), path("${meta.id}.positions.txt"),   emit: itsx_positions,   optional: true
      tuple val(meta), path("${meta.id}.problematic.txt"), emit: itsx_problematic, optional: true
      tuple val(meta), path("${meta.id}_no_detections.fasta.gz"), emit: itsx_nondetects,     optional: true
      tuple val(meta), path("${meta.id}_no_detections.txt"),      emit: itsx_nondetects_txt, optional: true
      tuple val(meta), path("${meta.id}.summary.txt"),            emit: itsx_summary, optional: true
      tuple val(meta), path("${meta.id}.extraction.results.gz"),  emit: itsx_details, optional: true
      tuple val(meta), path("${meta.id}.jsonl.gz"),               emit: itsx_jsonl,   optional: true
      tuple val(meta), path("${meta.id}.SSU.full_and_partial.fasta.gz"),  emit: itsx_ssu_part,  optional: true
      tuple val(meta), path("${meta.id}.ITS1.full_and_partial.fasta.gz"), emit: itsx_its1_part, optional: true
      tuple val(meta), path("${meta.id}.5_8S.full_and_partial.fasta.gz"), emit: itsx_58s_part,  optional: true
      tuple val(meta), path("${meta.id}.ITS2.full_and_partial.fasta.gz"), emit: itsx_its2_part, optional: true
      tuple val(meta), path("${meta.id}.LSU.full_and_partial.fasta.gz"),  emit: itsx_lsu_part,  optional: true
      tuple val(meta), path("parquet/*.parquet"), emit: parquet, optional: true
      tuple val("${task.process}"), val('duckdb'), eval('duckdb --version | cut -d" " -f1 | sed "s/^v//"'), topic: versions
      tuple val("${task.process}"), val('seqkit'), eval('seqkit version | sed "s/seqkit v//"'), topic: versions

    script:
    def sampID = "${meta.id}"
    """
    echo -e "Concatenating the extractor output from all chunks"
    echo -e "Input sample: " ${sampID}

    shopt -s nullglob

    ## Concatenate chunk files matching a suffix into a single per-sample file
    ## (a no-op if no chunk produced that file - e.g. ITSx2 has no `problematic` output)
    concat_chunks () {
      local suffix="\$1"     # chunk file suffix, e.g. ".ITS1.fasta.gz"
      local outfile="\$2"    # per-sample output file
      local label="\$3"      # human-readable label for the log

      local files=( chunks/${sampID}_chunk*"\$suffix" )
      echo -e "  - \$label: \${#files[@]}"
      if [ \${#files[@]} -gt 0 ]; then
        for f in "\${files[@]}"; do echo -e "        \$f"; done
        cat "\${files[@]}" > "\$outfile"
      fi
    }

    echo -e "Concatenating:"

    ## rRNA regions (full-length detections)
    concat_chunks ".full.fasta.gz" "${sampID}.full.fasta.gz" "full ITS sequences"
    concat_chunks ".SSU.fasta.gz"  "${sampID}.SSU.fasta.gz"  "SSU sequences"
    concat_chunks ".ITS1.fasta.gz" "${sampID}.ITS1.fasta.gz" "ITS1 sequences"
    concat_chunks ".5_8S.fasta.gz" "${sampID}.5_8S.fasta.gz" "5.8S sequences"
    concat_chunks ".ITS2.fasta.gz" "${sampID}.ITS2.fasta.gz" "ITS2 sequences"
    concat_chunks ".LSU.fasta.gz"  "${sampID}.LSU.fasta.gz"  "LSU sequences"

    ## rRNA regions (full + partial detections; ITSx v1.x only)
    concat_chunks ".SSU.full_and_partial.fasta.gz"  "${sampID}.SSU.full_and_partial.fasta.gz"  "SSU partial sequences"
    concat_chunks ".ITS1.full_and_partial.fasta.gz" "${sampID}.ITS1.full_and_partial.fasta.gz" "ITS1 partial sequences"
    concat_chunks ".5_8S.full_and_partial.fasta.gz" "${sampID}.5_8S.full_and_partial.fasta.gz" "5.8S partial sequences"
    concat_chunks ".ITS2.full_and_partial.fasta.gz" "${sampID}.ITS2.full_and_partial.fasta.gz" "ITS2 partial sequences"
    concat_chunks ".LSU.full_and_partial.fasta.gz"  "${sampID}.LSU.full_and_partial.fasta.gz"  "LSU partial sequences"

    ## Sequences with no rRNA detections
    concat_chunks "_no_detections.fasta.gz" "${sampID}_no_detections.fasta.gz" "no-detection sequences"
    concat_chunks "_no_detections.txt"      "${sampID}_no_detections.txt"      "no-detection IDs"

    ## Region coordinates and diagnostics
    concat_chunks ".positions.txt"          "${sampID}.positions.txt"          "positions"
    concat_chunks ".problematic.txt"        "${sampID}.problematic.txt"        "problematic sequences"
    concat_chunks ".extraction.results.gz"  "${sampID}.extraction.results.gz"  "extraction results"
    concat_chunks ".jsonl.gz"               "${sampID}.jsonl.gz"               "per-record calls (JSONL)"

    ## Summary reports cannot simply be concatenated - the per-chunk counts must be summed
    sum_files=( chunks/${sampID}_chunk*.summary.txt )
    echo -e "  - summary reports: \${#sum_files[@]}"
    if [ \${#sum_files[@]} -gt 0 ]; then
      for f in "\${sum_files[@]}"; do echo -e "        \$f"; done
      merge_itsx_summaries.sh -o ${sampID}.summary.txt "\${sum_files[@]}"
    fi

    echo -e "\\n"

    ## Convert the extracted regions to Parquet
    if [ ${params.ITSx_to_parquet} == true ]; then

      echo -e "\\nConverting the extracted regions to Parquet"
      mkdir -p parquet

      for region in full SSU ITS1 5_8S ITS2 LSU; do
        if [ -s "${sampID}.\${region}.fasta.gz" ]; then
          ITSx_to_DuckDB.sh \
            -i "${sampID}.\${region}.fasta.gz" \
            -o "parquet/${sampID}.\${region}.parquet"
        fi
      done

      echo -e "Parquet files created\\n"

    fi
    """
}


// Get near-full-length ITS from the extractor output (based on the positions file)
process get_its {

    label "main_container"

    publishDir "${params.outdir}/03_ITSx", mode: "${params.storagemode}"
    // cpus 1

    tag "${meta.id}"

    input:
      tuple val(meta), path(derep), path(positions)   // dereplicated primer-trimmed sequences + region coordinates

    output:
      tuple val(meta), path("${meta.id}_ITS1_58S_ITS2.fasta.gz"), emit: itsnf,  optional: true
      tuple val(meta), path("${meta.id}.extraction.tsv.gz"),      emit: report, optional: true
      tuple val("${task.process}"), val('R'), eval('Rscript -e "cat(R.version.string)" | sed "s/R version //" | cut -d" " -f1'), topic: versions
      tuple val("${task.process}"), val('Biostrings'), eval('Rscript -e "cat(as.character(packageVersion(\'Biostrings\')))"'), topic: versions

    script:
    def sampID = "${meta.id}"
    """
    echo -e "Extracting ITS1-5.8S-ITS2 region"
    echo -e "Input sample: " ${sampID}

    ## Run extraction (+ validation and exclusion of problematic sequences)
    ## NB. the same dereplicated sequences were used as the extractor input,
    ##     therefore the sequence IDs are guaranteed to match the positions file
    extract_itsx_regions.R \
      --fasta     ${derep} \
      --positions ${positions} \
      --region    ITS \
      --output    ${sampID}_ITS1_58S_ITS2.fasta.gz \
      --report    ${sampID}.extraction.tsv.gz
    """
}


  take:
    seqs

  main:

    // Add metadata to the channel (fetch sample ID from the FASTQ file name)
    ch_seqs = seqs.map { fastq ->
          def sample_id = fastq.getSimpleName().replaceAll(/_PrimerChecked/, '')
          def meta = [id: sample_id]
          [meta, fastq]
      }

    // Trim primers and dereplicate at sample level
    primer_trim(ch_seqs)

    // Size of dereplicated input for ITSx
    //   if null, use default value (currently, 10000)
    //   if 0, use all sequences in one chunk
    def chunk_size = (params.ITSx_chunk_size == null ? 10000 : params.ITSx_chunk_size as int)

    if( chunk_size == 0 ) {
      // Single-chunk workflow (no data splitting)
      // NB!  here, fasta will be gz-compressed -> will be handled in the itsx process
      chunks_ch = primer_trim.out.derep
        .map { meta, fasta -> [ meta + [chunk_id: null], fasta ] }
    }
    else {
      // Chunking mode: split the dereplicated primer-trimmed sequences (at sample level) into chunks while preserving metadata
      // NB!  here, fasta will be uncompressed
      chunks_ch = primer_trim.out.derep
        .flatMap { meta, fasta ->
          def chunks = fasta.splitFasta(by: chunk_size, file: true, decompress: true, compress: false)
          def result = []
          chunks.eachWithIndex { chunk_file, idx ->
            result << [ meta + [chunk_id: idx], chunk_file ]
          }
          return result
        }
    }

    // Run ITSx
    itsx(chunks_ch)

    // For single-chunk workflow, concatenate all chunks for each sample
    if(params.ITSx_chunk_size == 0){

        // Fetch results from the ITSx
        ch_res_itsx_full        = itsx.out.itsx_full
        ch_res_itsx_ssu         = itsx.out.itsx_ssu
        ch_res_itsx_its1        = itsx.out.itsx_its1
        ch_res_itsx_58s         = itsx.out.itsx_58s
        ch_res_itsx_its2        = itsx.out.itsx_its2
        ch_res_itsx_lsu         = itsx.out.itsx_lsu
        ch_res_itsx_positions   = itsx.out.itsx_positions
        ch_res_itsx_problematic = itsx.out.itsx_problematic
        ch_res_itsx_nondetects  = itsx.out.itsx_nondetects
        ch_res_itsx_summary     = itsx.out.itsx_summary
        ch_res_itsx_details     = itsx.out.itsx_details
        ch_res_itsx_ssu_part    = itsx.out.itsx_ssu_part
        ch_res_itsx_its1_part   = itsx.out.itsx_its1_part
        ch_res_itsx_58s_part    = itsx.out.itsx_58s_part
        ch_res_itsx_its2_part   = itsx.out.itsx_its2_part
        ch_res_itsx_lsu_part    = itsx.out.itsx_lsu_part
        
        if(params.ITSx_to_parquet == true ){
          itsx_to_parquet(
            itsx.out.itsx_full,
            itsx.out.itsx_ssu,
            itsx.out.itsx_its1,
            itsx.out.itsx_58s,
            itsx.out.itsx_its2,
            itsx.out.itsx_lsu
          )
          ch_res_parquet = itsx_to_parquet.out.parquet
        } else {
          ch_res_parquet = channel.empty()
        }

    } else {
    // For multi-chunk workflow, we need to pool the chunks per sample

      // Group all ITSx chunk outputs back by sample ID and concatenate
      itsx_all_chunks = itsx.out.itsx_full
          .mix(
            itsx.out.itsx_ssu,
            itsx.out.itsx_its1,
            itsx.out.itsx_58s,
            itsx.out.itsx_its2,
            itsx.out.itsx_lsu,
            itsx.out.itsx_nondetects,
            itsx.out.itsx_summary,
            itsx.out.itsx_details,
            itsx.out.itsx_positions,
            itsx.out.itsx_problematic,
            itsx.out.itsx_ssu_part,
            itsx.out.itsx_its1_part,
            itsx.out.itsx_58s_part,
            itsx.out.itsx_its2_part,
            itsx.out.itsx_lsu_part
          )

        concatenated_ch = itsx_all_chunks
          .map { meta, file ->
              [meta.id, meta, file]
          }
          .groupTuple(by: 0)
          .map { sample_id, metas, files ->
              [metas[0], files]
          }
    
        // Concatenate all chunks for each sample (no-op if channel is empty)
        itsx_concatenate(concatenated_ch)

        // Fetch results from the concatenated channel
        ch_res_itsx_full        = itsx_concatenate.out.itsx_full
        ch_res_itsx_ssu         = itsx_concatenate.out.itsx_ssu
        ch_res_itsx_its1        = itsx_concatenate.out.itsx_its1
        ch_res_itsx_58s         = itsx_concatenate.out.itsx_58s
        ch_res_itsx_its2        = itsx_concatenate.out.itsx_its2
        ch_res_itsx_lsu         = itsx_concatenate.out.itsx_lsu
        ch_res_itsx_positions   = itsx_concatenate.out.itsx_positions
        ch_res_itsx_problematic = itsx_concatenate.out.itsx_problematic
        ch_res_itsx_nondetects  = itsx_concatenate.out.itsx_nondetects
        ch_res_itsx_summary     = itsx_concatenate.out.itsx_summary
        ch_res_itsx_details     = itsx_concatenate.out.itsx_details
        ch_res_itsx_ssu_part    = itsx_concatenate.out.itsx_ssu_part
        ch_res_itsx_its1_part   = itsx_concatenate.out.itsx_its1_part
        ch_res_itsx_58s_part    = itsx_concatenate.out.itsx_58s_part
        ch_res_itsx_its2_part   = itsx_concatenate.out.itsx_its2_part
        ch_res_itsx_lsu_part    = itsx_concatenate.out.itsx_lsu_part
        ch_res_parquet          = itsx_concatenate.out.parquet
    }




  // Collect ITSx-extracted sequences
  if(params.its_region == "full" || params.its_region == "ITS1" || params.its_region == "ITS2" || params.its_region == "SSU" || params.its_region == "LSU" || params.its_region == "ITS1_5.8S_ITS2"){

    // Collect rRNA parts into separate channels (+ drop metadata)
    ch_cc_full = itsx.out.itsx_full.map { meta, fasta -> fasta }.flatten().collect().ifEmpty(file("NOFULL"))
    ch_cc_ssu  = itsx.out.itsx_ssu.map { meta, fasta -> fasta }.flatten().collect().ifEmpty(file("NOSSU"))
    ch_cc_its1 = itsx.out.itsx_its1.map { meta, fasta -> fasta }.flatten().collect().ifEmpty(file("NOITS1"))
    ch_cc_58s  = itsx.out.itsx_58s.map { meta, fasta -> fasta }.flatten().collect().ifEmpty(file("NO58S"))
    ch_cc_its2 = itsx.out.itsx_its2.map { meta, fasta -> fasta }.flatten().collect().ifEmpty(file("NOITS2"))
    ch_cc_lsu  = itsx.out.itsx_lsu.map { meta, fasta -> fasta }.flatten().collect().ifEmpty(file("NOLSU"))
    
    ch_cc_ssu_part  = itsx.out.itsx_ssu_part.map { meta, fasta -> fasta }.flatten().collect().ifEmpty(file("NOSSUPART"))
    ch_cc_its1_part = itsx.out.itsx_its1_part.map { meta, fasta -> fasta }.flatten().collect().ifEmpty(file("NOITS1PART"))
    ch_cc_58s_part  = itsx.out.itsx_58s_part.map { meta, fasta -> fasta }.flatten().collect().ifEmpty(file("NO58SPART"))
    ch_cc_its2_part = itsx.out.itsx_its2_part.map { meta, fasta -> fasta }.flatten().collect().ifEmpty(file("NOITS2PART"))
    ch_cc_lsu_part  = itsx.out.itsx_lsu_part.map { meta, fasta -> fasta }.flatten().collect().ifEmpty(file("NOLSUPART"))

    itsx_collect(
      ch_cc_full,
      ch_cc_ssu,
      ch_cc_its1,
      ch_cc_58s,
      ch_cc_its2,
      ch_cc_lsu,
      ch_cc_ssu_part,
      ch_cc_its1_part,
      ch_cc_58s_part,
      ch_cc_its2_part,
      ch_cc_lsu_part
      )

  } // end of collection of ITSx-extracted sequences


  emit:
    hashes           = primer_trim.out.hashes.map { meta, file -> file }
    uc               = primer_trim.out.uc.map { meta, file -> file }
    trimmed_seqs     = primer_trim.out.trimmed_seqs.map { meta, file -> file }
    itsx_full        = ch_res_itsx_full
    itsx_ssu         = ch_res_itsx_ssu
    itsx_its1        = ch_res_itsx_its1
    itsx_58s         = ch_res_itsx_58s
    itsx_its2        = ch_res_itsx_its2
    itsx_lsu         = ch_res_itsx_lsu
    itsx_positions   = ch_res_itsx_positions
    itsx_problematic = ch_res_itsx_problematic
    itsx_nondetects  = ch_res_itsx_nondetects
    itsx_summary     = ch_res_itsx_summary
    itsx_details     = ch_res_itsx_details
    itsx_ssu_part    = ch_res_itsx_ssu_part
    itsx_its1_part   = ch_res_itsx_its1_part
    itsx_58s_part    = ch_res_itsx_58s_part
    itsx_its2_part   = ch_res_itsx_its2_part
    itsx_lsu_part    = ch_res_itsx_lsu_part
    parquet          = ch_res_parquet


} // end of ITSx workflow
