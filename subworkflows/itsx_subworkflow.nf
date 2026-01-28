/*
============================================================================
  NextITS: Pipeline to process eukaryotic ITS amplicons
============================================================================
  License: Apache-2.0
  Github : https://github.com/vmikk/NextITS
  Website: https://Next-ITS.github.io/
----------------------------------------------------------------------------
*/

// Subworkflow for ITSx processing
// (which splits large dereplicated FASTAs into chunks).
// The workflow is as follows:
//   1. Trim primers and dereplicate at sample level
//   2. Split the dereplicated primer-trimmed sequences (at sample level) into chunks while preserving metadata
//   3. Run ITSx on each chunk
//   4. Group results back by original sample ID and concatenate + convert ITSx output to Parquet

// Path to the output results
out_3_itsx   = params.outdir + "/03_ITSx"
out_3_itsxp  = params.outdir + "/03_ITSx_PooledParts"


// Trim primers and dereplicate at sample level
process primer_trim {

    label "main_container"

    publishDir "${out_3_itsx}", mode: "${params.storagemode}"
    // cpus 2

    // Add sample ID to the log file
    tag "${meta.id}"

    input:
      tuple val(meta), path(fastq)

    output:
      tuple val(meta), path("${meta.id}_derep.fasta.gz"),    emit: derep,  optional: true
      tuple val(meta), path("${meta.id}_hash_table.txt.gz"), emit: hashes, optional: true
      tuple val(meta), path("${meta.id}_uc.uc.gz"),          emit: uc,     optional: true
      tuple val(meta), path("${meta.id}_primertrimmed_sorted.fq.gz"), emit: trimmed_seqs,   optional: true
      tuple val("${task.process}"), val('cutadapt'), eval('cutadapt --version'), topic: versions
      tuple val("${task.process}"), val('vsearch'), eval('vsearch --version 2>&1 | head -n 1 | sed "s/vsearch //g" | sed "s/,.*//g" | sed "s/^v//" | sed "s/_.*//"'), topic: versions
      tuple val("${task.process}"), val('seqkit'), eval('seqkit version | sed "s/seqkit v//"'), topic: versions
      tuple val("${task.process}"), val('phredsort'), eval('phredsort -v | sed "s/phredsort //"'), topic: versions
      tuple val("${task.process}"), val('seqhasher'), eval('seqhasher -v | sed "s/SeqHasher //"'), topic: versions
      tuple val("${task.process}"), val('parallel'), eval('parallel --version | head -n 1 | sed "s/GNU parallel //"'), topic: versions
      tuple val("${task.process}"), val('brename'), eval('brename --help | head -n 4 | tail -1 | sed "s/Version: //"'), topic: versions

    script:
    sampID="${meta.id}"
    """
    echo -e "Primer trimming and dereplication at sample level\\n"
    echo -e "Input sample: " ${sampID}

    ## Trim primers    
    echo -e "Trimming primers\\n"

    ## Reverse-complement rev primer
    RR=\$(rc.sh ${params.primer_reverse})

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
      exit 0
    fi
   
    ## Estimate sequence quality and sort sequences by quality
    echo -e "\\nSorting by sequence quality"
    seqkit replace -p "\\s.+" ${sampID}_primertrimmed.fq.gz \
      | phredsort -i - -o - --metric meep --header avgphred,maxee,meep \
      | gzip -1 > ${sampID}_primertrimmed_sorted.fq.gz
    echo -e "..Done"

    ## Hash sequences, add sample ID to the header
    ## columns: Sample ID - Hash - PacBioID - AvgPhredScore - MaxEE - MEEP - Sequence - Quality - Length
    ## Convert to Parquet format
    echo -e "\\nCreating hash table"
    seqhasher --hash sha1 --name ${sampID} ${sampID}_primertrimmed_sorted.fq.gz - \
      | seqkit fx2tab --length \
      | sed 's/;/\t/ ; s/;/\t/ ; s/ avgphred=/\t/ ; s/ maxee=/\t/ ; s/ meep=/\t/' \
      > ${sampID}_hash_table.txt
    echo -e "..Done"

    ## Check the number of fields per record (should be 9!)
    # awk '{print NF}' ${sampID}_hash_table.txt | sort | uniq -c
    # awk 'NF > 9 {print \$0 }' ${sampID}_hash_table.txt

    ## Dereplicate at sample level (use quality-sorted sequences to make sure that the representative sequence is with the highest quality)
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
      > ${sampID}_derep.fasta
    
    ## Remove temporary file
    rm ${sampID}_primertrimmed.fq.gz

    ## Compress results
    echo -e "\\nCompressing results"
    parallel -j${task.cpus} "gzip -${params.gzip_compression} {}" ::: \
      ${sampID}_hash_table.txt \
      ${sampID}_uc.uc \
      ${sampID}_derep.fasta

    echo -e "..Done"

    """
}


// Extract ITS region with ITSx
// NB. In input data, sequence header should not contain spaces!
process itsx {

    label "main_container"

    // No need to publish intermediate results for chunked workflow, as they will be concatenated later
    publishDir "${out_3_itsx}", 
               mode: "${params.storagemode}",
               enabled: params.ITSx_chunk_size == 0

    // cpus 2

    // Add sample ID to the log file
    tag { meta.chunk_id != null ? "${meta.id}__chunk${meta.chunk_id}" : "${meta.id}" }

    input:
      tuple val(meta), path(input)   // FASTA file with dereplicated sequences

    output:
      tuple val(meta), path( "${meta.id}*.full.fasta.gz"), emit: itsx_full, optional: true
      tuple val(meta), path( "${meta.id}*.SSU.fasta.gz"),  emit: itsx_ssu,  optional: true
      tuple val(meta), path( "${meta.id}*.ITS1.fasta.gz"), emit: itsx_its1, optional: true
      tuple val(meta), path( "${meta.id}*.5_8S.fasta.gz"), emit: itsx_58s,  optional: true
      tuple val(meta), path( "${meta.id}*.ITS2.fasta.gz"), emit: itsx_its2, optional: true
      tuple val(meta), path( "${meta.id}*.LSU.fasta.gz"),  emit: itsx_lsu,  optional: true
      tuple val(meta), path( "${meta.id}*.positions.txt"),   emit: itsx_positions, optional: true
      tuple val(meta), path( "${meta.id}*.problematic.txt"), emit: itsx_problematic, optional: true
      tuple val(meta), path( "${meta.id}*_no_detections.fasta.gz"), emit: itsx_nondetects, optional: true
      tuple val(meta), path( "${meta.id}*.summary.txt"),        emit: itsx_summary, optional: true
      tuple val(meta), path( "${meta.id}*.extraction.results.gz"), emit: itsx_details, optional: true
      tuple val(meta), path( "${meta.id}*.SSU.full_and_partial.fasta.gz"),  emit: itsx_ssu_part,  optional: true
      tuple val(meta), path( "${meta.id}*.ITS1.full_and_partial.fasta.gz"), emit: itsx_its1_part, optional: true
      tuple val(meta), path( "${meta.id}*.5_8S.full_and_partial.fasta.gz"), emit: itsx_58s_part,  optional: true
      tuple val(meta), path( "${meta.id}*.ITS2.full_and_partial.fasta.gz"), emit: itsx_its2_part, optional: true
      tuple val(meta), path( "${meta.id}*.LSU.full_and_partial.fasta.gz"),  emit: itsx_lsu_part,  optional: true
      tuple val("${task.process}"), val('ITSx'), eval('ITSx --help 2>&1 | head -n 3 | tail -n 1 | sed "s/Version: //"'), topic: versions
      tuple val("${task.process}"), val('cutadapt'), eval('cutadapt --version'), topic: versions
      tuple val("${task.process}"), val('vsearch'), eval('vsearch --version 2>&1 | head -n 1 | sed "s/vsearch //g" | sed "s/,.*//g" | sed "s/^v//" | sed "s/_.*//"'), topic: versions
      tuple val("${task.process}"), val('seqkit'), eval('seqkit version | sed "s/seqkit v//"'), topic: versions
      tuple val("${task.process}"), val('phredsort'), eval('phredsort -v | sed "s/phredsort //"'), topic: versions
      tuple val("${task.process}"), val('seqhasher'), eval('seqhasher -v | sed "s/SeqHasher //"'), topic: versions
      tuple val("${task.process}"), val('parallel'), eval('parallel --version | head -n 1 | sed "s/GNU parallel //"'), topic: versions
      tuple val("${task.process}"), val('brename'), eval('brename --help | head -n 4 | tail -1 | sed "s/Version: //"'), topic: versions
      tuple val("${task.process}"), val('duckdb'), eval('duckdb --version | cut -d" " -f1 | sed "s/^v//"'), topic: versions

    script:
    sampID="${meta.id}"
    chunkPrefix="${meta.id}_chunk${meta.chunk_id}"

    // Allow inclusion of sequences that only find a single domain, given that they meet the given E-value and score thresholds, on with parameters 1e-9,0 by default
    // singledomain = params.ITSx_singledomain ? "--allow_single_domain 1e-9,0" : ""

    """
    echo -e "Extraction of rRNA regions using ITSx\\n"
    echo -e "Input sample: " ${sampID}
    echo -e "Chunk ID: " ${meta.chunk_id}
    
    ## Check if input file is gz-compressed (by magic bytes `1f 8b`)
    tmp_created=0
    tmpfile=""
    if [[ -f "${input}" ]] && head -c 2 -- "${input}" | LC_ALL=C od -An -tx1 | tr -d ' \n' | grep -qi '^1f8b'; then
      echo -e "Input file is gz-compressed, decompressing..."
      tmpfile="\$(mktemp "tmp.decompressed.input.XXXXXX")"
      gunzip -c -- "${input}" > "\$tmpfile"
      tmp_created=1
      itsxinput="\$tmpfile"
      itsxoutput="${sampID}"
    else
      itsxinput="${input}"
      itsxoutput="${chunkPrefix}"
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
      --partial ${params.ITSx_partial} \
      --cpu ${task.cpus} \
      --preserve T \
      -o "\$itsxoutput"
    
    echo -e "..Done"

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
      echo -e "Partial files found, removing empty sequences\\n."

      find . -name "*.full_and_partial.fasta" \
        | parallel -j${task.cpus} "seqkit seq -m 1 -w 0 {} > {.}_tmp.fasta"

      rm *.full_and_partial.fasta
      brename -p "_tmp" -r "" -f "_tmp.fasta\$"

    fi

    ## Remove empty files (no sequences)
    echo -e "\\nRemoving empty files"
    find . -type f -name "*.fasta" -empty -print -delete
    echo -e "..Done"

    ## Remove temporary file (if input file was gz-compressed)
    if (( tmp_created )); then
      rm -f -- "\$tmpfile"
    fi

    ## Compress results
    echo -e "\\nCompressing files"

    ## ITSx results (no symlinked derep input)
    find . -type f -name "*.fasta" \
    | parallel -j${task.cpus} "gzip -${params.gzip_compression} {}"

    gzip -${params.gzip_compression} "\$itsxoutput".extraction.results

    echo -e "..Done"
    """
}


