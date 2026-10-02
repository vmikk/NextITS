
// Quality-score check
// Are Phred scores binned (e.g., NovaSeq, NextSeq) or continuous (e.g., MiSeq)?
// Is there an excess of 3' poly-G tails (two-colour chemistry)?
process illumina_qcheck {

    label "main_container"

    publishDir "${params.outdir}/01_Demux", mode: "${params.storagemode}"
    // cpus 1

    input:
      tuple path(input_R1), path(input_R2)

    output:
      path "Quality_check.tsv", emit: tsv
      tuple val("${task.process}"), val('seqkit'), eval('seqkit version | sed "s/seqkit v//"'), topic: versions

    script:
    """
    echo -e "Checking quality-score encoding\\n"
    echo -e "Input R1: " ${input_R1}
    echo -e "Input R2: " ${input_R2}

    ## Number of reads to inspect (per mate)
    NREADS=100000

    printf "Mate\\tNumReads\\tMaxReadLength\\tNumDistinctQ\\tQValues\\tMaxQ\\tPolyG_Percent\\tQualityType\\n" \\
      > Quality_check.tsv

    check_mate () {
      # \$1 = mate label, \$2 = FASTQ file

      ## Distinct Phred scores (offset 33) and read lengths
      seqkit head -n \$NREADS "\$2" \\
        | seqkit seq --qual \\
        | awk -v mate="\$1" '
          BEGIN { for(i = 33; i < 127; i++) ord[sprintf("%c", i)] = i - 33 }
          {
            n++
            if(length(\$0) > maxlen) maxlen = length(\$0)
            for(i = 1; i <= length(\$0); i++) q[ ord[substr(\$0, i, 1)] ] = 1
          }
          END {
            nq = 0; maxq = -1; vals = ""
            for(k = 0; k <= 93; k++) if(k in q){ nq++; maxq = k; vals = vals (vals == "" ? "" : ",") k }
            printf "%s\\t%d\\t%d\\t%d\\t%s\\t%d", mate, n, maxlen, nq, vals, maxq
          }' > tmp_quals.txt

      ## Reads with a 3-prime poly-G tail (10 or more G's)
      NPG=\$(seqkit head -n \$NREADS "\$2" | seqkit seq --seq | { grep -c -E 'G{10}\$' || true; })
      NR=\$(cut -f2 tmp_quals.txt)
      NDQ=\$(cut -f4 tmp_quals.txt)

      awk -v npg="\$NPG" -v nr="\$NR" -v ndq="\$NDQ" 'BEGIN{
        printf "\\t%.2f\\t%s\\n", (nr > 0 ? 100 * npg / nr : 0), (ndq <= 8 ? "binned" : "continuous") }' \\
        > tmp_pg.txt

      paste -d '' tmp_quals.txt tmp_pg.txt >> Quality_check.tsv
      rm tmp_quals.txt tmp_pg.txt
    }

    check_mate "R1" ${input_R1}
    check_mate "R2" ${input_R2}

    echo -e "\\nQuality-score summary:"
    cat Quality_check.tsv
    """
}



// Demultiplexing of Illumina paired-end reads with cutadapt
//
// Tag layouts (`params.illumina_barcodetype`):
//   "dual_symmetric"  = the same tag at the 5' end of both mates
//   "dual_asymmetric" = different tags at the 5' ends of the mates (`fwd...rev` FASTA; symmetric entries allowed)
//   "single"          = tag at the 5' end of one mate only
//
// For dual tags:
//   Pass 1  strict demultiplexing - both mates must carry the tags of the same sample (`--pair-adapters`)
//           (+ pass 1b for asymmetric tags, with swapped mates - amplicons are in mixed orientation)
//   Pass 2  pairs left over, in which both mates carry a known tag, are tag-jumps (or unknown tag combinations)
//   Pass 3  rescue pairs with a single readable tag (if the tag is used by a single sample only)
//
// NB. `--revcomp` can not be combined with `--pair-adapters`
process demux_illumina {

    label "main_container"

    publishDir "${params.outdir}/01_Demux", mode: "${params.storagemode}"
    // cpus 8

    input:
      tuple path(input_R1), path(input_R2)
      path barcodes                       // validated tags (single or symmetric dual tags)
      path(tags_dual, stageAs: "dual/*")  // `tags_fwd.fasta` + `tags_rev.fasta` (dual tags) or a dummy file

    output:
      path "Demux/*.fq.gz",     emit: samples_demux, optional: true   // `{sample}_R1.fq.gz` and `{sample}_R2.fq.gz`
      path "Demux_summary.tsv", emit: summary
      path "Demux_totals.tsv",  emit: totals
      path "logs/*",            emit: logs
      tuple val("${task.process}"), val('cutadapt'), eval('cutadapt --version'), topic: versions
      tuple val("${task.process}"), val('seqkit'), eval('seqkit version | sed "s/seqkit v//"'), topic: versions

    script:
    def rescue = params.illumina_demux_rescue ? "true" : "false"
    """
    echo -e "Demultiplexing Illumina paired-end reads with cutadapt\\n"
    echo -e "Input R1:     " ${input_R1}
    echo -e "Input R2:     " ${input_R2}
    echo -e "Barcodes:     " ${barcodes}
    echo -e "Tag layout:   " ${params.illumina_barcodetype}
    echo -e "Tag errors:   " ${params.barcode_errors}
    echo -e "Tag window:   " ${params.barcode_window}
    echo -e "Rescue:       " ${rescue}

    ## cutadapt keeps one file open per sample and mate
    ulimit -S -n 4096 2>/dev/null || true

    mkdir -p logs json Strict StrictB Rescued Demux

    ## Shared options for tag matching
    TAGOPTS="--errors ${params.barcode_errors} --no-indels --overlap ${params.barcode_overlap} --cores ${task.cpus}"

    ## Number of pairs with a match, from the cutadapt JSON report
    json_count () {
      python3 -c 'import json, sys; print(json.load(open(sys.argv[1]))["read_counts"]["read1_with_adapter"] or 0)' "\$1"
    }


    ## ---------------------------------------------------------------------------
    ## Tags
    ## ---------------------------------------------------------------------------

    HAS_DUAL=false
    if [ -s dual/tags_fwd.fasta ] && [ -s dual/tags_rev.fasta ]; then HAS_DUAL=true; fi

    case "${params.illumina_barcodetype}" in

      "single"|"dual_symmetric")
        if [[ \$HAS_DUAL == true ]]; then
          echo -e "\\nERROR: dual tags (in 'fwd...rev' format) were provided, but '--illumina_barcodetype ${params.illumina_barcodetype}' expects a single tag per sample."
          echo -e "Use '--illumina_barcodetype dual_asymmetric', or provide one tag per sample.\\n"
          exit 1
        fi

        prepare_illumina_tags.py \\
          --mode   ${params.illumina_barcodetype} \\
          --tags   ${barcodes} \\
          --window ${params.barcode_window} \\
          --outdir tags
        REVTAG="NA"
        ;;

      "dual_asymmetric")
        if [[ \$HAS_DUAL == false ]]; then
          echo -e "\\nERROR: '--illumina_barcodetype dual_asymmetric' requires tags in 'fwd...rev' format."
          echo -e "For the same tag on both mates use '--illumina_barcodetype dual_symmetric'.\\n"
          exit 1
        fi

        for ORI in forward revcomp; do
          prepare_illumina_tags.py \\
            --mode       dual_asymmetric \\
            --fwd        dual/tags_fwd.fasta \\
            --rev        dual/tags_rev.fasta \\
            --rev-orient \$ORI \\
            --window     ${params.barcode_window} \\
            --outdir     tags_\$ORI
        done

        ## Orientation of reverse tags
        ## (both mate orientations, as amplicons may be in mixed orientation)
        REVTAG="${params.illumina_revtag_orient}"
        if [[ \$REVTAG == "auto" ]]; then
          echo -e "\\nDetecting the orientation of reverse tags (first 20000 read pairs)"
          seqkit head -n 20000 ${input_R1} -o sub_R1.fq.gz
          seqkit head -n 20000 ${input_R2} -o sub_R2.fq.gz

          declare -A NPAIRS
          for ORI in forward revcomp; do
            N=0
            for SIDES in "T1 T2" "T2 T1"; do
              set -- \$SIDES
              cutadapt --pair-adapters --action=none \$TAGOPTS \\
                -g file:tags_\$ORI/\$1.fasta \\
                -G file:tags_\$ORI/\$2.fasta \\
                --json tmp.json \\
                -o /dev/null -p /dev/null \\
                sub_R1.fq.gz sub_R2.fq.gz > /dev/null
              N=\$(( N + \$(json_count tmp.json) ))
            done
            NPAIRS[\$ORI]=\$N
            echo -e "..reverse tags in \$ORI orientation: \$N pairs assigned"
          done
          rm -f sub_R1.fq.gz sub_R2.fq.gz tmp.json

          if (( NPAIRS[revcomp] > NPAIRS[forward] )); then REVTAG="revcomp"; else REVTAG="forward"; fi
          if (( NPAIRS[revcomp] == 0 && NPAIRS[forward] == 0 )); then
            echo -e "WARNING: no read pairs could be assigned in the test subset, using reverse tags as provided"
          fi
        fi
        echo -e "..Reverse tags are used in '\$REVTAG' orientation"
        mv tags_\$REVTAG tags
        rm -rf tags_forward tags_revcomp
        ;;

      *)
        echo -e "\\nERROR: unknown tag layout '${params.illumina_barcodetype}'\\n"
        exit 1
        ;;
    esac

    echo -e "\\nNumber of samples: " \$(grep -c '^>' tags/T1.fasta)


    ## ---------------------------------------------------------------------------
    ## Demultiplexing
    ## ---------------------------------------------------------------------------

    PASS3_JSON=""

    if [[ "${params.illumina_barcodetype}" == "single" ]]; then

      echo -e "\\n.. Single-tag demultiplexing (tag on either mate)\\n"
      cutadapt \\
        --revcomp --rename='{header}' \\
        \$TAGOPTS \\
        --minimum-length ${params.barcode_minlen} \\
        --discard-untrimmed \\
        -g file:tags/T1.fasta \\
        --json json/pass1.json \\
        -o "Strict/{name}_R1.fq.gz" \\
        -p "Strict/{name}_R2.fq.gz" \\
        ${input_R1} ${input_R2} \\
        > logs/cutadapt_pass1.log

      PASS1_JSONS="json/pass1.json"
      PASS2_JSON=""

    else

      ## Pass 1: strict demultiplexing, both mates must carry the tags of the same sample
      ## Non-matching pairs are left untrimmed (tags intact) and are used as input for the next passes
      echo -e "\\n.. Pass 1: strict demultiplexing (--pair-adapters)\\n"
      cutadapt \\
        --pair-adapters \\
        \$TAGOPTS \\
        --minimum-length ${params.barcode_minlen} \\
        -g file:tags/T1.fasta \\
        -G file:tags/T2.fasta \\
        --json json/pass1.json \\
        --untrimmed-output        unp1_R1.fq.gz \\
        --untrimmed-paired-output unp1_R2.fq.gz \\
        -o "Strict/{name}_R1.fq.gz" \\
        -p "Strict/{name}_R2.fq.gz" \\
        ${input_R1} ${input_R2} \\
        > logs/cutadapt_pass1.log

      PASS1_JSONS="json/pass1.json"

      ## Pass 1b: asymmetric tags in the opposite mate orientation
      if [[ "${params.illumina_barcodetype}" == "dual_asymmetric" ]]; then
        echo -e "\\n.. Pass 1b: strict demultiplexing, swapped mates\\n"
        cutadapt \\
          --pair-adapters \\
          \$TAGOPTS \\
          --minimum-length ${params.barcode_minlen} \\
          -g file:tags/T2.fasta \\
          -G file:tags/T1.fasta \\
          --json json/pass1b.json \\
          --untrimmed-output        unp_R1.fq.gz \\
          --untrimmed-paired-output unp_R2.fq.gz \\
          -o "StrictB/{name}_R1.fq.gz" \\
          -p "StrictB/{name}_R2.fq.gz" \\
          unp1_R1.fq.gz unp1_R2.fq.gz \\
          > logs/cutadapt_pass1b.log
        rm unp1_R1.fq.gz unp1_R2.fq.gz
        PASS1_JSONS="json/pass1.json json/pass1b.json"
      else
        mv unp1_R1.fq.gz unp_R1.fq.gz
        mv unp1_R2.fq.gz unp_R2.fq.gz
      fi

      ## Pass 2: discard tag-jumped pairs (both mates carry a known tag, but not a valid combination)
      ## `--action=none` keeps the tags, `--discard-trimmed --pair-filter=both` removes pairs with a tag on both mates
      echo -e "\\n.. Pass 2: discarding tag-jumped pairs\\n"
      cutadapt \\
        --action=none \\
        --discard-trimmed \\
        --pair-filter=both \\
        \$TAGOPTS \\
        -g file:tags/Tall.fasta \\
        -G file:tags/Tall.fasta \\
        --json json/pass2.json \\
        -o rescuable_R1.fq.gz \\
        -p rescuable_R2.fq.gz \\
        unp_R1.fq.gz unp_R2.fq.gz \\
        > logs/cutadapt_pass2.log
      rm unp_R1.fq.gz unp_R2.fq.gz

      PASS2_JSON="json/pass2.json"

      ## Pass 3: rescue pairs with a single readable tag
      ##   `-g` only (no `-G`) - that makes `--revcomp` meaningful (mates are swapped if the tag is found on R2)
      ##   only tags used by a single sample are considered
      if [[ ${rescue} == true ]] && [ -s tags/Tresc.fasta ]; then
        echo -e "\\n.. Pass 3: rescuing pairs with one readable tag\\n"
        cutadapt \\
          --revcomp --rename='{header}' \\
          \$TAGOPTS \\
          --minimum-length ${params.barcode_minlen} \\
          --discard-untrimmed \\
          -g file:tags/Tresc.fasta \\
          --json json/pass3.json \\
          -o "Rescued/{name}_R1.fq.gz" \\
          -p "Rescued/{name}_R2.fq.gz" \\
          rescuable_R1.fq.gz rescuable_R2.fq.gz \\
          > logs/cutadapt_pass3.log
        PASS3_JSON="json/pass3.json"
      else
        echo -e "\\n.. Pass 3 (rescue) skipped"
      fi
      rm rescuable_R1.fq.gz rescuable_R2.fq.gz

    fi


    ## ---------------------------------------------------------------------------
    ## Combine strict and rescued pairs per sample
    ## ---------------------------------------------------------------------------

    ## Remove empty outputs
    ##   an empty .gz is 20 bytes via igzip/pigz, but 37 bytes via Python's gzip,
    ##   while a single read is >= 250 bytes - so 100 is a safe threshold for both
    find Strict StrictB Rescued -type f -name "*.fq.gz" -size -100c -delete

    echo -e "\\n.. Combining strict and rescued pairs"
    grep '^>' tags/T1.fasta | sed 's/^>//' | sort -u > samples.txt
    while read -r s; do
      for r in R1 R2; do
        files=()
        for f in "Strict/\${s}_\${r}.fq.gz" "StrictB/\${s}_\${r}.fq.gz" "Rescued/\${s}~1_\${r}.fq.gz" "Rescued/\${s}~2_\${r}.fq.gz"; do
          if [ -s "\$f" ]; then files+=("\$f"); fi
        done
        if (( \${#files[@]} > 0 )); then
          cat "\${files[@]}" > "Demux/\${s}_\${r}.fq.gz"
        fi
      done
    done < samples.txt

    ## Per-sample read counts (R1 only = number of pairs)
    find Strict StrictB Rescued Demux -name "*_R1.fq.gz" | sort > r1_files.txt
    if [ -s r1_files.txt ]; then
      seqkit stats --tabular --threads ${task.cpus} --infile-list r1_files.txt > r1_counts.tsv
    else
      printf "file\\tformat\\ttype\\tnum_seqs\\n" > r1_counts.tsv
    fi

    ## Summary tables (fails if the combined counts differ from strict + rescued)
    summarize_demux_illumina.py \\
      --counts        r1_counts.tsv \\
      --samples       samples.txt \\
      --pass1         \$PASS1_JSONS \\
      --pass2         "\$PASS2_JSON" \\
      --pass3         "\$PASS3_JSON" \\
      --mode          ${params.illumina_barcodetype} \\
      --revtag-orient "\$REVTAG" \\
      --out-summary   Demux_summary.tsv \\
      --out-totals    Demux_totals.tsv

    cp -r json logs/
    cp tags/*.fasta logs/
    rm -rf Strict StrictB Rescued r1_files.txt r1_counts.tsv

    echo -e "\\nDemultiplexing totals:"
    cat Demux_totals.tsv
    echo -e "\\nDemultiplexing finished"
    """
}



// Reorient read pairs by primers (R1 = forward-primer strand), without trimming
// Pairs without any primer are discarded
// (the requirement of both primers is applied later, on merged reads, by `primer_check`)
//
// NB. for paired-end reads, `--revcomp` does not reverse-complement sequences:
//     it re-runs the adapter search with R1 and R2 swapped, and keeps the better-scoring orientation
process reorient_pe {

    label "main_container"

    // cpus 2

    tag "${sampID}"

    input:
      tuple val(sampID), path(reads, stageAs: "input/*")

    output:
      tuple val(sampID), path("Reoriented/${sampID}_R{1,2}.fq.gz"), emit: reads, optional: true
      path "${sampID}_reorient.tsv", emit: stats
      tuple val("${task.process}"), val('cutadapt'), eval('cutadapt --version'), topic: versions

    script:
    """
    echo -e "Reorienting read pairs\\n"
    echo -e "Sample: "  ${sampID}
    echo -e "Input R1: " ${reads[0]}
    echo -e "Input R2: " ${reads[1]}
    echo -e "Forward primer: " ${params.primer_forward}
    echo -e "Reverse primer: " ${params.primer_reverse}

    mkdir -p Reoriented

    ## Primers are not anchored, there could be a linker or a tag in front of them
    ## IUPAC codes are honoured by default
    cutadapt \\
      --revcomp \\
      --action=none \\
      --rename='{header}' \\
      -g "${params.primer_forward};min_overlap=${params.primer_foverlap}" \\
      -G "${params.primer_reverse};min_overlap=${params.primer_roverlap}" \\
      --errors ${params.primer_mismatches} \\
      --no-indels \\
      --discard-untrimmed \\
      --pair-filter=both \\
      --cores ${task.cpus} \\
      --json reorient.json \\
      -o Reoriented/${sampID}_R1.fq.gz \\
      -p Reoriented/${sampID}_R2.fq.gz \\
      ${reads[0]} ${reads[1]} \\
      > reorient.log

    ## Per-sample stats
    python3 - reorient.json ${sampID} > ${sampID}_reorient.tsv <<'PYEOF'
    import json, sys
    rc = json.load(open(sys.argv[1]))["read_counts"]
    print("SampleID\\tInput_Pairs\\tReoriented_Pairs\\tSwapped_Pairs")
    print(f"{sys.argv[2]}\\t{rc['input']}\\t{rc['output']}\\t{rc.get('reverse_complemented') or 0}")
    PYEOF

    cat ${sampID}_reorient.tsv

    ## Remove empty outputs
    find Reoriented -type f -name "*.fq.gz" -size -100c -delete
    if [ ! -f Reoriented/${sampID}_R1.fq.gz ] || [ ! -f Reoriented/${sampID}_R2.fq.gz ]; then
      echo -e "\\nNo read pairs with primers found"
      rm -f Reoriented/*.fq.gz
    fi
    """
}



// Merge paired-end reads (USEARCH or VSEARCH)
// + optional trimming of 3' poly-G tails (two-colour chemistry) prior to merging
process merge_pe {

    label "main_container"

    publishDir "${params.outdir}/01_Demux/Merged", mode: "${params.storagemode}", pattern: "*.fq.gz"
    // cpus 4

    tag "${sampID}"

    input:
      tuple val(sampID), path(reads, stageAs: "input/*")

    output:
      path "${sampID}.fq.gz", emit: merged, optional: true
      tuple val(sampID), path("NotMerged/${sampID}_R1.fq.gz"), path("NotMerged/${sampID}_R2.fq.gz"), emit: notmerged, optional: true
      path "${sampID}_merge.tsv", emit: stats
      tuple val("${task.process}"), val("${params.pe_merger}"), eval(params.pe_merger == "usearch" ? 'usearch --version | sed "s/usearch v//; s/_i86linux64//"' : 'vsearch --version 2>&1 | head -n 1 | sed "s/vsearch //g" | sed "s/,.*//g" | sed "s/^v//" | sed "s/_.*//"'), topic: versions
      tuple val("${task.process}"), val('fastp'), eval(params.qc_polyglen ? 'fastp --version 2>&1 | sed "s/fastp //"' : 'echo "not used"'), topic: versions
      tuple val("${task.process}"), val('seqkit'), eval('seqkit version | sed "s/seqkit v//"'), topic: versions

    script:
    def maxlen_u = params.pe_maxlen ? "-fastq_maxmergelen ${params.pe_maxlen}"   : ""
    def maxlen_v = params.pe_maxlen ? "--fastq_maxmergelen ${params.pe_maxlen}" : ""
    """
    echo -e "Merging paired-end reads\\n"
    echo -e "Sample: "   ${sampID}
    echo -e "Input R1: " ${reads[0]}
    echo -e "Input R2: " ${reads[1]}
    echo -e "Merger: "   ${params.pe_merger}

    mkdir -p NotMerged

    R1=${reads[0]}
    R2=${reads[1]}

    ## Trim poly-G tails (signal loss in two-colour chemistry is called as G)
    if [[ "${params.qc_polyglen ?: ''}" != "" ]]; then
      echo -e "\\nTrimming 3' poly-G tails (min length ${params.qc_polyglen})"
      fastp \\
        --in1 \$R1 --in2 \$R2 \\
        --out1 polyg_R1.fq.gz --out2 polyg_R2.fq.gz \\
        --trim_poly_g --poly_g_min_len ${params.qc_polyglen} \\
        --disable_adapter_trimming \\
        --disable_quality_filtering \\
        --length_required 20 \\
        --dont_eval_duplication \\
        --thread ${task.cpus} \\
        --json fastp.json --html fastp.html \\
        2> fastp.log
      R1=polyg_R1.fq.gz
      R2=polyg_R2.fq.gz
    fi

    if [[ "${params.pe_merger}" == "usearch" ]]; then

      if ! command -v usearch > /dev/null 2>&1; then
        echo -e "\\nERROR: USEARCH is not available in the environment."
        echo -e "Use '--pe_merger vsearch', or a container with USEARCH installed.\\n"
        exit 1
      fi

      ## USEARCH works only with uncompressed files
      gunzip -c \$R1 > unpacked_R1.fq
      gunzip -c \$R2 > unpacked_R2.fq

      ## Staggered pairs are merged, overhangs are trimmed (by default)
      usearch \\
        -fastq_mergepairs  unpacked_R1.fq \\
        -reverse           unpacked_R2.fq \\
        -fastq_maxdiffs    ${params.pe_maxdiffs} \\
        -fastq_pctid       ${params.pe_pctid} \\
        -fastq_minovlen    ${params.pe_minoverlap} \\
        -fastq_minmergelen ${params.pe_minlen} \\
        ${maxlen_u} \\
        -fastq_qmax 93 \\
        -threads ${task.cpus} \\
        -fastqout               merged.fq \\
        -fastqout_notmerged_fwd notmerged_R1.fq \\
        -fastqout_notmerged_rev notmerged_R2.fq \\
        -report                 merge_report.txt \\
        > merge.log 2>&1

      rm unpacked_R1.fq unpacked_R2.fq

    else

      vsearch \\
        --fastq_mergepairs \$R1 \\
        --reverse          \$R2 \\
        --fastq_maxdiffs   ${params.pe_maxdiffs} \\
        --fastq_maxdiffpct ${100 - (params.pe_pctid as Integer)} \\
        --fastq_minovlen   ${params.pe_minoverlap} \\
        --fastq_minmergelen ${params.pe_minlen} \\
        ${maxlen_v} \\
        --fastq_allowmergestagger \\
        --fastq_qmax 93 \\
        --threads ${task.cpus} \\
        --fastqout               merged.fq \\
        --fastqout_notmerged_fwd notmerged_R1.fq \\
        --fastqout_notmerged_rev notmerged_R2.fq \\
        --no_progress \\
        2> merge_report.txt

    fi

    cat merge_report.txt

    ## Keep only sequence IDs in headers (no spaces)
    seqkit seq --only-id merged.fq          | pigz -p ${task.cpus} -${params.gzip_compression} > ${sampID}.fq.gz
    seqkit seq --only-id notmerged_R1.fq    | pigz -p ${task.cpus} -${params.gzip_compression} > NotMerged/${sampID}_R1.fq.gz
    seqkit seq --only-id notmerged_R2.fq    | pigz -p ${task.cpus} -${params.gzip_compression} > NotMerged/${sampID}_R2.fq.gz

    ## Per-sample stats
    NIN=\$(( \$(gunzip -c ${reads[0]} | wc -l) / 4 ))
    NMERGED=\$(( \$(wc -l < merged.fq) / 4 ))
    NNOT=\$(( \$(wc -l < notmerged_R1.fq) / 4 ))
    printf "SampleID\\tMerge_Input_Pairs\\tMerged_Reads\\tNotMerged_Pairs\\n" > ${sampID}_merge.tsv
    printf "%s\\t%s\\t%s\\t%s\\n" ${sampID} \$NIN \$NMERGED \$NNOT >> ${sampID}_merge.tsv
    cat ${sampID}_merge.tsv

    rm merged.fq notmerged_R1.fq notmerged_R2.fq
    rm -f polyg_R1.fq.gz polyg_R2.fq.gz

    ## Remove empty outputs
    if [ \$NMERGED -eq 0 ]; then rm ${sampID}.fq.gz; fi
    if [ \$NNOT -eq 0 ];    then rm NotMerged/${sampID}_R1.fq.gz NotMerged/${sampID}_R2.fq.gz; fi
    """
}







// Demultiplexing with cutadapt - for Illumina PE reads (only not merged)
// NB. it's possible to use anchored adapters (e.g., -g ^file:barcodes.fa),
//     but there could be a preceding nucleotides before the barcode,
//     therefore, modified barcodes would be used here (e.g., XN{30})
process demux_pe {

    label "main_container"

    publishDir "${params.outdir}/01_Demux", mode: 'symlink'
    // cpus 20

    input:
      tuple path(input_R1), path(input_R2)
      path barcodes   // barcodes_modified.fa (e.g., XN{30})

    output:
      tuple path("Combined/*.R1.fastq.gz"), path("Combined/*.R2.fastq.gz"), emit: demux_pe, optional: true
      path "NonMerged_samples.txt", emit: samples_nonm_pe, optional: true

    script:
    """
    echo -e "\nDemultiplexing not-merged reads"

    echo -e "Input R1: " ${input_R1}
    echo -e "Input R2: " ${input_R2}
    echo -e "Barcodes: " ${barcodes}

    ## First round
    echo -e "\nRound 1:"

    cutadapt -g file:${barcodes} \
      -o round1-{name}.R1.fastq.gz \
      -p round1-{name}.R2.fastq.gz \
      --errors ${params.barcode_errors} \
      --overlap ${params.barcode_overlap} \
      --no-indels \
      --cores ${task.cpus} \
      ${input_R1} ${input_R2} \
      > cutadapt_round_1.log

    echo -e ".. round 1 finished"
    
    ## Second round
    echo -e "\nRound 2:"

    cutadapt -g file:${barcodes} \
      -o round2-{name}.R2.fastq.gz \
      -p round2-{name}.R1.fastq.gz \
      --errors ${params.barcode_errors} \
      --overlap ${params.barcode_overlap} \
      --no-indels \
      --cores ${task.cpus} \
      round1-unknown.R2.fastq.gz round1-unknown.R1.fastq.gz \
      > cutadapt_round_2.log

    echo -e ".. round 2 finished"

    ## Remove empty files (no sequences)
    echo -e "\nRemoving empty files"
    find . -type f -name "round*.fastq.gz" -size -29c -print -delete
    echo -e "..Done"

    ## Remove unknowns
    echo -e "Removing unknowns"
    rm round1-unknown.R{1,2}.fastq.gz
    rm round2-unknown.R{1,2}.fastq.gz

    ## Combine sequences from round 1 and round 2 for each sample
    echo -e "\nCombining sequences from round 1 and round 2 for each sample"

    if test -n "\$(find . -maxdepth 1 -name 'round*.fastq.gz' -print -quit)"
    then

      mkdir -p Combined

      find . -name "round*.R1.fastq.gz" | sort | parallel -j1 \
        "cat {} >> Combined/{= s/round1-//; s/round2-// =}"

      find . -name "round*.R2.fastq.gz" | sort | parallel -j1 \
        "cat {} >> Combined/{= s/round1-//; s/round2-// =}"

      ## Write sample names to the file
      find Combined -name "*.R1.fastq.gz" \
        | sed 's/Combined\\///; s/\\.R1\\.fastq\\.gz//' \
        > NonMerged_samples.txt

      echo -e "..Done"

    else
      echo "..No files"
    fi

    ## Clean up
    echo -e "..Removing temporary files"
    find . -type f -name "round*.fastq.gz" -print -delete



    echo -e "\nDemultiplexing finished"
    """
}


// Trim primers of nonmerged PE reads
// + Estimate sequence qualities
process trim_primers_pe {

    label "main_container"

    publishDir "${params.outdir}/03_PrimerTrim", mode: 'symlink'
    // cpus 2

    // Add sample ID to the log file
    tag "${input.getSimpleName()}"

    input:
      path input   // tuple of size 2

    output:
      path "${input.getSimpleName()}_R1.fa.gz", emit: primertrimmed_fa_R1, optional: true
      path "${input.getSimpleName()}_R2.fa.gz", emit: primertrimmed_fa_R2, optional: true
      path "${input.getSimpleName()}_hash_table_R1.txt.gz", emit: hashes_R1, optional: true
      path "${input.getSimpleName()}_hash_table_R2.txt.gz", emit: hashes_R2, optional: true
      path "${input.getSimpleName()}_R1.fq.gz", emit: primertrimmed_fq_R1, optional: true
      path "${input.getSimpleName()}_R2.fq.gz", emit: primertrimmed_fq_R2, optional: true
      path "${input.getSimpleName()}_uc_R1.uc.gz", emit: ucR1, optional: true
      path "${input.getSimpleName()}_uc_R2.uc.gz", emit: ucR2, optional: true

    script:
    sampID="${input.getSimpleName()}"
    
    """
    echo -e "Input sample: " ${sampID}
    echo -e "Forward primer: " ${params.primer_forward}
    echo -e "Reverse primer: " ${params.primer_reverse}

    ## Reverse-complement primers
    FR=\$(rc.sh ${params.primer_forward})
    RR=\$(rc.sh ${params.primer_reverse})
    
    echo -e "Forward primer RC: " "\$FR"
    echo -e "Reverse primer RC: " "\$RR"

    ## Discard sequences without both primers
    echo -e "\nChecking primers"
    
    echo -e "..Forward strain"

    cutadapt \
      -a ${params.primer_forward}";required;min_overlap=${params.primer_foverlap}"..."\$RR"";required;min_overlap=${params.primer_roverlap}" \
      --errors ${params.primer_mismatches} \
      --cores ${task.cpus} \
      --action=none \
      --discard-untrimmed \
      -o for_R1.fastq.gz -p for_R2.fastq.gz \
      ${input[0]} ${input[1]} \
      > cutadapt_1.log

    
    echo -e "..Reverse strain"

    cutadapt \
      -a ${params.primer_reverse}";required;min_overlap=${params.primer_roverlap}"..."\$FR"";required;min_overlap=${params.primer_foverlap}" \
      --errors ${params.primer_mismatches} \
      --cores ${task.cpus} \
      --action=none \
      --discard-untrimmed \
      -p rev_R1.fastq.gz -o rev_R2.fastq.gz \
      ${input[0]} ${input[1]} \
      > cutadapt_2.log
    
    # cutadapt \
    #   -a FWDPRIMER...RCREVPRIMER \
    #   -A REVPRIMER...RCFWDPRIMER \
    #   --discard-untrimmed \
    #   -o out.1.fastq.gz -p out.2.fastq.gz \
    #   in.1.fastq.gz in.2.fastq.gz


    echo -e "\nReorienting"

    if [ -s for_R1.fastq.gz ]; then
      zcat for_R1.fastq.gz | seqkit replace -p "\\s.+" | gzip -7 > OK_R1.fastq.gz
      zcat for_R2.fastq.gz | seqkit replace -p "\\s.+" | gzip -7 > OK_R2.fastq.gz
    fi

    if [ -s rev_R1.fastq.gz ]; then
      echo -e "..Adding sequences to the main pool"
      zcat rev_R1.fastq.gz | seqkit replace -p "\\s.+" | gzip -7 >> OK_R1.fastq.gz
      zcat rev_R2.fastq.gz | seqkit replace -p "\\s.+" | gzip -7 >> OK_R2.fastq.gz

    else
      echo -e "..Probably all sequences are in forward orientation"
    fi


    echo -e "\nTrimming primers"
    if [ -s OK_R1.fastq.gz]; then

      cutadapt \
        -a ${params.primer_forward}";required;min_overlap=${params.primer_foverlap}"..."\$RR"";required;min_overlap=${params.primer_roverlap}" \
        --errors ${params.primer_mismatches} \
        --cores ${task.cpus} \
        --action=trim \
        --discard-untrimmed \
        --minimum-length ${params.trim_minlen} \
        --output ${sampID}_R1.fq.gz --paired-output ${sampID}_R2.fq.gz \
        OK_R1.fastq.gz OK_R2.fastq.gz

    fi


    ## Quality estimation and dereplication

    if [ -s ${sampID}_R1.fq.gz ]; then 

      ## Estimate sequence quality (for the extracted region)
      ## Sequence ID - Hash - Length - Average Phred score
      echo -e "\nCreating sequence hash table with average sequence quality"
      
      seqkit fx2tab --length --avg-qual ${sampID}_R1.fq.gz \
        | hash_sequences.sh \
        | awk '{print \$1 "\t" \$6 "\t" \$4 "\t" \$5}' \
        > tmp_hash_table_R1.txt

      seqkit fx2tab --length --avg-qual ${sampID}_R2.fq.gz \
        | hash_sequences.sh \
        | awk '{print \$1 "\t" \$6 "\t" \$4 "\t" \$5}' \
        > tmp_hash_table_R2.txt
      
      echo -e "..Done"


      ## Estimating MaxEE
      echo -e "\nEstimating maximum number of expected errors per sequence"

      vsearch \
          --fastx_filter ${sampID}_R1.fq.gz \
          --fastq_qmax 93 \
          --eeout \
          --fastaout - \
        | seqkit seq --name \
        | sed 's/;ee=/\t/g' \
        > tmp_ee_R1.txt

      vsearch \
          --fastx_filter ${sampID}_R2.fq.gz \
          --fastq_qmax 93 \
          --eeout \
          --fastaout - \
        | seqkit seq --name \
        | sed 's/;ee=/\t/g' \
        > tmp_ee_R2.txt

      echo -e "..Done"


      echo -e "\nMerging quality estimates"

      max_ee.R \
        tmp_hash_table_R1.txt \
        tmp_ee_R1.txt \
        ${sampID}_hash_table_R1.txt

      max_ee.R \
        tmp_hash_table_R2.txt \
        tmp_ee_R2.txt \
        ${sampID}_hash_table_R2.txt

      echo -e "..Done"


      ## Independent dereplication of pair-end reads
      echo -e "\nDereplicating R1 and R2 (independently)"
      
      seqkit fq2fa -w 0 ${sampID}_R1.fq.gz \
      | vsearch \
        --derep_fulllength - \
        --output - \
        --strand both \
        --fasta_width 0 \
        --threads 1 \
        --relabel_sha1 \
        --sizein --sizeout \
        --uc ${sampID}_uc_R1.uc \
        --quiet \
      | gzip -7 \
      > ${sampID}_R1.fa.gz

      seqkit fq2fa -w 0 ${sampID}_R2.fq.gz \
      | vsearch \
        --derep_fulllength - \
        --output - \
        --strand both \
        --fasta_width 0 \
        --threads 1 \
        --relabel_sha1 \
        --sizein --sizeout \
        --uc ${sampID}_uc_R2.uc \
        --quiet \
      | gzip -7 \
      > ${sampID}_R2.fa.gz


      echo -e "..Done"

      ## Compress results
      echo -e "Compressing result"
      gzip -7 ${sampID}_hash_table_R1.txt
      gzip -7 ${sampID}_hash_table_R2.txt
      gzip -7 ${sampID}_uc_R1.uc
      gzip -7 ${sampID}_uc_R2.uc


    else
      echo -e "\nNo sequences found after primer removal"
    fi

    ## Clean up
    if [ -f for_R1.fastq.gz ]; then rm for_R1.fastq.gz; fi
    if [ -f for_R2.fastq.gz ]; then rm for_R2.fastq.gz; fi
    if [ -f rev_R1.fastq.gz ]; then rm rev_R1.fastq.gz; fi
    if [ -f rev_R2.fastq.gz ]; then rm rev_R2.fastq.gz; fi
    if [ -f OK_R1.fastq.gz  ]; then rm OK_R1.fastq.gz; fi
    if [ -f OK_R2.fastq.gz  ]; then rm OK_R2.fastq.gz; fi
  
    echo -e "..Done"

    """
}




// Combine paired reads into single sequences 
// by reverse-complementing the reverse read and inserting poly-N padding
// + Estimate sequence qualities (without N pads!)
process join_pe {

    label "main_container"

    // publishDir "${out_1_joinPE}", mode: "${params.storagemode}"
    // cpus 2

    // Add sample ID to the log file
    tag "${input}"

    input:
      val input       // Sample name "(e.g., Barcode07_1__IS859)"
      path all_samples

    output:
      path "${input}_JoinedPE.fq.gz", emit: jj_FQ, optional: true
      path "${input}_JoinedPE_hash_table.txt.gz", emit: jj_hashes, optional: true

    script:
    sampID="${input}"
    
    """
    echo -e "Joining non-merged Illumina reads"
    echo -e "Input sample: " ${sampID}

    echo -e "\\nJoining with N-pads"
    vsearch \
      --fastq_join ${input}.R1.fastq.gz \
      --reverse ${input}.R2.fastq.gz \
      --join_padgap ${params.illumina_joinpadgap} \
      --join_padgapq ${params.illumina_joinpadqual} \
      --fastqout - \
    | seqkit replace -p "\\s.+" \
    | gzip -${params.gzip_compression} \
    > ${sampID}_JoinedPE.fq.gz

    ## Check if there are some sequences in the file
    if [ -n "\$(find . -name ${sampID}_JoinedPE.fq.gz -prune -size +29c)" ]; then

      echo -e "\\nJoining without N-pads (for quality estimation)"
      vsearch \
        --fastq_join ${input}.R1.fastq.gz \
        --reverse ${input}.R2.fastq.gz \
        --join_padgap "" \
        --join_padgapq "" \
        --fastqout - \
      | seqkit replace -p "\\s.+" \
      | gzip -${params.gzip_compression} \
      > tmp_for_qual.fq.gz


      ## Estimate sequence quality (without N pads!)
      ## Sequence ID - Hash - Length - Average Phred score
      echo -e "\\nCreating sequence hash table with average sequence quality"
        
      seqkit fx2tab --length --avg-qual tmp_for_qual.fq.gz \
        | hash_sequences.sh \
        | awk '{print \$1 "\t" \$6 "\t" \$4 "\t" \$5}' \
        > tmp_hash_table.txt
      
      echo -e "..Done"

      ## Estimating MaxEE
      echo -e "\\nEstimating maximum number of expected errors per sequence"

      vsearch \
          --fastx_filter tmp_for_qual.fq.gz \
          --fastq_qmax 93 \
          --eeout \
          --fastaout - \
        | seqkit seq --name \
        | sed 's/;ee=/\t/g' \
        > tmp_ee.txt

      echo -e "..Done"

      echo -e "\\nMerging quality estimates"

      max_ee.R \
        tmp_hash_table.txt \
        tmp_ee.txt \
        ${sampID}_JoinedPE_hash_table.txt

      echo -e "..Done"

      ## Compress results
      gzip -${params.gzip_compression} ${sampID}_JoinedPE_hash_table.txt

      ## Clean up
      rm tmp_for_qual.fq.gz
      rm tmp_hash_table.txt
      rm tmp_ee.txt
    
    else
      echo -e "\\nIt looks like there are no joined reads"
    fi

    ## Remove redundant symlinks
    find -L . -name "*.fastq.gz" | grep -v ${input} | parallel -j1 "rm {}"

    """
}

