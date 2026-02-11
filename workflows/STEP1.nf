/*
============================================================================
  NextITS: Pipeline to process eukaryotic ITS amplicons
============================================================================
  License: Apache-2.0
  Github : https://github.com/vmikk/NextITS
  Website: https://Next-ITS.github.io/
----------------------------------------------------------------------------
*/

// ---- Step-1 workflow ----


// Include functions
include { software_versions_to_yaml } from '../modules/version_parser.nf'
include { dumpParamsTsv }             from '../modules/dump_parameters.nf'
include { CHIMERA_REMOVAL }           from '../subworkflows/chimera_removal_subworkflow.nf'

if ( params.seqplatform == "Illumina" ){
  include { qc_pe; merge_pe; demux_pe; trim_primers_pe; join_pe } from '../modules/Illumina_pe.nf'
}


// Define output paths for different steps
out_0_bam    = params.outdir + "/00_BAM2FASTQ"
out_1_demux  = params.outdir + "/01_Demux"
out_1_joinPE = params.outdir + "/01_JoinedPE"
out_2_primer = params.outdir + "/02_PrimerCheck"
out_3_itsx   = params.outdir + "/03_ITSx"
out_3_itsxp  = params.outdir + "/03_ITSx_PooledParts"
out_3_trim   = params.outdir + "/03_PrimerTrim"
out_3_trimPE = params.outdir + "/03_PrimerTrim_NotMerged"
out_4_homop  = params.outdir + "/04_Homopolymer"
// out_5_chim   = params.outdir + "/05_Chimera"
out_6_tj     = params.outdir + "/06_TagJumpFiltration"
out_7_seq    = params.outdir + "/07_SeqTable"
out_8_smr    = params.outdir + "/08_RunSummary"
out_9_db     = params.outdir + "/09_DB"
out_tracedir = params.tracedir

// Sub-workflow-specific outputs
out_3_quickstats = params.outdir + "/03_Stats"


// Convert BAM to FASTQ
process bam2fastq {

    label "main_container"
    publishDir "${out_0_bam}", mode: "${params.storagemode}"

    // cpus 2

    input:
      path input
      path bam_index

    output:
      path "*.fastq.gz", emit: fastq, optional: false
      tuple val("${task.process}"), val('bam2fastq'), eval('bam2fastq --version | head -n 1 | sed "s/bam2fastq //"'), topic: versions

    script:
    """
    echo -e "Converting BAM to FASTQ\\n"
    echo -e "Input file: " ${input}
    echo -e "BAM index: "  ${bam_index}

    bam2fastq \
      -c ${params.gzip_compression} \
      --num-threads ${task.cpus} \
      ${input}

    echo -e "\\nConvertion finished"
    """
}


// Quality filtering for single-end reads
process qc_se {

    label "main_container"

    // cpus 10

    // Add file ID to the log file
    tag "${input.getSimpleName()}"

    input:
      path input

    output:
      path "${input.getSimpleName()}.fq.gz", emit: filtered, optional: true
      tuple val("${task.process}"), val('vsearch'), eval('vsearch --version 2>&1 | head -n 1 | sed "s/vsearch //g" | sed "s/,.*//g" | sed "s/^v//" | sed "s/_.*//"'), topic: versions
      tuple val("${task.process}"), val('seqkit'), eval('seqkit version | sed "s/seqkit v//"'), topic: versions

    script:
    filter_maxee      = params.qc_maxee      ? "--fastq_maxee ${params.qc_maxee}"          : ""
    filter_maxeerate  = params.qc_maxeerate  ? "--fastq_maxee_rate ${params.qc_maxeerate}" : ""
    """
    echo -e "QC\\n"
    echo -e "Input file: " ${input}

    ## We do not need to change the file name (output name should be the same as input)
    ## Therefore, temporary rename input
    mv ${input} inp.fq.gz

    vsearch \
      --fastq_filter inp.fq.gz \
      --fastq_qmax 93 \
      ${filter_maxee} \
      ${filter_maxeerate} \
      --fastq_maxns ${params.qc_maxn} \
      --threads ${task.cpus} \
      --fastqout - \
    | seqkit grep \
      --by-seq --ignore-case --invert-match --only-positive-strand --use-regexp -w 0 \
      --pattern '"(A{${params.qc_maxhomopolymerlen},}|C{${params.qc_maxhomopolymerlen},}|T{${params.qc_maxhomopolymerlen},}|G{${params.qc_maxhomopolymerlen},})"' \
    | gzip -${params.gzip_compression} \
    > "${input.getSimpleName()}.fq.gz"

    ## qc_maxhomopolymerlen
    # e.g. "(A{26,}|C{26,}|T{26,}|G{26,})"

    echo -e "\\nQC finished"
    """
}


// Validate tags for demultiplexing
process tag_validation {

    label "main_container"
    // cpus 1

    publishDir "${out_1_demux}", pattern: "tag_names_renamed.tsv", mode: "${params.storagemode}"

    input:
      path barcodes

    output:
      path "barcodes_validated.fasta", emit: fasta
      path "biosamples_asym.csv",      emit: biosamples_asym, optional: true
      path "biosamples_sym.csv",       emit: biosamples_sym,  optional: true
      path "file_renaming.tsv",        emit: file_renaming,   optional: true
      path "unknown_combinations.tsv", emit: unknown_combinations, optional: true
      path "tag_names_renamed.tsv",    emit: tag_names_renamed, optional: true

    script:
    """
    echo -e "Valdidating demultiplexing tags\\n"
    echo -e "Input file: " ${barcodes}

    ## Convert Windows-style line endings (CRLF) to Unix-style (LF)
    LC_ALL=C sed -i 's/\r\$//g' ${barcodes}

    ## Perform tag validation
    validate_tags.R \
      --tags   ${barcodes} \
      --output barcodes_validated.fasta

    echo -e "\\nTag validation finished"
    """
}



// Demultiplexing with LIMA - for PacBio reads
process demux {

    label "main_container"

    publishDir "${out_1_demux}", mode: "${params.storagemode}"  // , saveAs: { filename -> "foo_$filename" }
    // cpus 10

    input:
      path input_fastq
      path barcodes
      path biosamples_sym       // for dual or asymmetric barcodes
      path biosamples_asym      // for dual or asymmetric barcodes
      path file_renaming        // for dual or asymmetric barcodes
      path unknown_combinations // for dual or asymmetric barcodes

    output:
      path "LIMA/*.fq.gz",             emit: samples_demux
      path "LIMA/lima.lima.report.gz", emit: lima_report
      path "LIMA/lima.lima.counts",    emit: lima_counts
      path "LIMA/lima.lima.summary",   emit: lima_summary
      tuple val("${task.process}"), val('lima'), eval('lima --version | head -n 1 | sed "s/lima //"'), topic: versions
      tuple val("${task.process}"), val('brename'), eval('brename --help | head -n 4 | tail -1 | sed "s/Version: //"'), topic: versions

    script:
    """
    echo -e "Input file: " ${input_fastq}
    echo -e "Barcodes: "   ${barcodes}

    ## Directory for the results
    mkdir -p LIMA
    
    echo -e "Validating data\n"

    ## Check if symmetric barcodes were provided in the `...` format
    ## (if `biosamples_sym` does not exists, it means that it is a dummy file)
    ## (if exists, it means that tags were split into sym and asym at the tag validation step)
    if [[ ${params.lima_barcodetype} = "dual_symmetric" ]] && [ -e ${biosamples_sym} ] ; then
        echo -e "\\nERROR: Symmetric tags are provided in '...' format.\\n"
        echo -e "In the FASTA file, please include only one tag per sample, since these tags are identical.\\n"
        exit 1
    fi

    ## Count the number of samples in Biosample files - only for `dual` and `dual_asymmetric` barcodes
    if [[ ${params.lima_barcodetype} == "dual_asymmetric" ]] || [[ ${params.lima_barcodetype} == "dual" ]]; then

      if [ ! -e ${biosamples_asym} ]; then
        
        echo -e "\\nERROR: Tags are specified in wrong format"
        echo -e "Use the '...' format in FASTA file.\\n"
        exit 1
      
      else
        line_count_sym=\$(wc  -l < ${biosamples_sym})
        line_count_asym=\$(wc -l < ${biosamples_asym})

        echo -e "..Number of lines in symmetric file: "  \$line_count_sym
        echo -e "..Number of lines in asymmetric file: " \$line_count_asym

        ## Check the presence of dual barcode combinations
        ## If line count is less than 2, it means there are no samples specified
        if [[ ${params.lima_barcodetype} == "dual_asymmetric" ]] && [[ \$line_count_asym -lt 2 ]]; then
          echo -e "\\nERROR: No asymmetric barcodes detected for demultiplexing.\\n"
          return 1
        fi

        if [[ ${params.lima_barcodetype} == "dual" ]] && [[ \$line_count_asym -lt 2 ]]; then
          echo -e "\\nWARNING: No asymmetric barcodes detected, consider using '--lima_barcodetype dual_symmetric'.\\n"
        fi

      fi  # end of missing asym biosamples

    fi    # end of dual/asym validation



    ## Combine shared arguments into a single variable
    ## Note the array syntax - that's because of LIMA parser error messages
    ## (note also that it works in bash, but not in zsh)
    common_args=("--ccs \
      --window-size  ${params.lima_windowsize} \
      --min-length   ${params.lima_minlen} \
      --min-score    ${params.lima_minscore} \
      --min-ref-span ${params.lima_minrefspan} \
      --split-named \
      --num-threads ${task.cpus} \
      --log-level INFO \
      ${input_fastq} \
      ${barcodes}")


    ## Demultiplex, depending on the barcode type selected
    case ${params.lima_barcodetype} in

      "single")
        echo -e "\\nDemultiplexing with LIMA (single barcode)"
        lima --same --single-side \
          --log-file LIMA/_log.txt \
          \$common_args \
          "LIMA/lima.fq.gz"
        ;;

      "dual_symmetric")
        echo -e "\\nDemultiplexing with LIMA (dual symmetric barcodes)"
        lima --same \
          --min-end-score       ${params.lima_minendscore} \
          --min-scoring-regions ${params.lima_minscoringregions} \
          --log-file LIMA/_log.txt \
          \$common_args \
          "LIMA/lima.fq.gz"
        ;;

      "dual_asymmetric")
        echo -e "\\nDemultiplexing with LIMA (dual asymmetric barcodes)"
        lima --different \
          --min-end-score       ${params.lima_minendscore} \
          --min-scoring-regions ${params.lima_minscoringregions} \
          --biosample-csv       ${biosamples_asym} \
          --log-file LIMA/_log.txt \
          \$common_args \
          "LIMA/lima.fq.gz"
        ;;

      "dual")
        mkdir -p LIMAs LIMAd

        if [[ \$line_count_sym -ge 2 ]]; then
        echo -e "\\nDemultiplexing with LIMA (dual symmetric barcodes)"
        lima --same \
          --min-end-score       ${params.lima_minendscore} \
          --min-scoring-regions ${params.lima_minscoringregions} \
          --biosample-csv       ${biosamples_sym} \
          --log-file LIMAs/_log.txt \
          \$common_args \
          "LIMAs/lima.fq.gz"
        fi

        if [[ \$line_count_asym -ge 2 ]]; then
        echo -e "\\nDemultiplexing with LIMA (dual asymmetric barcodes)"
        lima --different \
          --min-end-score       ${params.lima_minendscore} \
          --min-scoring-regions ${params.lima_minscoringregions} \
          --biosample-csv       ${biosamples_asym} \
          --log-file LIMAd/_log.txt \
          \$common_args \
          "LIMAd/lima.fq.gz"
        fi
        ;;
    esac


    ## Combining symmetric and asymmetric files
    if [ ${params.lima_barcodetype} = "dual" ]; then

      echo -e "\\nPooling of symmetric and asymmetric barcodes"
      cd LIMA
      find ../LIMAd -name "*.fq.gz" | parallel -j1 "ln -s {} ."
      find ../LIMAs -name "*.fq.gz" | parallel -j1 "ln -s {} ."
      cd ..

    fi


    ## Rename barcode combinations into sample names
    ## Only user-provided combinations whould be kept (based on `lima --biosample-csv`)
    if [[ ${params.lima_barcodetype} == "dual_asymmetric" ]] || [[ ${params.lima_barcodetype} == "dual" ]]; then

      echo -e "\\n..Renaming files from tag IDs to sample names"
      brename -p "(.+)" -r "{kv}" -k ${file_renaming} LIMA/

      echo -e "\\n..Checking for unknown tag combinations"
      echo -e "\\n...Number of unknowns detected:"
      find LIMA -name "lima.*.fq.gz" | wc -l
      
      if [[ ${params.lima_remove_unknown} == "false" ]]; then

        if [ -s ${unknown_combinations} ]; then
          echo -e "\\n...Renaming unknown combinations"
          brename -p "(.+)" -r "{kv}" -k ${unknown_combinations} LIMA/
        else
          echo -e "\\n...No unknown combinations require renaming"
        fi

        echo -e "\\n...Number of unknowns remained:"
        find LIMA -name "lima.*.fq.gz" | wc -l

      fi

      echo -e "\\n...Removing unknowns:"
      find LIMA -name "lima.*.fq.gz" | parallel -j1 "echo {} && rm {}"

    fi  # end of dual/asym renaming

    if [[ ${params.lima_barcodetype} == "dual_symmetric" ]] || [[ ${params.lima_barcodetype} == "single" ]]; then

      echo -e "\\n..Renaming demultiplexed files"
      rename --filename \
        's/^lima.//g; s/--.*\$/.fq.gz/' \
        \$(find LIMA -name "*.fq.gz")
    
    fi


    ## Combine summary stats for dual barcodes (two LIMA runs)
    if [[ ${params.lima_barcodetype} == "dual" ]]; then

      echo -e "\\n..Combining dual-barcode log files"

      if [ -f "LIMAd/lima.lima.summary" ]; then
        echo -e "Asymmetric barcodes summary\\n\\n" >> LIMA/lima.lima.summary
        cat LIMAd/lima.lima.summary >> LIMA/lima.lima.summary

        echo -e "Asymmetric barcodes counts\\n\\n" >> LIMA/lima.lima.counts
        cat LIMAd/lima.lima.counts >> LIMA/lima.lima.counts

        echo -e "Asymmetric barcodes report\\n\\n" >> LIMA/lima.lima.report
        cat LIMAd/lima.lima.report >> LIMA/lima.lima.report
      fi

      if [ -f "LIMAs/lima.lima.summary" ]; then
        echo -e "\\n\\nSymmetric barcodes summary\\n\\n" >> LIMA/lima.lima.summary
        cat LIMAs/lima.lima.summary >> LIMA/lima.lima.summary

        echo -e "\\n\\nSymmetric barcodes counts\\n\\n" >> LIMA/lima.lima.counts
        cat LIMAs/lima.lima.counts >> LIMA/lima.lima.counts

        ## Reports should be identical for symmetric and asymmetric barcodes, so no need to combine them
        # echo -e "\\n\\nSymmetric barcodes report\\n\\n" >> LIMA/lima.lima.report
        # cat LIMAs/lima.lima.report >> LIMA/lima.lima.report
      fi

    fi  # end of dual logs pooling


    ## Compress logs
    echo -e "..Compressing log file"
    gzip -${params.gzip_compression} LIMA/lima.lima.report


    ## LIMA defaults:
    # SYMMETRIC  : --ccs --min-score 0 --min-end-score 80 --min-ref-span 0.75 --same --single-end
    # ASYMMETRIC : --ccs --min-score 80 --min-end-score 50 --min-ref-span 0.75 --different --min-scoring-regions 2

    echo -e "\\nDemultiplexing finished"
    """
}


// Merge Illumina PE reads
process merge_pe {

    label "main_container"

    // publishDir "${out_1_demux}", mode: "${params.storagemode}"
    // cpus 10

    input:
      path input_R1
      path input_R2

    output:
      path "Merged.fq.gz", emit: r12
      tuple path("NotMerged_R1.fq.gz"), path("NotMerged_R2.fq.gz"), emit: nm, optional: true

    script:
    """
    echo -e "Merging Illumina pair-end reads\\n"

    ## By default, fastp modifies sequences header
    ## e.g., `merged_150_15` means that 150bp are from read1, and 15bp are from read2
    ## But we'll preserve only sequence ID

    fastp \
      --in1 ${input_R1} \
      --in2 ${input_R2} \
      --merge --correction \
      --overlap_len_require ${params.pe_minoverlap} \
      --overlap_diff_limit ${params.pe_difflimit} \
      --overlap_diff_percent_limit ${params.pe_diffperclimit} \
      --length_required ${params.pe_minlen} \
      --disable_quality_filtering \
      --disable_adapter_trimming \
      --dont_eval_duplication \
      --compression 6 \
      --thread ${task.cpus} \
      --out1 NotMerged_R1.fq.gz \
      --out2 NotMerged_R2.fq.gz \
      --json log.json \
      --html log.html \
      --stdout \
    | seqkit seq --only-id \
    | gzip -${params.gzip_compression} \
    > Merged.fq.gz

    ##  --merged_out Merged.fq.gz \
    ##  --n_base_limit ${params.pe_nlimit} \

    # --overlap_len_require         the minimum length to detect overlapped region of PE reads
    # --overlap_diff_limit          the maximum number of mismatched bases to detect overlapped region of PE reads
    # --overlap_diff_percent_limit  the maximum percentage of mismatched bases to detect overlapped region of PE reads
    ## NB: reads should meet these three conditions simultaneously!

    echo -e "..done"
    """
}


// Modify barcodes for cutadapt (restrict the search window)
process prep_barcodes {

    label "main_container"

    // publishDir "${out_1_demux}", mode: "${params.storagemode}"
    // cpus 1

    input:
      path barcodes

    output:
      path "barcodes_modified.fa", emit: barcodesm

    script:
    """
    echo -e "Restricting the search window for barcode lookup"
    echo -e "Provided barcodes: " ${barcodes}

    ## Add `XN{30}` to the barcodes

    sed -e '/^>/! s/^/XN{${params.barcode_window}}/' \
      ${barcodes} \
      > barcodes_modified.fa

    echo -e "..Done"
    """
}


// Demultiplexing with cutadapt - for Illumina SE reads
process demux_illumina {

    label "main_container"

    publishDir "${out_1_demux}", mode: "${params.storagemode}"
    // cpus 10

    input:
      path input_fastq
      path barcodes

    output:
      path "*.fq.gz", emit: samples_demux

    script:
    """
    echo -e "Input file: " ${input_fastq}
    echo -e "Barcodes: " ${barcodes}

    echo -e "\nDemultiplexing with cutadapt:"

    ## Demultiplex with cutadapt
    cutadapt \
      -g file:${barcodes} \
      --revcomp --rename "{header}" \
      --errors ${params.barcode_errors} \
      --overlap ${params.barcode_overlap} \
      --no-indels \
      --cores ${task.cpus} \
      --discard-untrimmed \
      --action none \
      -o "{name}.fq.gz" \
      ${input_fastq} \
      > cutadapt.log

    echo -e "\\n..done"

    ## Remove empty files (no sequences)
    echo -e "\\nRemoving empty files"
    find . -type f -name "*.fq.gz" -size -29c -print -delete
    echo -e "..Done"

    echo -e "\\nDemultiplexing finished"
    """
}



// Primer disambiguation
process disambiguate {

    label "main_container"

    // publishDir "${out_2_primer}", mode: "${params.storagemode}"
    // cpus 1

    output:
      path "primer_F.fasta",  emit: F
      path "primer_R.fasta",  emit: R
      path "primer_Fr.fasta", emit: Fr
      path "primer_Rr.fasta", emit: Rr

    script:

    """

    ## Disambiguate forward primer
    echo -e "Disambiguating forward primer"
    disambiguate_primers.R \
      ${params.primer_forward} \
      primer_F.fasta

    ## Disambiguate reverse primer
    echo -e "\\nDisambiguating reverse primer"
    disambiguate_primers.R \
      ${params.primer_reverse} \
      primer_R.fasta

    ## Reverse-complement primers
    echo -e "\\nReverse-complementing primers"
    seqkit seq -r -p --seq-type dna primer_F.fasta > primer_Fr.fasta
    seqkit seq -r -p --seq-type dna primer_R.fasta > primer_Rr.fasta

    """
}


// Check primers + QC + Reorient sequences
// Count number of primer occurrences withnin a read,
// discard reads with > 1 primer occurrence
// NB. read names should not contain spaces! (because of bedtools)
process primer_check {

    label "main_container"

    publishDir "${out_2_primer}", mode: "${params.storagemode}"

    // cpus 1

    // Add sample ID to the log file
    tag "${input.getSimpleName()}"

    input:
      path input
      path primer_F
      path primer_R
      path primer_Fr
      path primer_Rr

    output:
      path "${input.getSimpleName()}_PrimerChecked.fq.gz", emit: fq_primer_checked, optional: true
      path "${input.getSimpleName()}_PrimerArtefacts.fq.gz", emit: primerartefacts, optional: true
      tuple val("${task.process}"), val('seqkit'), eval('seqkit version | sed "s/seqkit v//"'), topic: versions
      tuple val("${task.process}"), val('runiq'), eval('runiq --version | sed "s/runiq //"'), topic: versions
      tuple val("${task.process}"), val('mlr'), eval('mlr --version | sed "s/mlr //"'), topic: versions
      tuple val("${task.process}"), val('bedtools'), eval('bedtools --version | sed "s/bedtools v//"'), topic: versions
      tuple val("${task.process}"), val('csvtk'), eval('csvtk version | sed "s/csvtk v//"'), topic: versions
      tuple val("${task.process}"), val('cutadapt'), eval('cutadapt --version'), topic: versions

    script:
    """
    echo -e "Input file: " ${input}
    echo -e "Forward primer: " ${params.primer_forward}
    echo -e "Reverse primer: " ${params.primer_reverse}

    ### Count number of pattern occurrences for each sequence
    count_primers (){
      # \$1 = file with primers

      seqkit replace -p "\\s.+" ${input} \
      | seqkit locate \
          --max-mismatch ${params.primer_mismatches} \
          --only-positive-strand \
          --pattern-file "\$1" \
          --threads ${task.cpus} \
      | awk -vOFS='\\t' 'NR > 1 { print \$1 , \$5 , \$6 }' \
      | runiq - \
      | mlr --tsv \
          --implicit-tsv-header \
          --headerless-tsv-output \
          sort -f 1 -n 2 \
      | bedtools merge -i stdin
    }

    echo -e "\\nCounting primers"
    echo -e "..forward primer"
    count_primers ${primer_F}  >  PF.txt

    echo -e "..rc-forward primer"
    count_primers ${primer_Fr} >> PF.txt
    
    echo -e "..reverse primer"
    count_primers ${primer_R}  >  PR.txt

    echo -e "..rc-reverse primer"
    count_primers ${primer_Rr} >> PR.txt

    ## Sort by seqID and start position, remove overlapping regions,
    ## Find duplicated records
    echo -e "\\nLooking for multiple primer occurrences"
    
    echo -e "..Processing forward primers"
    if [ -s PF.txt ]; then

      csvtk sort \
        -t -T -H -k 1:N -k 2:n \
        --num-cpus ${task.cpus} \
        PF.txt \
      | bedtools merge -i stdin \
      | awk '{ print \$1 }' \
      | runiq -i - \
      > multiprimer.txt

    else
      echo -e "...No forward primer matches found (in both orientations)"
    fi

    echo -e "..Processing reverse primers"
    if [ -s PR.txt ]; then

      csvtk sort \
        -t -T -H -k 1:N -k 2:n \
        --num-cpus ${task.cpus} \
        PR.txt \
      | bedtools merge -i stdin \
      | awk '{ print \$1 }' \
      | runiq -i - \
      >> multiprimer.txt

    else
      echo -e "...No reverse primer matches found (in both orientations)"
    fi


    ## If some artefacts are found
    if [ -s multiprimer.txt ]; then

      ## Keep only uinque seqIDs
      runiq multiprimer.txt > multiprimers.txt
      rm multiprimer.txt

      echo -e "\\nNumber of artefacts found: " \$(wc -l < multiprimers.txt)

      echo -e "..Removing artefacts"
      ## Remove primer artefacts
      seqkit grep --invert-match \
        --threads ${task.cpus} \
        --pattern-file multiprimers.txt \
        --out-file no_multiprimers.fq.gz \
        ${input}

      ## Extract primer artefacts
      echo -e "..Extracting artefacts"
      seqkit grep \
        --threads ${task.cpus} \
        --pattern-file multiprimers.txt \
        --out-file "${input.getSimpleName()}_PrimerArtefacts.fq.gz" \
        ${input}

      echo -e "..done"

    else

      echo -e "\\nNo primer artefacts found"
      ln -s ${input} no_multiprimers.fq.gz
    
    fi
    echo -e "..Done"

    echo -e "\\nReorienting sequences"

    ## Reverse-complement rev primer
    RR=\$(rc.sh ${params.primer_reverse})

    ## Reorient sequences, discard sequences without both primers
    cutadapt \
      -a ${params.primer_forward}";required;min_overlap=${params.primer_foverlap}"..."\$RR"";required;min_overlap=${params.primer_roverlap}" \
      --errors ${params.primer_mismatches} \
      --revcomp --rename "{header}" \
      --discard-untrimmed \
      --cores ${task.cpus} \
      --action none \
      --output ${input.getSimpleName()}_PrimerChecked.fq.gz \
      no_multiprimers.fq.gz

    echo -e "\\nAll done"

    ## Clean up
    if [ -f no_multiprimers.fq.gz ]; then rm no_multiprimers.fq.gz; fi

    ## Remove empty file (no valid sequences)
    echo -e "\\nRemoving empty files"
    find . -type f -name ${input.getSimpleName()}_PrimerChecked.fq.gz -size -29c -print -delete
    echo -e "..Done"
    
    """
}


// Extract ITS region with ITSx
// NB. sequence header should not contain spaces!
process itsx {

    label "main_container"

    publishDir "${out_3_itsx}", mode: "${params.storagemode}"
    // cpus 2

    // Add sample ID to the log file
    tag "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}"

    input:
      path input

    output:
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}_hash_table.txt.gz", emit: hashes, optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}_uc.uc.gz",      emit: uc,        optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.full.fasta.gz", emit: itsx_full, optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.SSU.fasta.gz",  emit: itsx_ssu,  optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.ITS1.fasta.gz", emit: itsx_its1, optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.5_8S.fasta.gz", emit: itsx_58s,  optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.ITS2.fasta.gz", emit: itsx_its2, optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.LSU.fasta.gz",  emit: itsx_lsu,  optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.positions.txt",   optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.problematic.txt", optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}_no_detections.fasta.gz", emit: itsx_nondetects, optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.summary.txt",        emit: itsx_summary, optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.extraction.results.gz", emit: itsx_details, optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.SSU.full_and_partial.fasta.gz",  emit: itsx_ssu_part,  optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.ITS1.full_and_partial.fasta.gz", emit: itsx_its1_part, optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.5_8S.full_and_partial.fasta.gz", emit: itsx_58s_part,  optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.ITS2.full_and_partial.fasta.gz", emit: itsx_its2_part, optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.LSU.full_and_partial.fasta.gz",  emit: itsx_lsu_part,  optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}_primertrimmed_sorted.fq.gz",     emit: trimmed_seqs,   optional: true
      path "parquet/*.parquet", emit: parquet, optional: true
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
    
    sampID="${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}"

    // Allow inclusion of sequences that only find a single domain, given that they meet the given E-value and score thresholds, on with parameters 1e-9,0 by default
    // singledomain = params.ITSx_singledomain ? "--allow_single_domain 1e-9,0" : ""

    """
    echo -e "Extraction of rRNA regions using ITSx\\n"

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
      ${input}

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
      > derep.fasta
    echo -e "..Done"

    ## ITSx extraction
    echo -e "\\nITSx extraction"
    ITSx \
      -i derep.fasta \
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
      -o ${sampID}
    
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

    ## Remove temporary file
    rm derep.fasta
    rm ${sampID}_primertrimmed.fq.gz

    ## Compress results
    echo -e "\\nCompressing files"

    parallel -j${task.cpus} "gzip -${params.gzip_compression} {}" ::: \
      ${sampID}_hash_table.txt \
      ${sampID}_uc.uc \
      *.fasta \
      ${sampID}.extraction.results

    ## Convert ITSx output to Parquet
    if [ ${params.ITSx_to_parquet} == true ]; then

      echo -e "\\nConverting ITSx output to Parquet"
      mkdir -p parquet

      if [ -f ${sampID}.full.fasta.gz ]; then
        ITSx_to_DuckDB.sh -i ${sampID}.full.fasta.gz -o parquet/${sampID}.full.parquet
      fi

      if [ -f ${sampID}.SSU.fasta.gz ]; then
        ITSx_to_DuckDB.sh -i ${sampID}.SSU.fasta.gz -o parquet/${sampID}.SSU.parquet
      fi

      if [ -f ${sampID}.ITS1.fasta.gz ]; then
        ITSx_to_DuckDB.sh -i ${sampID}.ITS1.fasta.gz -o parquet/${sampID}.ITS1.parquet
      fi

      if [ -f ${sampID}.5_8S.fasta.gz ]; then
        ITSx_to_DuckDB.sh -i ${sampID}.5_8S.fasta.gz -o parquet/${sampID}.5_8S.parquet
      fi

      if [ -f ${sampID}.ITS2.fasta.gz ]; then
        ITSx_to_DuckDB.sh -i ${sampID}.ITS2.fasta.gz -o parquet/${sampID}.ITS2.parquet
      fi

      if [ -f ${sampID}.LSU.fasta.gz ]; then
        ITSx_to_DuckDB.sh -i ${sampID}.LSU.fasta.gz -o parquet/${sampID}.LSU.parquet
      fi

      echo -e "Parquet files created\\n"

    fi

    echo -e "..Done"
    """
}

// Collect all ITS parts extracted by ITSx
process itsx_collect {

    label "main_container"

    publishDir "${out_3_itsxp}", mode: "${params.storagemode}"
    // cpus 1

    input:
      path(itsx_full, stageAs: "full/*")
      path(itsx_ssu,  stageAs: "ssu/*")
      path(itsx_its1, stageAs: "its1/*")
      path(itsx_58s,  stageAs: "58s/*")
      path(itsx_its2, stageAs: "its2/*")
      path(itsx_lsu,  stageAs: "lsu/*")
      path(itsx_ssu_part,  stageAs: "ssu_partial/*")
      path(itsx_its1_part, stageAs: "its1_partial/*")
      path(itsx_58s_part,  stageAs: "58s_partial/*")
      path(itsx_its2_part, stageAs: "its2_partial/*")
      path(itsx_lsu_part,  stageAs: "lsu_partial/*")

    output:
      path "ITS_Full.fasta.gz", emit: full, optional: true 
      path "SSU.fasta.gz",      emit: ssu,  optional: true 
      path "ITS1.fasta.gz",     emit: its1, optional: true 
      path "5_8S.fasta.gz",     emit: s58,  optional: true 
      path "ITS2.fasta.gz",     emit: its2, optional: true 
      path "LSU.fasta.gz",      emit: lsu,  optional: true 
      path "SSU_full_and_partial.fasta.gz",  emit: ssu_part,  optional: true 
      path "ITS1_full_and_partial.fasta.gz", emit: its1_part, optional: true 
      path "5_8S_full_and_partial.fasta.gz", emit: s58_part,  optional: true 
      path "ITS2_full_and_partial.fasta.gz", emit: its2_part, optional: true 
      path "LSU_full_and_partial.fasta.gz",  emit: lsu_part,  optional: true 
    
    script:
    """
    # Check if each sub-dir has files with rRNA regions, then concatenate

    if [[ ! -n \$(find ./full -name NOFULL) ]]; then
      echo -e "Pooling full ITS"
      find full -name "*.fasta.gz" \
        | parallel -j1 "cat {}" >> ITS_Full.fasta.gz
    fi

    if [[ ! -n \$(find ./ssu -name NOSSU) ]]; then
      echo -e "Pooling SSU"
      find ssu -name "*.fasta.gz" \
        | parallel -j1 "cat {}" >> SSU.fasta.gz
    fi

    if [[ ! -n \$(find ./its1 -name NOITS1) ]]; then
      echo -e "Pooling ITS1"
      find its1 -name "*.fasta.gz" \
        | parallel -j1 "cat {}" >> ITS1.fasta.gz
    fi
    
    if [[ ! -n \$(find ./58s -name NO58S) ]]; then
      echo -e "Pooling 5.8S"
      find 58s -name "*.fasta.gz" \
        | parallel -j1 "cat {}" >> 5_8S.fasta.gz
    fi

    if [[ ! -n \$(find ./its2 -name NOITS2) ]]; then
      echo -e "Pooling ITS2"
      find its2 -name "*.fasta.gz" \
        | parallel -j1 "cat {}" >> ITS2.fasta.gz
    fi

    if [[ ! -n \$(find ./lsu -name NOLSU) ]]; then
      echo -e "Pooling LSU"
      find lsu -name "*.fasta.gz" \
        | parallel -j1 "cat {}" >> LSU.fasta.gz
    fi

    ##### Full and partial sequences #####

    if [[ ! -n \$(find ./ssu_partial -name NOSSUPART) ]]; then
      echo -e "Pooling SSU partial sequences"
      find ssu_partial -name "*.fasta.gz" \
        | parallel -j1 "cat {}" >> SSU_full_and_partial.fasta.gz
    fi

    if [[ ! -n \$(find ./its1_partial -name NOITS1PART) ]]; then
      echo -e "Pooling ITS1 partial sequences"
      find its1_partial -name "*.fasta.gz" \
        | parallel -j1 "cat {}" >> ITS1_full_and_partial.fasta.gz
    fi

    if [[ ! -n \$(find ./58s_partial -name NO58SPART) ]]; then
      echo -e "Pooling 5.8S partial sequences"
      find 58s_partial -name "*.fasta.gz" \
        | parallel -j1 "cat {}" >> 5_8S_full_and_partial.fasta.gz
    fi

    if [[ ! -n \$(find ./its2_partial -name NOITS2PART) ]]; then
      echo -e "Pooling ITS2 partial sequences"
      find its2_partial -name "*.fasta.gz" \
        | parallel -j1 "cat {}" >> ITS2_full_and_partial.fasta.gz
    fi

    if [[ ! -n \$(find ./lsu_partial -name NOLSUPART) ]]; then
      echo -e "Pooling LSU partial sequences"
      find lsu_partial -name "*.fasta.gz" \
        | parallel -j1 "cat {}" >> LSU_full_and_partial.fasta.gz
    fi

    echo -e "\\n..Done"
    """
}


// Trim primers (do not extract ITS)
// + Estimate sequence qualities
process trim_primers {

    label "main_container"

    publishDir "${out_3_trim}", mode: "${params.storagemode}"
    // cpus 2

    // Add sample ID to the log file
    tag "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}"

    input:
      path input

    output:
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}_hash_table.txt.gz",          emit: hashes, optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}_primertrimmed_sorted.fq.gz", emit: primertrimmed_fq, optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.fa.gz",    emit: primertrimmed_fa, optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}_uc.uc.gz", emit: uc, optional: true
      tuple val("${task.process}"), val('cutadapt'), eval('cutadapt --version'), topic: versions
      tuple val("${task.process}"), val('vsearch'), eval('vsearch --version 2>&1 | head -n 1 | sed "s/vsearch //g" | sed "s/,.*//g" | sed "s/^v//" | sed "s/_.*//"'), topic: versions
      tuple val("${task.process}"), val('seqkit'), eval('seqkit version | sed "s/seqkit v//"'), topic: versions
      tuple val("${task.process}"), val('phredsort'), eval('phredsort -v | sed "s/phredsort //"'), topic: versions
      tuple val("${task.process}"), val('seqhasher'), eval('seqhasher -v | sed "s/SeqHasher //"'), topic: versions

    script:
    sampID="${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}"

    """
    echo -e "Trimming primers\\n"
    echo -e "Input sample: "   ${sampID}
    echo -e "Forward primer: " ${params.primer_forward}
    echo -e "Reverse primer: " ${params.primer_reverse}

    ## Reverse-complement rev priver
    RR=\$(rc.sh ${params.primer_reverse})
    echo -e "Reverse primer RC: " "\$RR"

    echo -e "\\nTrimming primers"
    cutadapt \
      -a ${params.primer_forward}";required;min_overlap=${params.primer_foverlap}"..."\$RR"";required;min_overlap=${params.primer_roverlap}" \
      --errors ${params.primer_mismatches} \
      --revcomp --rename "{header}" \
      --cores ${task.cpus} \
      --action=trim \
      --discard-untrimmed \
      --minimum-length ${params.trim_minlen} \
      --output ${sampID}_primertrimmed.fq.gz \
      ${input}


    if [ -n "\$(find . -name ${sampID}_primertrimmed.fq.gz -prune -size +29c)" ]; then 

      ## Estimate sequence quality and sort sequences by quality
      echo -e "\\nSorting by sequence quality"
      seqkit replace -p "\\s.+" ${sampID}_primertrimmed.fq.gz \
        | phredsort -i - -o - --metric meep --header avgphred,maxee,meep \
        | gzip -${params.gzip_compression} \
        > ${sampID}_primertrimmed_sorted.fq.gz
      echo -e "..Done"

      rm ${sampID}_primertrimmed.fq.gz

      ## Hash sequences, add sample ID to the header
      ## columns: Sample ID - Hash - PacBioID - AvgPhredScore - MaxEE - MEEP - Sequence - Quality - Length
      ## Convert to Parquet format
      echo -e "\\nCreating hash table"
      seqhasher --hash sha1 --name ${sampID} ${sampID}_primertrimmed_sorted.fq.gz - \
        | seqkit fx2tab --length \
        | sed 's/;/\t/ ; s/;/\t/ ; s/ avgphred=/\t/ ; s/ maxee=/\t/ ; s/ meep=/\t/' \
        > ${sampID}_hash_table.txt
      echo -e "..Done"

      ## Compress results
      echo -e "Compressing result"
      gzip -${params.gzip_compression} ${sampID}_hash_table.txt

      ## Dereplicate at sample level
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
        --uc ${sampID}_uc.uc \
        --quiet \
      | gzip -${params.gzip_compression} \
      > ${sampID}.fa.gz

      echo -e "..Done"

      ## Compress UC file
      gzip -${params.gzip_compression} ${sampID}_uc.uc

    else

      echo -e "\\nNo sequences found after primer removal"
      if [ -f ${sampID}_primertrimmed.fq.gz ]; then rm ${sampID}_primertrimmed.fq.gz; fi

    fi

    echo -e "..Done"

    """
}



// Assemble near-full-length ITS from ITSx output
process assemble_its {

    label "main_container"

    publishDir "${out_3_itsx}", mode: "${params.storagemode}"
    // cpus 1

    // Add sample ID to the log file
    tag "${ITS1.getSimpleName()}"

    input:
      path ITS1
      path S58
      path ITS2

    output:
      path "${ITS1.getSimpleName()}_ITS1_58S_ITS2.fasta.gz", emit: itsnf, optional: true
      // path "ITS1_58S.fasta.gz",   emit: its1p, optional: true
      // path "58S_ITS2.fasta.gz",   emit: its2p, optional: true

    script:
    sampID="${ITS1.getSimpleName()}"

    """
    echo -e "Checking if ITS1, 5.8S, and ITS2 parts are available"

    if [[ -f ${ITS1} && ${S58} && ${ITS2} ]]; then

      echo -e "\\n..All parts found"

      ## Prepare tables for ID matching
      echo -e "\\n..Converting data to tabular format"
      seqkit fx2tab ${ITS1} | sed 's/\t\$//g' | csvtk add-header -t -n id,ITS1 > tmp_1_ITS1.txt
      seqkit fx2tab ${S58}  | sed 's/\t\$//g' | csvtk add-header -t -n id,58S  > tmp_1_s58.txt
      seqkit fx2tab ${ITS2} | sed 's/\t\$//g' | csvtk add-header -t -n id,ITS2 > tmp_1_ITS2.txt

      ## Join ITS fragments
      echo -e "\\n..Joining ITS fragments"
      csvtk join -t -f   "id" tmp_1_ITS1.txt tmp_1_s58.txt tmp_1_ITS2.txt > tmp_2_ITS1_58S_ITS2.txt

      ## Check joining results
      NUMSEQS=\$(wc -l < tmp_2_ITS1_58S_ITS2.txt)
      echo "...Number of joined sequences: " \$((NUMSEQS - 1))

      if [ "\$NUMSEQS" -gt 1 ]; then

        ## Convert table back to fasta
        ## Remove leading and trailing Ns
        echo -e "\\n..Preparing fasta"
        awk 'NR>1 { print \$1 "\t" \$2\$3\$4 }' tmp_2_ITS1_58S_ITS2.txt \
          | seqkit tab2fx -w 0  \
          | seqkit replace -p "^n+|n+\$" -r "" -is -w 0 \
          | gzip -${params.gzip_compression} > ${sampID}_ITS1_58S_ITS2.fasta.gz

      else
        echo "...There are no sequences with all ITS parts present\\n"
        echo -e "\\n..Skipping ITS assembly for this sample"
      fi

    else 
      echo -e "\\n..Some or all parts are missing"
      echo -e "\\n..Skipping ITS assembly for this sample"
    fi

    """
}



// Merge tables with sequence qualities
process seq_qual {

    label "main_container"

    publishDir "${out_9_db}", mode: "${params.storagemode}"
    // cpus 4

    input:
      path(input, stageAs: "hash_tables/*")

    output:
      path "SeqQualities.parquet", emit: quals
      tuple val("${task.process}"), val('duckdb'), eval('duckdb --version | cut -d" " -f1  | sed "s/^v//"'), topic: versions

    script:
    def memoryArg = task.memory ? "-m ${task.memory.toMega()}.MB" : ""
    """
    echo -e "Aggregating sequence qualities"

    merge_hash_tables.sh \
      -i ./hash_tables \
      -o SeqQualities.parquet \
      -t ${task.cpus} \
      ${memoryArg}

    echo -e "..Done"
    """
}


// Homopolymer compression
process homopolymer {

    label "main_container"

    publishDir "${out_4_homop}", mode: "${params.storagemode}"
    // cpus 1

    // Add sample ID to the log file
    tag "${input.getSimpleName().replaceAll(/_ITS1_58S_ITS2/, '')}"

    input:
      path input

    output:
      path "${input.getSimpleName().replaceAll(/_ITS1_58S_ITS2/, '')}_Homopolymer_compressed.fa.gz", emit: hc, optional: true
      path "${input.getSimpleName().replaceAll(/_ITS1_58S_ITS2/, '')}_uch.uc.gz", emit: uch, optional: true
      tuple val("${task.process}"), val('vsearch'), eval('vsearch --version 2>&1 | head -n 1 | sed "s/vsearch //g" | sed "s/,.*//g" | sed "s/^v//" | sed "s/_.*//"'), topic: versions
      tuple val("${task.process}"), val('seqkit'), eval('seqkit version | sed "s/seqkit v//"'), topic: versions
      tuple val("${task.process}"), val('R'), eval('Rscript -e "cat(R.version.string)" | sed "s/R version //" | cut -d" " -f1'),  topic: versions
      tuple val("${task.process}"), val('data.table'), eval('Rscript -e "cat(as.character(packageVersion(\'data.table\')))"'),  topic: versions
  
    script:
    sampID="${input.getSimpleName().replaceAll(/_ITS1_58S_ITS2/, '')}"

    """

    ## Homopolymer compression
    echo -e "Homopolymer compression"

    zcat ${input} \
      | homopolymer_compression.sh - \
      > homo_compressed.fa
    
    echo -e "..Done"

    ## Re-cluster homopolymer-compressed data
    echo -e "\\nRe-clustering homopolymer-compressed data"
    vsearch \
      --cluster_size homo_compressed.fa \
      --id ${params.hp_similarity} \
      --iddef ${params.hp_iddef} \
      --qmask "dust" \
      --strand "both" \
      --fasta_width 0 \
      --threads ${task.cpus} \
      --sizein --sizeout \
      --minseqlength 20 \
      --centroids homo_clustered.fa \
      --uc ${sampID}_uch.uc
    echo -e "..Done"

    ## Check if clustering was succeful
    ## (e.g., if all compressed sequences were too short, the file with be empty)
    if [ -s homo_clustered.fa ]; then

      ## Compress UC file
      gzip -${params.gzip_compression} ${sampID}_uch.uc

      ## Substitute homopolymer-comressed sequences with uncompressed ones
      ## (update size annotaions)
      echo -e "\\nExtracting sequences"

      seqkit fx2tab ${input} > inp_tab.txt
      seqkit fx2tab homo_clustered.fa > clust_tab.txt

      if [ -s inp_tab.txt ]; then
        substitute_compressed_seqs.R \
          inp_tab.txt clust_tab.txt res.fa

        echo -e "..Done"
      else
        echo -e "..Input data looks empty, nothing to proceed with"
      fi

      if [ -s res.fa ]; then
        gzip -c res.fa > ${sampID}_Homopolymer_compressed.fa.gz
      fi

      ## Remove temporary files
      rm homo_compressed.fa
      rm homo_clustered.fa
      rm inp_tab.txt
      rm clust_tab.txt
      rm res.fa

    else
      echo -e "Clustering homopolymer-compressed sequences returned to results"
      echo -e "(most likely, sequences were too short)\\n"
    fi

    """
}


// If no homopolymer compression is required, just dereplicate the samples
process just_derep {

    label "main_container"

    // publishDir "${out_4_homop}", mode: "${params.storagemode}"
    // cpus 1

    // Add sample ID to the log file
    tag "${input.getSimpleName()}"

    input:
      path input

    output:
      path "${input.getSimpleName()}.fa.gz", emit: nhc, optional: true
      path "${input.getSimpleName()}_uch.uc.gz", emit: ucnh, optional: true
      tuple val("${task.process}"), val('vsearch'), eval('vsearch --version 2>&1 | head -n 1 | sed "s/vsearch //g" | sed "s/,.*//g" | sed "s/^v//" | sed "s/_.*//"'), topic: versions

    script:
    sampID="${input.getSimpleName()}"

    """
    echo -e "Dereplicating sequences\\n"

    vsearch \
        --derep_fulllength ${input} \
        --output - \
        --strand both \
        --fasta_width 0 \
        --threads 1 \
        --sizein --sizeout \
        --uc ${sampID}_uc.uc \
      | gzip -${params.gzip_compression} \
      > ${sampID}.fa.gz

    """
}


// Pool sequences from all samples and add sample ID into header (for OTU and "ASV" table creation)
process pool_seqs {

    label "main_container"
    
    // publishDir "${out_6_tj}", mode: "${params.storagemode}"
    // cpus 2

    input:
      path(input, stageAs: 'sequences/*')

    output:
      path "Seq_tab_not_filtered.txt.gz", emit: seqtabnf
      path "Seq_not_filtered.fa.gz",      emit: seqsnf
      tuple val("${task.process}"), val('seqkit'), eval('seqkit version | sed "s/seqkit v//"'), topic: versions
      tuple val("${task.process}"), val('parallel'), eval('parallel --version | head -n 1 | sed "s/GNU parallel //"'), topic: versions

    script:
    """

    echo -e "\\nPooling and renaming sequences"

    ## If there is a sample ID in the header already, remove it
    parallel -j 1 --group \
      --rpl '{/:} s:(.*/)?([^/.]+)(\\.[^/]+)*\$:\$2:' \
      "zcat {} \
        | sed -r '/^>/ s/;sample=[^;]*/;/g ; s/;;/;/g' \
        | sed 's/>.*/&;sample='{/:}';/ ; s/_NoChimera//g ; s/_RescuedChimera//g  ; s/_JoinedPE//g ; s/_Homopolymer_compressed//g' \
        | sed 's/Rescued_Chimeric_sequences.part_//g' \
        | sed -r '/^>/ s/;;/;/g'" \
      ::: sequences/*.fa.gz \
      | vsearch --sortbysize - --sizein --sizeout --fasta_width 0 --output - \
      | sed -r '/^>/ s/;;/;/g' \
      | gzip -${params.gzip_compression} \
      > Seq_not_filtered.fa.gz

    echo "..Done"

    echo -e "\\nExtracting sequence count table"
    seqkit seq --name Seq_not_filtered.fa.gz \
      | sed 's/;/\t/g; s/size=//; s/sample=// ; s/\t*\$//' \
      | csvtk -t cut -f 2,1,3 \
      | csvtk -t add-header -n "SampleID,SeqID,Abundance" \
      | gzip -${params.gzip_compression} \
      > Seq_tab_not_filtered.txt.gz

    echo "..Done"

    """
}


// De-novo clustering of sequences for tag-jump removal
process tj_preclust {

    label "main_container"

    // publishDir "${out_6_tj}", mode: "${params.storagemode}"
    // cpus 10

    input:
      path input

    output:
      path "TJPreclust.uc.parquet", emit: preclust_uc_parquet
      tuple val("${task.process}"), val('vsearch'), eval('vsearch --version 2>&1 | head -n 1 | sed "s/vsearch //g" | sed "s/,.*//g" | sed "s/^v//" | sed "s/_.*//"'), topic: versions

    script:
    def derep = (params.tj_id as BigDecimal).compareTo(1G) == 0    // to handle floating point comparisons too
    """
    echo -e "Pre-clustering sequences prior to tag-jump removal\\n"

    echo -e "Running dereplication\\n"
  
    vsearch \
      --derep_fulllength ${input} \
      --sizein --sizeout \
      --strand both \
      --fasta_width 0 \
      --threads 1 \
      --uc     Dereplicated.uc \
      --output Dereplicated.fa

    echo -e "\\nCompressing files"
    pigz -p ${task.cpus} -${params.gzip_compression} Dereplicated.uc
    pigz -p ${task.cpus} -${params.gzip_compression} Dereplicated.fa

    ## Additional clustering (e.g., at 99% similarity)
    if [[ ${derep} == false ]]; then

      echo -e "\\nAdditional clustering at ${params.tj_id} similarity threshold\\n"

      vsearch \
        --cluster_size Dereplicated.fa.gz \
        --id    ${params.tj_id} \
        --iddef ${params.tj_iddef} \
        --sizein --sizeout \
        --qmask dust --strand plus \
        --maxrejects 128 --maxaccepts 1 \
        --fasta_width 0 \
        --threads   ${task.cpus} \
        --uc        Clustered.uc \
        --centroids Clustered.fa

      echo -e "\\nCompressing files"
      pigz -p ${task.cpus} -${params.gzip_compression} Clustered.uc
      pigz -p ${task.cpus} -${params.gzip_compression} Clustered.fa

    fi


    ## Parse UC file
    if [[ ${derep} == true ]]; then

      echo -e "\\nParsing UC file"
      ucs --map-only --split-id --rm-dups \
        -i Dereplicated.uc.gz \
        -o TJPreclust.uc.parquet

    else

      echo -e "\\nParsing dereplicated UC file"
      ucs --map-only --split-id --rm-dups \
        -i Dereplicated.uc.gz \
        -o Dereplicated.parquet

      echo -e "\\nParsing clustered UC file"
      ucs --map-only --split-id --rm-dups \
        -i Clustered.uc.gz \
        -o Clustered.parquet

      echo -e "\\nCombining dereplication and clustering UC files"
      merge_tj_memberships.sh \
        -d Dereplicated.parquet \
        -c Clustered.parquet \
        -o TJPreclust.uc.parquet \
        -t ${task.cpus}

    fi

    echo -e "\\n..Done"
    """
}



// Tag-jump removal
process tj {

    label "main_container"

    publishDir "${out_6_tj}", mode: "${params.storagemode}"
    // cpus 1

    input:
      path seqtab  // seq table in long format
      path precls  // pre-clustered membership 

    output:
      path "Seq_tab_TagJumpFiltered.txt.gz", emit: seqtabtj
      path "TagJump_scores.qs",              emit: tjs
      path "TagJump_plot.pdf"
      tuple val("${task.process}"), val('R'), eval('Rscript -e "cat(R.version.string)" | sed "s/R version //" | cut -d" " -f1'),  topic: versions
      tuple val("${task.process}"), val('data.table'), eval('Rscript -e "cat(as.character(packageVersion(\'data.table\')))"'),  topic: versions
      tuple val("${task.process}"), val('ggplot2'), eval('Rscript -e "cat(as.character(packageVersion(\'ggplot2\')))"'),  topic: versions

    script:
    """

    echo -e "Tag-jump removal"
    
    tag_jump_removal_longtab.R \
      --seqtab ${seqtab} \
      --precls ${precls} \
      -f       ${params.tj_f} \
      -p       ${params.tj_p}

    echo "..Done"

    """
}




// Prepare a table with non-tag-jumped sequences
// Add quality estimate to singletons
// Add chimera-scores for putative de novo chimeras
process prep_seqtab {

    label "main_container"

    publishDir "${out_7_seq}", mode: "${params.storagemode}"
    // cpus 4

    input:
      path seqtab  // tag-jump filtered sequence table (long format)
      path seqsnf  // sequences in FASTA
      path denovos // de novo chimera scores
      path quals   // quality scores

    output:
      path "Seqs.parquet",      emit: seq_pq
      path "Seqs.txt.gz",       emit: seq_tl   // long table
      path "Seqs.fa.gz",        emit: seq_fa
      // path "Seqs.RData",     emit: seq_rd   // deprecated
      // path "Seq_tab.txt.gz", emit: seq_tw   // wide table
      tuple val("${task.process}"), val('R'), eval('Rscript -e "cat(R.version.string)" | sed "s/R version //" | cut -d" " -f1'),  topic: versions
      tuple val("${task.process}"), val('data.table'), eval('Rscript -e "cat(as.character(packageVersion(\'data.table\')))"'),  topic: versions
      tuple val("${task.process}"), val('arrow'), eval('Rscript -e "cat(as.character(packageVersion(\'arrow\')))"'),  topic: versions
      tuple val("${task.process}"), val('Biostrings'), eval('Rscript -e "cat(as.character(packageVersion(\'Biostrings\')))"'),  topic: versions

    script:
    """

    echo -e "Sequence table creation"
    
    seq_table_assembly.R \
      --seqtab  ${seqtab} \
      --fasta   ${seqsnf}   \
      --chimera ${denovos}  \
      --quality ${quals} \
      --threads ${task.cpus}

    echo "..Done"

    """
}







// Run summary - count number of reads in the output of different processes
process read_counts {

    label "main_container"

    publishDir "${out_8_smr}",                 mode: "${params.storagemode}", pattern: "*.xlsx"
    publishDir "${out_8_smr}/PerProcessStats", mode: "${params.storagemode}", pattern: "*.txt"
    // cpus 4

    input:
      path(input_fastq, stageAs: "1_input/*")
      path(qc, stageAs: "2_qc/*")
      path(samples_demux, stageAs: "3_demux/*")
      path(samples_primerch, stageAs: "4_primerch/*")
      path(samples_primermult, stageAs: "4_primerartefacts/*")
      path(samples_itsx_or_primertrim, stageAs: "5_itsxtrim/*")
      path(homopolymers, stageAs: "5_homopolymers/*")
      path(samples_chimref, stageAs: "6_chimref/*")
      path(samples_chimdenovo, stageAs: "7_chimdenov/*")
      path(chimera_recovered, stageAs: "8_chimrecov/*")
      path(samples_tj)
      path(seqtab)

    output:
      path "Run_summary.xlsx",                  emit: xlsx
      path "Counts_1.RawData.txt",              emit: counts_1_raw
      path "Counts_2.QC.txt",                   emit: counts_2_qc
      path "Counts_3.Demux.txt",                emit: counts_3_demux,        optional: true
      path "Counts_4.PrimerCheck.txt",          emit: counts_4_primer,       optional: true
      path "Counts_4.PrimerArtefacts.txt",      emit: counts_4_primerartef,  optional: true
      path "Counts_5.ITSx_or_PrimTrim.txt",     emit: counts_5_itsx_ptrim,   optional: true
      path "Counts_5.Homopolymers.txt",         emit: counts_5_homopolymers, optional: true
      path "Counts_6.ChimRef_reads.txt",        emit: counts_6_chimref_r,    optional: true
      path "Counts_6.ChimRef_uniqs.txt",        emit: counts_6_chimref_u,    optional: true
      path "Counts_7.ChimDenov.txt",            emit: counts_7_chimdenov,    optional: true
      path "Counts_8.ChimRecov_reads.txt",      emit: counts_8_chimrecov_r,  optional: true
      path "Counts_8.ChimRecov_uniqs.txt",      emit: counts_8_chimrecov_u,  optional: true
      tuple val("${task.process}"), val('seqkit'), eval('seqkit version | sed "s/seqkit v//"'), topic: versions
      tuple val("${task.process}"), val('parallel'), eval('parallel --version | head -n 1 | sed "s/GNU parallel //"'), topic: versions
      tuple val("${task.process}"), val('R'), eval('Rscript -e "cat(R.version.string)" | sed "s/R version //" | cut -d" " -f1'),  topic: versions
      tuple val("${task.process}"), val('data.table'), eval('Rscript -e "cat(as.character(packageVersion(\'data.table\')))"'),  topic: versions

    script:

    """
    echo -e "Summarizing run statistics\\n"
    echo -e "Counting the number of reads in:\\n"


    ## Count raw reads
    echo -e "\\n..Raw data"
    seqkit stat --basename --tabular --threads ${task.cpus} --quiet \
      1_input/* > Counts_1.RawData.txt
    
    ## Count number of reads passed QC
    echo -e "\\n..Sequenced passed QC"
    seqkit stat --basename --tabular --threads ${task.cpus} --quiet \
      2_qc/* > Counts_2.QC.txt
    
    ## Count demultiplexed reads
    echo -e "\\n..Demultiplexed data"
    seqkit stat --basename --tabular --threads ${task.cpus} --quiet \
      3_demux/* > Counts_3.Demux.txt
    

    ## Count primer-checked reads
    echo -e "\\n..Primer-checked data"
    if [ `find 4_primerch -name no_primerchecked 2>/dev/null` ]
    then
      echo -e "... No files found"
      touch Counts_4.PrimerCheck.txt
    else
      seqkit stat --basename --tabular --threads ${task.cpus} --quiet \
        4_primerch/* > Counts_4.PrimerCheck.txt
    fi


    ## Count primer-artefacts
    echo -e "\\n..Primer-artefacts"
    if [ `find 4_primerartefacts -name no_multiprimer 2>/dev/null` ]
    then
      echo -e "... No files found"
      touch Counts_4.PrimerArtefacts.txt
    else
      seqkit stat --basename --tabular --threads ${task.cpus} --quiet \
        4_primerartefacts/* > Counts_4.PrimerArtefacts.txt
    fi


    ## Count ITSx reads or primer-trimmed reads (if ITSx was not used)
    ## Take number of reads into account (--sizein)
    echo -e "\\n..ITSx- or primer-trimmed data"
    if [ `find 5_itsxtrim \\( -name no_itsx -o -name no_primertrim \\) 2>/dev/null` ]
    then
      echo -e "... No files found"
      touch Counts_5.ITSx_or_PrimTrim.txt
    else
      find 5_itsxtrim -name "*.fasta.gz" \
        | parallel -j ${task.cpus} "count_number_of_reads.sh {} {/.}" \
        | sed '1i SampleID\tNumReads' \
        > Counts_5.ITSx_or_PrimTrim.txt
    fi


    ## Count homopolymer-correction results
    echo -e "\\n..Counting homopolymer-corrected reads"
    if [ `find 5_homopolymers \\( -name no_homopolymer \\) 2>/dev/null` ]
    then
      echo -e "... No files found"
      touch Counts_5.Homopolymers.txt
    else
      find 5_homopolymers -name "*.uc.gz" \
        | parallel -j ${task.cpus} "count_homopolymer_stats.sh {} {/.}" \
        | sed '1i SampleID\tQuery\tTarget' \
        > Counts_5.Homopolymers.txt
    fi
    
    ## Count number of reads for reference-based chimeras
    echo -e "\\n..Reference-based chimeras"
    if [ `find 6_chimref -name no_chimref 2>/dev/null` ]
    then
      echo -e "... No files found"
      touch Counts_6.ChimRef_reads.txt
      touch Counts_6.ChimRef_uniqs.txt
    else
      
      ## Count number of reads
      find 6_chimref -name "*.fa.gz" \
        | parallel -j ${task.cpus} "count_number_of_reads.sh {} {/.}" \
        | sed '1i SampleID\tNumReads' \
        > Counts_6.ChimRef_reads.txt

      ## Count number of unique sequences
      seqkit stat --basename --tabular --threads ${task.cpus} --quiet \
        6_chimref/* > Counts_6.ChimRef_uniqs.txt

    fi


    ## Number of de novo chimeras (read counts are not taken into account!)
    echo -e "\\n..De novo chimeras"
    if [ `find 7_chimdenov -name no_chimdenovo 2>/dev/null` ]
    then
      echo -e "... No files found"
      touch Counts_7.ChimDenov.txt
    else
      cat 7_chimdenov/* > Counts_7.ChimDenov.txt
    fi


    ## Rescued chimeras
    echo -e "\\n..Rescued chimeric sequences"
    if [ `find 8_chimrecov -name no_chimrescued 2>/dev/null` ]
    then
      echo -e "... No files found"
      touch Counts_8.ChimRecov_reads.txt
      touch Counts_8.ChimRecov_uniqs.txt
    else
      
      ## Count number of reads
      echo -e "...Reads"
      find 8_chimrecov -name "*.fa.gz" \
        | parallel -j ${task.cpus} "count_number_of_reads.sh {} {/.}" \
        | sed '1i SampleID\tNumReads' \
        > Counts_8.ChimRecov_reads.txt

      ## Count number of unique sequences
      echo -e "...Unique sequences"
      seqkit stat --basename --tabular --threads ${task.cpus} --quiet \
        8_chimrecov/* > Counts_8.ChimRecov_uniqs.txt

    fi
    
    ## Summarize read counts
    read_count_summary.R \
      --raw          Counts_1.RawData.txt \
      --qc           Counts_2.QC.txt \
      --demuxed      Counts_3.Demux.txt \
      --primer       Counts_4.PrimerCheck.txt \
      --primerartef  Counts_4.PrimerArtefacts.txt \
      --itsx         Counts_5.ITSx_or_PrimTrim.txt \
      --homopolymer  Counts_5.Homopolymers.txt \
      --chimrefn     Counts_6.ChimRef_reads.txt \
      --chimrefu     Counts_6.ChimRef_uniqs.txt \
      --chimdenovo   Counts_7.ChimDenov.txt \
      --chimrecovn   Counts_8.ChimRecov_reads.txt \
      --chimrecovu   Counts_8.ChimRecov_uniqs.txt \
      --tj           ${samples_tj} \
      --seqtab       ${seqtab} \
      --maxchim      ${params.max_ChimeraScore} \
      --threads      ${task.cpus}

    """
}

// Quick stats of demultiplexing and primer checking steps
// (for the `seqstats` sub-workflow)
process quick_stats {

    label "main_container"

    publishDir "${out_3_quickstats}",                 mode: "${params.storagemode}", pattern: "*.xlsx"
    publishDir "${out_3_quickstats}/PerProcessStats", mode: "${params.storagemode}", pattern: "*.txt"
    // cpus 5

    input:
      path(input_fastq, stageAs: "1_input/*")
      path(qc, stageAs: "2_qc/*")
      path(samples_demux, stageAs: "3_demux/*")
      path(samples_primerch, stageAs: "4_primerch/*")
      path(samples_primermult, stageAs: "4_primerartefacts/*")

    output:
      path "Run_summary.xlsx",                  emit: xlsx
      path "Counts_1.RawData.txt",              emit: counts_1_raw
      path "Counts_2.QC.txt",                   emit: counts_2_qc
      path "Counts_3.Demux.txt",                emit: counts_3_demux,       optional: true
      path "Counts_4.PrimerCheck.txt",          emit: counts_4_primer,      optional: true
      path "Counts_4.PrimerArtefacts.txt",      emit: counts_4_primerartef, optional: true
      tuple val("${task.process}"), val('seqkit'), eval('seqkit version | sed "s/seqkit v//"'), topic: versions
      tuple val("${task.process}"), val('parallel'), eval('parallel --version | head -n 1 | sed "s/GNU parallel //"'), topic: versions
      tuple val("${task.process}"), val('R'), eval('Rscript -e "cat(R.version.string)" | sed "s/R version //" | cut -d" " -f1'),  topic: versions
      tuple val("${task.process}"), val('data.table'), eval('Rscript -e "cat(as.character(packageVersion(\'data.table\')))"'),  topic: versions

    script:

    """
    echo -e "Summarizing run statistics\\n"
    echo -e "Counting the number of reads in:\\n"


    ## Count raw reads
    echo -e "\\n..Raw data"
    seqkit stat --basename --tabular --threads ${task.cpus} \
      1_input/* > Counts_1.RawData.txt
    
    ## Count number of reads passed QC
    echo -e "\\n..Sequenced passed QC"
    seqkit stat --basename --tabular --threads ${task.cpus} \
      2_qc/* > Counts_2.QC.txt
    
    ## Count demultiplexed reads
    echo -e "\\n..Demultiplexed data"
    seqkit stat --basename --tabular --threads ${task.cpus} \
      3_demux/* > Counts_3.Demux.txt

    ## Count primer-checked reads
    echo -e "\\n..Primer-checked data"
    if [ `find 4_primerch -name no_primerchecked 2>/dev/null` ]
    then
      echo -e "... No files found"
      touch Counts_4.PrimerCheck.txt
    else
      seqkit stat --basename --tabular --threads ${task.cpus} \
        4_primerch/* > Counts_4.PrimerCheck.txt
    fi

    ## Count primer-artefacts
    echo -e "\\n..Primer-areifacts"
    if [ `find 4_primerartefacts -name no_multiprimer 2>/dev/null` ]
    then
      echo -e "... No files found"
      touch Counts_4.PrimerArtefacts.txt
    else
      seqkit stat --basename --tabular --threads ${task.cpus} \
        4_primerartefacts/* > Counts_4.PrimerArtefacts.txt
    fi
    
    ## Summarize read counts
    quick_stats.R \
      --raw          Counts_1.RawData.txt \
      --qc           Counts_2.QC.txt \
      --demuxed      Counts_3.Demux.txt \
      --primer       Counts_4.PrimerCheck.txt \
      --primerartef  Counts_4.PrimerArtefacts.txt \
      --threads      ${task.cpus}

    """
}

// Auto documentation of analysis procedures
// (generate narrative description of methods with references)
process document_analysis_s1 {

    label "main_container"

    publishDir "${out_tracedir}", mode: 'copy', overwrite: true
    // cpus 1

    input:
      path versions       // "software_versions.yml"
      path params         // "pipeline_params.tsv"

    output:
      path "README_Step1_Methods.txt",  emit: docs


    script:
    """
    echo -e "Descriptive summary generation\\n"

    document_s1.R \
      ${versions} \
      ${params} \
      README_Step1_Methods.txt

    """
}




//  The default workflow - Step-1
workflow S1 {

  // Primer disambiguation
  disambiguate()


  /*
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Demultiplex data
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  */

  // Run demultiplexing
  if( params.demultiplexed == false ){
    
    // Input file with barcodes (FASTA)
    ch_barcodes = Channel.value(params.barcodes)

    // Validate tags
    tag_validation(ch_barcodes)

    // PacBio
    if ( params.seqplatform == "PacBio" ) {
      
      // Input file with multiplexed reads (FASTQ.gz or BAM)
      ch_input = Channel.value(params.input)

      // Check the extension of input
      input_type = file(params.input).getExtension() =~ /bam|BAM/ ? "bam" : "oth"
      // println("${input_type}")

      // If BAM is provided as input, convert it to FASTQ
      if ( input_type == 'bam'){

        // Add BAM index file
        ch_input_pbi = ch_input + ".pbi"

        bam2fastq(ch_input, ch_input_pbi)
        qc_se(bam2fastq.out.fastq)

      } else {

        // Initial QC
        qc_se(ch_input)

      }

      // Demultiplexing with dual barcodes requires 4 additional files:
      //  - "biosamples" with symmertic/asymmetirc tag combinations
      //  - table for assigning sample names to demuxed files
      //  - and a table for renaming unknown combinations (if params.lima_remove_unknown == true)
      // Create dummy files (for single or symmetic tags) if neccesary
      ch_biosamples_sym  = tag_validation.out.biosamples_sym.flatten().collect().ifEmpty(file("biosamples_sym"))
      ch_biosamples_asym = tag_validation.out.biosamples_asym.flatten().collect().ifEmpty(file("biosamples_asym"))
      ch_file_renaming   = tag_validation.out.file_renaming.flatten().collect().ifEmpty(file("file_renaming"))
      ch_unknown_combs   = tag_validation.out.unknown_combinations.flatten().collect().ifEmpty(file("unknown_combinations"))

      // Demultiplexing
      demux(
        qc_se.out.filtered,
        tag_validation.out.fasta,
        ch_biosamples_sym, 
        ch_biosamples_asym,
        ch_file_renaming,
        ch_unknown_combs)

      // Check primers
      primer_check(
        demux.out.samples_demux.flatten(),
        disambiguate.out.F,
        disambiguate.out.R,
        disambiguate.out.Fr,
        disambiguate.out.Rr
        )

    } // end of PacBio-specific tasks

    // Illumina 
    if ( params.seqplatform == "Illumina" ) {
      
      // Input file with multiplexed pair-end reads (FASTQ.gz)
      ch_inputR1 = Channel.value(params.input_R1)
      ch_inputR2 = Channel.value(params.input_R2)

      // Initial QC
      qc_pe(ch_inputR1, ch_inputR2)

      // PE assembly
      merge_pe(
        qc_pe.out.filtered_R1,
        qc_pe.out.filtered_R2)

      // Modify barcodes (restict search window)
      prep_barcodes(tag_validation.out.fasta)

      // Demultiplexing
      demux_illumina(
        merge_pe.out.r12,
        prep_barcodes.out.barcodesm)

      ch_demux_merged = demux_illumina.out.samples_demux.flatten()

      // Illumina nonmerged PE reads sub-workflow (optional)
      if(params.illumina_keep_notmerged == true){

        // Demultiplexing non-merged reads
        demux_pe(
          merge_pe.out.nm,
          prep_barcodes.out.barcodesm)

        // Channel of non-merged reads by sample (split into sample tuples)
        // ch_R1 = demux_pe.out.demux_pe....

        // Non-merged sample list
        ch_nonmerged = demux_pe.out.samples_nonm_pe.splitText().map{it -> it.trim()}

        // Trim primers of nonmerged PE reads
        // Estimate sequence qualities
        // Dereplicate R1 and R2 independently
        // trim_primers_pe(demux_pe.out.demux_pe.flatten())

        // Join nonmerged reads with poly-N pads
        join_pe(
          ch_nonmerged,
          demux_pe.out.demux_pe.flatten().collect()   // all non-merged R1 and R2 files
          )

        // Add joined reads to the merged reads
        ch_joined = join_pe.out.jj_FQ.flatten()
        ch_demuxed = ch_demux_merged.concat(ch_joined)

      } else { // end of Illumina non-merged reads

        // Channel with demultiplexed reads
        ch_demuxed = ch_demux_merged

      }

      // Check primers
      primer_check(
        ch_demuxed,
        disambiguate.out.F,
        disambiguate.out.R,
        disambiguate.out.Fr,
        disambiguate.out.Rr
        )

    } // end of Illumina-specific tasks

  }   // end of demultiplexing



  // If samples were already demuliplexed
  if( params.demultiplexed == true ){

    // Input files with demultiplexed reads (FASTQ.gz)
    ch_input = Channel.fromPath( params.input + '/*.{fastq.gz,fastq,fq.gz,fq}' )

    // Check if the input channel is empty
    ch_input
      .ifEmpty {
          error("ERROR: No FASTQ files found in the input directory: ${params.input}")
          exit 1
      }   

    // QC
    qc_se(ch_input)

    // Check primers
    primer_check(
      qc_se.out.filtered,
      disambiguate.out.F,
      disambiguate.out.R,
      disambiguate.out.Fr,
      disambiguate.out.Rr
      )

  }  // end of pre-demultiplexed branch


  /*
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ITS extraction or primer trimming
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  */

  // Extract ITS
  if(params.its_region == "full" || params.its_region == "ITS1" || params.its_region == "ITS2" || params.its_region == "SSU" || params.its_region == "LSU"){

    // Run ITSx
    itsx(primer_check.out.fq_primer_checked)

    // Merge tables with sequence qualities
    seq_qual(itsx.out.hashes.collect())
  }

  // Trim the primers (instead of ITS extraction)
  if(params.its_region == "none"){
      
    // Trim primers with cutadapt
    trim_primers(primer_check.out.fq_primer_checked)

    // Merge tables with sequence qualities
    seq_qual(trim_primers.out.hashes.collect())
  }

  // Trim the primers, run ITSx, and assemble near-full-length ITS
  if(params.its_region == "ITS1_5.8S_ITS2"){

    // Run ITSx
    itsx(primer_check.out.fq_primer_checked)

    // Assemble ITS1-5.8S-ITS2 from ITSx-extracted parts
    if (params.ITSx_partial == 0) {
      assemble_its(
        itsx.out.itsx_its1,
        itsx.out.itsx_58s,
        itsx.out.itsx_its2)
    } else {
      assemble_its(
        itsx.out.itsx_its1_part,
        itsx.out.itsx_58s,
        itsx.out.itsx_its2_part)
    }

    // Merge tables with sequence qualities
    seq_qual(itsx.out.hashes.collect())
  }


  // Collect ITSx-extracted sequences
  if(params.its_region == "full" || params.its_region == "ITS1" || params.its_region == "ITS2" || params.its_region == "SSU" || params.its_region == "LSU" || params.its_region == "ITS1_5.8S_ITS2"){

    // Collect rRNA parts into separate channels
    ch_cc_full = itsx.out.itsx_full.flatten().collect().ifEmpty(file("NOFULL"))
    ch_cc_ssu  = itsx.out.itsx_ssu.flatten().collect().ifEmpty(file("NOSSU"))
    ch_cc_its1 = itsx.out.itsx_its1.flatten().collect().ifEmpty(file("NOITS1"))
    ch_cc_58s  = itsx.out.itsx_58s.flatten().collect().ifEmpty(file("NO58S"))
    ch_cc_its2 = itsx.out.itsx_its2.flatten().collect().ifEmpty(file("NOITS2"))
    ch_cc_lsu  = itsx.out.itsx_lsu.flatten().collect().ifEmpty(file("NOLSU"))
    
    ch_cc_ssu_part  = itsx.out.itsx_ssu_part.flatten().collect().ifEmpty(file("NOSSUPART"))
    ch_cc_its1_part = itsx.out.itsx_its1_part.flatten().collect().ifEmpty(file("NOITS1PART"))
    ch_cc_58s_part  = itsx.out.itsx_58s_part.flatten().collect().ifEmpty(file("NO58SPART"))
    ch_cc_its2_part = itsx.out.itsx_its2_part.flatten().collect().ifEmpty(file("NOITS2PART"))
    ch_cc_lsu_part  = itsx.out.itsx_lsu_part.flatten().collect().ifEmpty(file("NOLSUPART"))

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

  } // end of ITSx-extracted sequences
  

  /*
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Homopolymer compression & chimera removal
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  */

  // Homopolymer compression
  if(params.hp == true){

    // --Full-length ITS sequences
    if(params.its_region == "full"){
      homopolymer(itsx.out.itsx_full)
    }
    // --ITS1 sequences
    if(params.its_region == "ITS1"){
      if (params.ITSx_partial == 0) {
        homopolymer(itsx.out.itsx_its1)
      } else {
        homopolymer(itsx.out.itsx_its1_part)
      }
    }
    // --ITS2 sequences
    if(params.its_region == "ITS2"){
      if (params.ITSx_partial == 0) {
        homopolymer(itsx.out.itsx_its2)
      } else {
        homopolymer(itsx.out.itsx_its2_part)
      }
    }
    // --SSU sequences
    if(params.its_region == "SSU"){
      if (params.ITSx_partial == 0) {
        homopolymer(itsx.out.itsx_ssu)
      } else {
        homopolymer(itsx.out.itsx_ssu_part)
      }
    }
    // --LSU sequences
    if(params.its_region == "LSU"){
      if (params.ITSx_partial == 0) {
        homopolymer(itsx.out.itsx_lsu)
      } else {
        homopolymer(itsx.out.itsx_lsu_part)
      }
    }

    // --Primer-trimmed sequences
    if(params.its_region == "none"){
      homopolymer(trim_primers.out.primertrimmed_fa)
    }
    // Near-full-length ITS
    if(params.its_region == "ITS1_5.8S_ITS2"){
      homopolymer(assemble_its.out.itsnf)
    }


  } else {
  // No homopolymer comression is required,
  // Just dereplicate the data

    if(params.its_region == "full" || params.its_region == "ITS1" || params.its_region == "ITS2" || params.its_region == "SSU" || params.its_region == "LSU"){
      
      // --Full-length ITS sequences
      if(params.its_region == "full"){
        just_derep(itsx.out.itsx_full)
      }
      // --ITS1 sequences
      if(params.its_region == "ITS1"){
        if (params.ITSx_partial == 0) {
          just_derep(itsx.out.itsx_its1)
        } else {
          just_derep(itsx.out.itsx_its1_part)
        }
      }
      // --ITS2 sequences
      if(params.its_region == "ITS2"){
        if (params.ITSx_partial == 0) {
          just_derep(itsx.out.itsx_its2)
        } else {
          just_derep(itsx.out.itsx_its2_part)
        }
      }
      // --SSU sequences
      if(params.its_region == "SSU"){
        if (params.ITSx_partial == 0) {
          just_derep(itsx.out.itsx_ssu)
        } else {
          just_derep(itsx.out.itsx_ssu_part)
        }
      }
      // --LSU sequences
      if(params.its_region == "LSU"){
        if (params.ITSx_partial == 0) {
          just_derep(itsx.out.itsx_lsu)
        } else {
          just_derep(itsx.out.itsx_lsu_part)
        }
      }
    
    }  // end of ITS


  } // end of homopolymer correction condition


  /*
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Chimera removal
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  */


  // Chimera removal (optional)
  ch_chimerabd = Channel.value(params.chimera_db)

  // Input depends on the selected workflow
  if(params.hp == true){

    ch_input_for_chim = homopolymer.out.hc

  } else {

    if(params.its_region == "none"){
      ch_input_for_chim = trim_primers.out.primertrimmed_fa
    } else if(params.its_region == "ITS1_5.8S_ITS2"){
      ch_input_for_chim = assemble_its.out.itsnf
    } else {
      ch_input_for_chim = just_derep.out.nhc
    }
    
  }
  
  CHIMERA_REMOVAL(ch_input_for_chim, ch_chimerabd)



  /*
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Data aggregation
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  */

  // Pool sequences (for a final sequence table)
  pool_seqs(CHIMERA_REMOVAL.out.filtered)

  // Tag-jump removal
  if(params.tj == true){

    // Pre-clustering prior to tag-jump removal
    tj_preclust(pool_seqs.out.seqsnf)

    // Tag-jump removal
    tj(
      pool_seqs.out.seqtabnf,
      tj_preclust.out.preclust_uc_parquet)

    ch_seqtab_after_tj = tj.out.seqtabtj
    ch_tj_scores       = tj.out.tjs

  } else {

    // Skip tag-jump removal
    ch_seqtab_after_tj = pool_seqs.out.seqtabnf
    ch_tj_scores = file("no_tj")

  }

  // Check optional channel with de novo chimera scores
  ch_denovoscores = CHIMERA_REMOVAL.out.denovo_agg.ifEmpty(file('DeNovo_Chimera.txt'))

  // Create sequence table
  prep_seqtab(
    ch_seqtab_after_tj,    // (optionally) tag-jump-filtered sequence table (long format)
    pool_seqs.out.seqsnf,  // Sequences in FASTA format
    ch_denovoscores,       // de novo chimera scores
    seq_qual.out.quals     // sequence qualities
    )



  /*
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Read count summary
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  */
 
  // Initial data - Per-sample input channels
  if( params.demultiplexed == false ){

    if(params.seqplatform == "PacBio"){

      // Input data and QC = single multiplexed file
      ch_counts_1 = ch_input
      ch_counts_2 = qc_se.out.filtered

      ch_all_demux = demux.out.samples_demux.flatten().collect()
    }

    if(params.seqplatform == "Illumina"){
      ch_all_demux = demux_illumina.out.samples_demux.flatten().collect()
    }

  } else {
  
    // Input data and QC = several demultiplexed files
    ch_counts_1 = ch_input.flatten().collect()
    ch_counts_2 = qc_se.out.filtered.flatten().collect()

    ch_all_demux = Channel.fromPath( params.input + '/*.{fastq.gz,fastq,fq.gz,fq}' ).flatten().collect()
  }
  

  // Primer-checked and multiprimer sequences
  ch_all_primerchecked = primer_check.out.fq_primer_checked.flatten().collect().ifEmpty(file("no_primerchecked"))
  ch_all_primerartefacts = primer_check.out.primerartefacts.flatten().collect().ifEmpty(file("no_multiprimer"))
      
  // ITSx and primer trimming channel
  if(params.its_region == "full"){
    ch_all_trim = itsx.out.itsx_full.flatten().collect().ifEmpty(file("no_itsx"))
  }
  if(params.its_region == "ITS1"){
    if (params.ITSx_partial == 0) {
      ch_all_trim = itsx.out.itsx_its1.flatten().collect().ifEmpty(file("no_itsx"))
    } else {
      ch_all_trim = itsx.out.itsx_its1_part.flatten().collect().ifEmpty(file("no_itsx"))
    }
  }
  if(params.its_region == "ITS2"){
    if (params.ITSx_partial == 0) {
      ch_all_trim = itsx.out.itsx_its2.flatten().collect().ifEmpty(file("no_itsx"))
    } else {
      ch_all_trim = itsx.out.itsx_its2_part.flatten().collect().ifEmpty(file("no_itsx"))
    }
  }
  if(params.its_region == "SSU"){
    if (params.ITSx_partial == 0) {
      ch_all_trim = itsx.out.itsx_ssu.flatten().collect().ifEmpty(file("no_itsx"))
    } else {
      ch_all_trim = itsx.out.itsx_ssu_part.flatten().collect().ifEmpty(file("no_itsx"))
    }
  }
  if(params.its_region == "LSU"){
    if (params.ITSx_partial == 0) {
      ch_all_trim = itsx.out.itsx_lsu.flatten().collect().ifEmpty(file("no_itsx"))
    } else {
      ch_all_trim = itsx.out.itsx_lsu_part.flatten().collect().ifEmpty(file("no_itsx"))
    }
  }
  if(params.its_region == "ITS1_5.8S_ITS2"){
    ch_all_trim = assemble_its.out.itsnf.flatten().collect().ifEmpty(file("no_itsx"))
  }
  if(params.its_region == "none"){
    ch_all_trim = trim_primers.out.primertrimmed_fq.flatten().collect().ifEmpty(file("no_primertrim"))
  }

  // Homopolymer-correction channel
  if(params.hp == true){
    ch_homopolymers = homopolymer.out.uch.flatten().collect().ifEmpty(file("no_homopolymer"))
  } else {
    ch_homopolymers = file("no_homopolymer")
  }

  // Chimeric channels
  ch_chimref     = CHIMERA_REMOVAL.out.chimeric.flatten().collect().ifEmpty(file("no_chimref"))
  ch_chimdenovo  = CHIMERA_REMOVAL.out.denovo_agg.flatten().collect().ifEmpty(file("no_chimdenovo"))
  ch_chimrescued = CHIMERA_REMOVAL.out.rescued.flatten().collect().ifEmpty(file("no_chimrescued"))

  // Count reads and prepare summary stats for the run
  // Currently, implemented only for PacBio
  // For Illumina, need replace:
  //   `ch_input` -> `ch_inputR1` & `ch_inputR2`
  //   `qc_se`    -> `qc_pe`

  if(params.seqplatform == "PacBio"){

    read_counts(
      ch_counts_1,             // input data (single multiplexed file or several demultiplexed files)
      ch_counts_2,             // data that passed QC (single or several demuxed files)
      ch_all_demux,            // demultiplexed sequences per sample
      ch_all_primerchecked,    // primer-cheched sequences
      ch_all_primerartefacts,  // multiprimer artefacts
      ch_all_trim,             // ITSx-extracted or primer-trimmed sequences
      ch_homopolymers,         // Homopolymer stats
      ch_chimref,              // Reference-based chimeras
      ch_chimdenovo,           // De novo chimeras
      ch_chimrescued,          // Rescued chimeras
      ch_tj_scores,            // Tag-jump filtering scores
      prep_seqtab.out.seq_pq   // Final table with sequences (in Parquet format)
      )

  } // end of read_counts for PacBio


  
  // Dump the software versions to a file
  ch_versions_yml = software_versions_to_yaml(Channel.topic('versions'))
      .collectFile(
          storeDir: "${params.tracedir}",
          name:     'software_versions.yml',
          sort:     true,
          newLine:  true
      )

  // Dump the parameters to a file
  ch_params_tsv = dumpParamsTsv()
    .collectFile(
        storeDir: "${params.tracedir}",
        name:     "pipeline_params.tsv",
        sort:     true,
        newLine:  true
    )

  // Document the analysis procedures
  document_analysis_s1(
    ch_versions_yml,
    ch_params_tsv)

}






// Quick workflow for demultiplexing and estimation of the number of reads per sample
// Only PacBio non-demultiplexed reads are supported
workflow seqstats {

  // Primer disambiguation
  disambiguate()

  // Input file with barcodes (FASTA)
  ch_barcodes = Channel.value(params.barcodes)

  // Input file with multiplexed reads (FASTQ.gz)
  ch_input = Channel.value(params.input)

  // Initial QC
  qc_se(ch_input)

  // Validate tags
  tag_validation(ch_barcodes)

  // Tag-validation channels
  ch_biosamples_sym  = tag_validation.out.biosamples_sym.flatten().collect().ifEmpty(file("biosamples_sym"))
  ch_biosamples_asym = tag_validation.out.biosamples_asym.flatten().collect().ifEmpty(file("biosamples_asym"))
  ch_file_renaming   = tag_validation.out.file_renaming.flatten().collect().ifEmpty(file("file_renaming"))
  ch_unknown_combs   = tag_validation.out.unknown_combinations.flatten().collect().ifEmpty(file("unknown_combinations"))

  // Demultiplexing
  demux(
    qc_se.out.filtered,
    tag_validation.out.fasta,
    ch_biosamples_sym, 
    ch_biosamples_asym,
    ch_file_renaming,
    ch_unknown_combs)

  // Check primers
  primer_check(
    demux.out.samples_demux.flatten(),
    disambiguate.out.F,
    disambiguate.out.R,
    disambiguate.out.Fr,
    disambiguate.out.Rr
    )

  // Prepare input channels
  ch_all_demux = demux.out.samples_demux.flatten().collect()
  ch_all_primerchecked = primer_check.out.fq_primer_checked.flatten().collect().ifEmpty(file("no_primerchecked"))
  ch_all_primerartefacts = primer_check.out.primerartefacts.flatten().collect().ifEmpty(file("no_multiprimer"))

  // Count reads and prepare summary stats for the run
  quick_stats(
      ch_input,                // input data
      qc_se.out.filtered,      // data that passed QC
      ch_all_demux,            // demultiplexed sequences per sample
      ch_all_primerchecked,    // primer-cheched sequences
      ch_all_primerartefacts   // primer artefacts
      )

} // end of `seqstats` subworkflow

