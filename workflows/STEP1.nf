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
include { ITS_EXTRACTION }            from '../subworkflows/itsx_subworkflow.nf'

// Illumina paired-end reads (demultiplexing, reorientation, read merging)
include { ILLUMINA_PE }               from '../subworkflows/illumina_subworkflow.nf'


// Convert BAM to FASTQ
process bam2fastq {

    label "main_container"
    publishDir "${params.outdir}/00_BAM2FASTQ", mode: "${params.storagemode}"

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

    bam2fastq \\
      -c ${params.gzip_compression} \\
      --num-threads ${task.cpus} \\
      ${input}

    echo -e "\\nConvertion finished"
    """
}


// Quality filtering for single-end reads
process qc_se {

    label "main_container"

    // cpus 8

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
    echo -e "Homopolymer length to remove: " ${params.qc_maxhomopolymerlen}

    ## Split the number of CPUs into two (for seqkit and pigz)
    ## NB! vsearch currently does not suppot multithreading for `--fastq_filter`
    ## see https://github.com/torognes/vsearch/issues/466
    total_cpus=${task.cpus}
    half_cpus=\$(( (total_cpus + 1) / 2 ))
    (( half_cpus < 1 )) && half_cpus=1
    (( half_cpus > total_cpus )) && half_cpus=\$total_cpus

    ## Instead of using regex (e.g., "(A{25,}|C{25,}|T{25,}|G{25,})"), create fixed patterns
    A_run=\$(printf '%*s' "${params.qc_maxhomopolymerlen}" '' | tr ' ' 'A')
    C_run=\$(printf '%*s' "${params.qc_maxhomopolymerlen}" '' | tr ' ' 'C')
    G_run=\$(printf '%*s' "${params.qc_maxhomopolymerlen}" '' | tr ' ' 'G')
    T_run=\$(printf '%*s' "${params.qc_maxhomopolymerlen}" '' | tr ' ' 'T')

    ## We do not need to change the file name (output name should be the same as input)
    ## Therefore, temporary rename input
    mv ${input} inp.fq.gz

    vsearch \\
      --fastq_filter inp.fq.gz \\
      --fastq_qmax 93 \\
      ${filter_maxee} \\
      ${filter_maxeerate} \\
      --fastq_maxns ${params.qc_maxn} \\
      --threads 1 \\
      --no_progress \\
      --fastqout - \\
    | seqkit grep \\
      --by-seq --ignore-case --invert-match --only-positive-strand -w 0 \\
      --threads \$half_cpus \\
      --pattern "\$A_run" \\
      --pattern "\$C_run" \\
      --pattern "\$G_run" \\
      --pattern "\$T_run" \\
    | pigz -p \$half_cpus -${params.gzip_compression} \\
    > "${input.getSimpleName()}.fq.gz"

    echo -e "\\nQC finished"
    """
}


// Validate tags for demultiplexing
process tag_validation {

    label "main_container"
    // cpus 1

    publishDir "${params.outdir}/01_Demux", pattern: "tag_names_renamed.tsv", mode: "${params.storagemode}"

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
    validate_tags.R \\
      --tags   ${barcodes} \\
      --output barcodes_validated.fasta

    echo -e "\\nTag validation finished"
    """
}



// Demultiplexing with LIMA - for PacBio reads
process demux {

    label "main_container"

    publishDir "${params.outdir}/01_Demux", mode: "${params.storagemode}"  // , saveAs: { filename -> "foo_$filename" }
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


// Primer disambiguation
process disambiguate {

    label "main_container"

    // publishDir "${params.outdir}/02_PrimerCheck", mode: "${params.storagemode}"
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

    publishDir "${params.outdir}/02_PrimerCheck", mode: "${params.storagemode}"

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



// Merge tables with sequence qualities
process seq_qual {

    label "main_container"

    publishDir "${params.outdir}/09_DB", mode: "${params.storagemode}"
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

    publishDir "${params.outdir}/04_Homopolymer", mode: "${params.storagemode}"
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

    // publishDir "${params.outdir}/04_Homopolymer", mode: "${params.storagemode}"
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
    
    // publishDir "${params.outdir}/06_TagJumpFiltration", mode: "${params.storagemode}"
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

    // publishDir "${params.outdir}/06_TagJumpFiltration", mode: "${params.storagemode}"
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

    publishDir "${params.outdir}/06_TagJumpFiltration", mode: "${params.storagemode}"
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

    publishDir "${params.outdir}/07_SeqTable", mode: "${params.storagemode}"
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

    publishDir "${params.outdir}/08_RunSummary",                 mode: "${params.storagemode}", pattern: "*.{xlsx,tsv}"
    publishDir "${params.outdir}/08_RunSummary/PerProcessStats", mode: "${params.storagemode}", pattern: "*.txt"
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
      path "per_sample.tsv",                    emit: per_sample
      path "per_run.tsv",                       emit: per_run
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
      find 5_itsxtrim \\( -name "*.fasta.gz" -o -name "*.fa.gz" \\) \
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

    publishDir "${params.outdir}/03_Stats",                 mode: "${params.storagemode}", pattern: "*.xlsx"
    publishDir "${params.outdir}/03_Stats/PerProcessStats", mode: "${params.storagemode}", pattern: "*.txt"
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

    publishDir "${params.tracedir}", mode: 'copy', overwrite: true
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

  is_demultiplexed = params.demultiplexed
  is_illumina = params.seqplatform == "Illumina"
  run_hp = params.hp
  run_tj = params.tj

  // Primer disambiguation
  disambiguate()


  /*
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Demultiplex data
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  */

  // Illumina paired-end reads (multiplexed or per-sample)
  if( is_illumina ){

    if( !is_demultiplexed ){

      // Validate tags
      tag_validation(channel.value(params.barcodes))

      // Multiplexed read pairs
      ch_illumina_multiplexed = channel.of( tuple(file(params.input_R1), file(params.input_R2)) )
      ch_illumina_persample   = channel.empty()

      // Tags: single or symmetric dual tags (FASTA), dual tags (`tags_fwd.fasta` + `tags_rev.fasta`)
      ch_illumina_tags      = tag_validation.out.fasta
      ch_illumina_tags_dual = tag_validation.out.tags_dual.ifEmpty(file("no_dual_tags"))

    } else {

      // Per-sample read pairs, tuple(sampleID, [R1, R2])
      // Illumina-style suffixes (`_S1_L001`) are removed from sample names
      ch_illumina_multiplexed = channel.empty()
      ch_illumina_persample   = channel
        .fromFilePairs( params.input + '/' + params.illumina_pe_pattern, size: 2 )
        .ifEmpty { error("ERROR: No paired FASTQ files matching `${params.illumina_pe_pattern}` found in the input directory: ${params.input}") }
        .map { id, reads ->
          def sampID = id.replaceAll(/_S\d+_L\d{3}$/, '')
          if( sampID.contains('.') ){
            error("ERROR: sample name `${sampID}` contains a dot, please rename the input files")
          }
          tuple(sampID, reads)
        }

      ch_illumina_tags      = file("no_tags")
      ch_illumina_tags_dual = file("no_dual_tags")
    }

    // Demultiplexing, reorientation, read merging (+ optional joining of non-merged reads)
    ILLUMINA_PE(
      ch_illumina_multiplexed,
      ch_illumina_persample,
      ch_illumina_tags,
      ch_illumina_tags_dual)

    // QC of merged reads
    qc_se(ILLUMINA_PE.out.merged)

    // Channel to use for primer checking
    // (joined reads are quality-filtered prior to joining)
    ch_for_primer_check = qc_se.out.filtered.mix(ILLUMINA_PE.out.joined)



  } else if( !is_demultiplexed ){
  // PacBio multiplexed reads

    // Input file with barcodes (FASTA)
    ch_barcodes = channel.value(params.barcodes)

    // Validate tags
    tag_validation(ch_barcodes)

    // Input file with multiplexed reads (FASTQ.gz or BAM)
    ch_input = channel.value(params.input)

    // Check the extension of input
    input_type = file(params.input).getExtension() =~ /bam|BAM/ ? "bam" : "oth"

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

    // Channel to use for primer checking
    ch_for_primer_check = demux.out.samples_demux.flatten()

  } else {
  // If samples were already demuliplexed (single-end reads, any platform)

    // Input files with demultiplexed reads (FASTQ.gz)
    ch_input = channel.fromPath( params.input + '/*.{fastq.gz,fastq,fq.gz,fq}' )

    // Check if the input channel is empty
    ch_input
      .ifEmpty {
          error("ERROR: No FASTQ files found in the input directory: ${params.input}")
          exit(1)
      }   

    // QC
    qc_se(ch_input)

    // Channel to use for primer checking
    ch_for_primer_check = qc_se.out.filtered

  }  // end of pre-demultiplexed branch

  // Check primers
  primer_check_out = primer_check(
    ch_for_primer_check,
    disambiguate.out.F,
    disambiguate.out.R,
    disambiguate.out.Fr,
    disambiguate.out.Rr
    )


  /*
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ITS extraction or primer trimming
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  */

  // Trim primers, dereplicate, and (optionally) extract the target rRNA region
  //   `params.its_region` selects the region ("none" = primer trimming only)
  //   `params.itsx_tool`  selects the extractor (ITSx v1.x or ITSx2)
  ITS_EXTRACTION(primer_check_out.fq_primer_checked)

  // Merge tables with sequence qualities
  seq_qual(ITS_EXTRACTION.out.hashes.collect())
  

  /*
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Homopolymer compression & chimera removal
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  */

  // The sequences to work with are already selected by the subworkflow,
  // according to `params.its_region` (and `params.ITSx_partial` for ITSx v1.x)
  ch_region_seqs = ITS_EXTRACTION.out.region_seqs

  // "none" data is already dereplicated by `primer_trim`, and the near-full-length ITS
  // assembled by `get_its` inherits the dereplicated IDs - only the extracted regions
  // need to be re-dereplicated (extraction can collapse distinct amplicons into
  // identical rRNA regions)
  run_derep = !(params.its_region in ["none", "ITS1_5.8S_ITS2"])

  // Homopolymer compression
  if(run_hp){

    homopolymer(ch_region_seqs)

  } else if(run_derep){

    // No homopolymer compression is required - just dereplicate the data
    just_derep(ch_region_seqs)

  } // end of homopolymer correction condition


  /*
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Chimera removal
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  */


  // Chimera removal (optional)
  ch_chimerabd = channel.value(params.chimera_db)

  // Input depends on the selected workflow
  if(run_hp){
    ch_input_for_chim = homopolymer.out.hc
  } else if(run_derep){
    ch_input_for_chim = just_derep.out.nhc
  } else {
    ch_input_for_chim = ch_region_seqs
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
  if(run_tj){

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
 

  if( is_illumina ){

    // Raw data = read pairs (R1 only, to count pairs)
    ch_counts_1 = is_demultiplexed
      ? ch_illumina_persample.map { _id, reads -> reads[0] }.collect()
      : file(params.input_R1)

    // QC = per-sample merged reads that passed QC
    ch_counts_2 = qc_se.out.filtered.flatten().collect().ifEmpty(file("no_qc"))

    // Per-sample demultiplexed pairs are taken from the reorientation stats
    ch_all_demux = file("no_demux")

    // Illumina-specific stats
    ch_pe_reorient = ILLUMINA_PE.out.reorient_stats
    ch_pe_merge    = ILLUMINA_PE.out.merge_stats
    ch_pe_joined   = ILLUMINA_PE.out.joined.flatten().collect().ifEmpty(file("no_joined"))

  } else if( !is_demultiplexed ){

    // Input data and QC = single multiplexed file
    ch_counts_1 = ch_input
    ch_counts_2 = qc_se.out.filtered

    ch_all_demux = demux.out.samples_demux.flatten().collect()

  } else {
  
    // Input data and QC = several demultiplexed files
    ch_counts_1 = ch_input.flatten().collect()
    ch_counts_2 = qc_se.out.filtered.flatten().collect()

    ch_all_demux = channel.fromPath( params.input + '/*.{fastq.gz,fastq,fq.gz,fq}' ).flatten().collect()
  }

  if( !is_illumina ){
    ch_pe_reorient = file("no_reorient_stats")
    ch_pe_merge    = file("no_merge_stats")
    ch_pe_joined   = file("no_joined")
  }
  

  // Primer-checked and multiprimer sequences
  ch_all_primerchecked = primer_check_out.fq_primer_checked.flatten().collect().ifEmpty(file("no_primerchecked"))
  ch_all_primerartefacts = primer_check_out.primerartefacts.flatten().collect().ifEmpty(file("no_multiprimer"))
      
  // Did the ITS extraction run? Used below for the report inputs
  run_itsx = params.its_region != "none"

  // ITSx-extracted or primer-trimmed sequences (for the read count summary)
  // NB. these counts are read from the `;size=` annotations of the dereplicated sequences, 
  //     so the FASTA is used here and not the quality-sorted FASTQ
  ch_all_trim = run_itsx
    ? ITS_EXTRACTION.out.region_seqs.flatten().collect().ifEmpty(file("no_itsx"))
    : ITS_EXTRACTION.out.region_seqs.flatten().collect().ifEmpty(file("no_primertrim"))

  // Homopolymer-correction channel
  if(run_hp){
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
  ch_versions_yml = software_versions_to_yaml(channel.topic('versions'))
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

  // Record the exact `nextflow run` invocation, for the report and for the record
  ch_command_txt = channel.of(workflow.commandLine)
    .collectFile(
        storeDir: "${params.tracedir}",
        name:     "execution_command.txt",
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
  ch_barcodes = channel.value(params.barcodes)

  // Input file with multiplexed reads (FASTQ.gz)
  ch_input = channel.value(params.input)

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
