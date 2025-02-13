#!/usr/bin/env nextflow
/*

============================================================================
  NextITS: Pipeline to process fungal ITS amplicons
============================================================================
  License: Apache-2.0
  Github : https://github.com/vmikk/NextITS
  Website: https://next-its.github.io/
----------------------------------------------------------------------------
*/

// NB!!:
// - provide absolute paths to the input data (e.g. --input and --barcodes)
// - File names should not contain period (.) characters (except for extensions)

// Databases:
//  - UDB for chimera identification


// Enable DSL2 syntax
nextflow.enable.dsl = 2

// Enable topic channels
nextflow.preview.topic = true

// Include the software version parser function
include { software_versions_to_yaml } from './modules/version_parser.nf'


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  VALIDATE INPUTS

// Validate & print parameter summary
// NB! works only with old schema (`everit-json-schema` library doesn't support JSON Schema draft-2020-12)
WorkflowMain.initialise(workflow, params, log)

// Include the pipeline initialisation subworkflow
// requires newer nf-core template and schema
// include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_NextITS_pipeline'



// Pipeline help message
def helpMsg() {
    log.info"""
    =====================================================================
    NextITS v.${workflow.manifest.version}
    =====================================================================
    
    Pipeline Usage:
    To run the pipeline, enter the following in the command line:
        nextflow run vmikk/nextits -r ${workflow.manifest.version} --input ... --outdir ...
    
    Options:
    REQUIRED:
        --input           File with single-end input sequences, PacBio (FASTQ or BAM) or a directory with pre-demultiplexed files
        --input_R1        Files with paired-end input sequences, Illumina (FASTQ)
        --input_R2
        --barcodes        Barcodes for demultiplexing (FASTA)
        --outdir          The output directory where the results will be saved

    OPTIONAL:
        --demultiplexed     Boolean, input is multiplexed (true, single FASTQ file) or pre-demultiplexed (multiple FASTQ files)
        --seqplatform       Sequencing platform type - "PacBio" (default) or "Illumina"
        --its_region        ITS part selector - "full" (defalut), "ITS1", "ITS2", "none" (trims primers only), or "ITS1_5.8S_ITS2"
        --primer_forward    Forward primer sequence (default, ITS9mun)
        --primer_reverse    Reverse primer sequence (default, ITS4ngsUni)
        --primer_mismatches
        --primer_foverlap   Min primer overlap (default, F primer length - 2)
        --primer_roverlap   Min primer overlap (default, R primer length - 2)
        --qc_maxn           Discard sequences with more than the specified number of N’s 
        --trim_minlen       Min sequence length after primer trimming (default, 10)
        --ITSx_tax          ITSx taxonomy profile (default, "all")
        --ITSx_evalue       ITSx E-value cutoff threshold (default, 1e-1)
        --ITSx_partial      Keep partial ITS sequences (defalt, off), otherwise specify min length cutoff
        --hp                Homopolymer compression (default, true)
        --hp_similarity     Allowed sequence similarity for homopolymer compression (default, 0.999)
        --hp_iddef          Sequence similarity definition for homopolymer compression (default, 2)
        
      # Chimera identification
        --chimera_db               Database for reference-based chimera removal
        --chimera_rescueoccurrence Min occurrence of chimeric sequences required to rescue them (default, 2)
        --chimeranov_abskew   De novo chimera identification `abskew` parameter (default, 2.0)
        --chimeranov_dn       De novo chimera identification `dn` parameter (default, 1.4)
        --chimeranov_mindiffs De novo chimera identification `mindiffs` parameter (default, 3)
        --chimeranov_mindiv   De novo chimera identification `mindiv` parameter (default, 0.8)
        --chimeranov_minh     De novo chimera identification `minh` parameter (default, 0.28)
        --chimeranov_xn       De novo chimera identification `xn` parameter (default, 8.0)

      # Tag-jump removal
        --tj_f      Tag-jump filtering, UNCROSS parameter `f` (default, 0.01)
        --tj_p      Tag-jump filtering parameter `p` (default, 1)
        --otu_id    Sequence similarity for OTU clustering (default, 0.98)
        --otu_iddef Sequence similarity definition for tag-jump removal step (default, 2)

      # PacBio-specific parameters
        --lima_barcodetype Tag type ("single", "dual", "dual_symmetric", "dual_asymmetric")
        --lima_minscore    Minimum barcode score for demultiplexing (default, 93)
        --lima_minendscore Minimum second barcode score (only for asymmetric and dual barcoding scheme; default, 50)
        --lima_minrefspan  Minimum read span relative to the barcode length (0-1; default, 0.75)
        --lima_minscoringregions Number of barcodes scored required for demultiplexing using dual barcodes (default, 2 = requires both barcodes)
        --lima_windowsize  Window size for barcode lookup (default, 70 bp)
        --lima_minlen      Minimum sequence length after clipping barcodes (default, 40)
        --qc_maxee         Maximum number of expected errors (default, false)
        --qc_maxeerate     Maximum number of expected errors per base (default, 0.01)
        --qc_maxhomopolymerlen  Threshold for a homopolymer region length in a sequence (default, 25)

      # Illumina-specific parameters
        --qc_avgphred      Average Phred score for QC (default, false)
        --qc_twocolor      Enable two-color chemistry mode, e.g. for Illumina NovaSeq (default, false)
        --qc_phredmin      Two-color mode: min Phred score of qualified bases (default, 24)
        --qc_phredperc     Two-color mode: Percentage of bases allowed to be unqualified (default, 30)
        --qc_polyglen      Two-color mode: minimum length of polyG tail (default, 8)
        --barcode_window   Window size for barcode lookup (default, 30 bp)
        --barcode_errors   Maximum allowed number of errors in barcodes (default, 1)
        --barcode_overlap  Min overlap between read and barcode (default, 11)
        --pe_minoverlap    Min length to detect overlapped region of PE reads (default, 20)
        --pe_difflimit     Max number of mismatched bases in PE overlap (default, 5)
        --pe_diffperclimit Max percentage of mismatched bases in PE overlap (default, 20)
        --pe_minlen        Min length of merged sequences (default, 30)
        --illumina_keep_notmerged  Keep not merged Illumina reads (default, true)
        --illumina_joinpadgap      Join not merged reads into one sequence using padding sequence string (default, NNNNNNNNNN)
        --illumina_joinpadqual     Join not merged reads into one sequence using padding quality string (default, IIIIIIIIII)

      # Miscellaneous parameters
        --gzip_compression Compression level for GZIP (default, 7; 1 = fastest (worst compression), 9 = slowest (best))

    NEXTFLOW-SPECIFIC:
        -profile       Configuration profile
        -resume        Execute the pipeline using the cached results (e.g., in case of )
        -work-dir      Path to the directory where intermediate result files are stored
        -qs            Queue size (max number of processes that can be executed in parallel); e.g., 8
        -r             Pipeline version to run (GitHub branch, tag, or SHA number)
    """.stripIndent()
}
// Show help msg
if (params.helpMsg){
    helpMsg()
    exit(0)
}

// Additional parameter validation
if (params.input == false && params.seqplatform == "PacBio") {
    println( "ERROR: Please provide the input file with sequences in FASTQ.gz or BAM format with `--input` parameter.")
    exit(1)
}
if (params.input_R1 == false && params.input_R2 == false && params.seqplatform == "Illumina") {
    println( "ERROR: Please provide input files with sequences in FASTQ.gz format with `--input_R1` and `--input_R2` parameters.")
    exit(1)
}
if (params.barcodes == false && params.demultiplexed == false) {
    println( "ERROR: Please provide the file with sample barcodes in FASTA format with `--barcodes` parameter.")
    exit(1)
}
if (!params.chimera_db || !file(params.chimera_db).exists()) {
    println( "ERROR: Please provide the UDB file with reference sequences for chimera removal with `--chimera_db` parameter.")
    exit(1)
}
if (!(params.chimera_db.toLowerCase().endsWith('.udb'))) {
    println( "ERROR: The reference database file specified with `--chimera_db` parameter must be in UDB format." )
    exit 1
}
if (params.hp == true && params.seqplatform == "Illumina" && params.illumina_keep_notmerged == true) {
    println( "ERROR: Homopolymer compression is not implemented for Illumina non-merged reads.")
    exit(1)
}
if (params.seqplatform == "Illumina" && params.demultiplexed == true) {
    println( "ERROR: Handling demultiplexed data for Illumina is not implemented yet.")
    exit(1)
}


if (params.seqplatform == "Illumina" && params.illumina_keep_notmerged == true && params.its_region != "none") {
    println( "WARNING: Unmerged Illumina reads are not compatible with ITSx. Amplicons will be primer-trimmed.")
}


// Print the parameters to the console and to the log
/*
log.info """
    =======================================================================
    NextITS v.${workflow.manifest.version}
    =======================================================================
    Input data path: ${params.input}
    Barcodes:        ${params.barcodes}
    Output path:     ${params.outdir}
    ITS region:      ${params.its_region}
    """
    .stripIndent()

log.info """
        Pipeline info:
          Pipeline profile:       ${workflow.profile}
          Config file used:       ${workflow.configFiles}
          Container engine:       ${workflow.containerEngine}
        """
        .stripIndent()

log.info """
        Core Nextflow options:
          launchDir:              ${workflow.launchDir}
          workDir:                ${workflow.workDir}
          projectDir:             ${workflow.projectDir}
        """
        .stripIndent()

log.info "======================================================================="
log.info "\n"
*/


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
out_5_chim   = params.outdir + "/05_Chimera"
out_6_tj     = params.outdir + "/06_TagJumpFiltration"
out_7_seq    = params.outdir + "/07_SeqTable"
out_8_smr    = params.outdir + "/08_RunSummary"
out_9_db     = params.outdir + "/09_DB"

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
    echo -e "Converting BAM to FASTQ\n"
    echo -e "Input file: " ${input}
    echo -e "BAM index: "  ${bam_index}

    bam2fastq \
      -c ${params.gzip_compression} \
      --num-threads ${task.cpus} \
      ${input}

    echo -e "\nConvertion finished"
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
    echo -e "QC"
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

    echo -e "\nQC finished"
    """
}


// Quality filtering for pair-end reads
process qc_pe {

    label "main_container"

    // cpus 10

    input:
      path input_R1
      path input_R2

    output:
      path "QC_R1.fq.gz", emit: filtered_R1
      path "QC_R2.fq.gz", emit: filtered_R2

    script:
    filter_avgphred  = params.qc_avgphred  ? "--average_qual ${params.qc_avgphred}"                 : "--average_qual 0"
    filter_phredmin  = params.qc_phredmin  ? "--qualified_quality_phred ${params.qc_phredmin}"      : ""
    filter_phredperc = params.qc_phredperc ? "--unqualified_percent_limit ${params.qc_phredperc}"   : ""
    filter_polyglen  = params.qc_polyglen  ? "--trim_poly_g --poly_g_min_len ${params.qc_polyglen}" : ""
    """
    echo -e "QC"
    echo -e "Input R1: " ${input_R1}
    echo -e "Input R2: " ${input_R2}

    ## If `filter_phredmin` && `filter_phredperc` are specified,
    # Filtering based on percentage of unqualified bases
    # how many percents of bases are allowed to be unqualified (Q < 24)
        
    fastp \
      --in1 ${input_R1} \
      --in2 ${input_R2} \
      --disable_adapter_trimming \
      --n_base_limit ${params.qc_maxn} \
      ${filter_avgphred} \
      ${filter_phredmin} \
      ${filter_phredperc} \
      ${filter_polyglen} \
      --length_required 100 \
      --thread ${task.cpus} \
      --html qc.html \
      --json qc.json \
      --out1 QC_R1.fq.gz \
      --out2 QC_R2.fq.gz

    echo -e "\nQC finished"
    """
}



// Validate tags for demultiplexing
process tag_validation {

    label "main_container"
    // cpus 1

    input:
      path barcodes

    output:
      path "barcodes_validated.fasta", emit: fasta
      path "biosamples_asym.csv",      emit: biosamples_asym, optional: true
      path "biosamples_sym.csv",       emit: biosamples_sym,  optional: true
      path "file_renaming.tsv",        emit: file_renaming,   optional: true

    script:
    """
    echo -e "Valdidating demultiplexing tags\n"
    echo -e "Input file: " ${barcodes}

    ## Convert Windows-style line endings (CRLF) to Unix-style (LF)
    LC_ALL=C sed -i 's/\r\$//g' ${barcodes}

    ## Perform tag validation
    validate_tags.R \
      --tags   ${barcodes} \
      --output barcodes_validated.fasta

    echo -e "Tag validation finished"
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
      path biosamples_sym    // for dual or asymmetric barcodes
      path biosamples_asym   // for dual or asymmetric barcodes
      path file_renaming     // for dual or asymmetric barcodes

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
        echo -e "\nERROR: Symmetric tags are provided in '...' format.\n"
        echo -e "In the FASTA file, please include only one tag per sample, since these tags are identical.\n"
        exit 1
    fi

    ## Count the number of samples in Biosample files - only for `dual` and `dual_asymmetric` barcodes
    if [[ ${params.lima_barcodetype} == "dual_asymmetric" ]] || [[ ${params.lima_barcodetype} == "dual" ]]; then

      if [ ! -e ${biosamples_asym} ]; then
        
        echo -e "\nERROR: Tags are specified in wrong format"
        echo -e "Use the '...' format in FASTA file.\n"
        exit 1
      
      else
        line_count_sym=\$(wc  -l < ${biosamples_sym})
        line_count_asym=\$(wc -l < ${biosamples_asym})

        echo -e "..Number of lines in symmetric file: "  \$line_count_sym
        echo -e "..Number of lines in asymmetric file: " \$line_count_asym

        ## Check the presence of dual barcode combinations
        ## If line count is less than 2, it means there are no samples specified
        if [[ ${params.lima_barcodetype} == "dual_asymmetric" ]] && [[ \$line_count_asym -lt 2 ]]; then
          echo -e "\nERROR: No asymmetric barcodes detected for demultiplexing.\n"
          return 1
        fi

        if [[ ${params.lima_barcodetype} == "dual" ]] && [[ \$line_count_asym -lt 2 ]]; then
          echo -e "\nWARNING: No asymmetric barcodes detected, consider using '--lima_barcodetype dual_symmetric'.\n"
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
        echo -e "\nDemultiplexing with LIMA (single barcode)"
        lima --same --single-side \
          --log-file LIMA/_log.txt \
          \$common_args \
          "LIMA/lima.fq.gz"
        ;;

      "dual_symmetric")
        echo -e "\nDemultiplexing with LIMA (dual symmetric barcodes)"
        lima --same \
          --min-end-score       ${params.lima_minendscore} \
          --min-scoring-regions ${params.lima_minscoringregions} \
          --log-file LIMA/_log.txt \
          \$common_args \
          "LIMA/lima.fq.gz"
        ;;

      "dual_asymmetric")
        echo -e "\nDemultiplexing with LIMA (dual asymmetric barcodes)"
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
        echo -e "\nDemultiplexing with LIMA (dual symmetric barcodes)"
        lima --same \
          --min-end-score       ${params.lima_minendscore} \
          --min-scoring-regions ${params.lima_minscoringregions} \
          --biosample-csv       ${biosamples_sym} \
          --log-file LIMAs/_log.txt \
          \$common_args \
          "LIMAs/lima.fq.gz"
        fi

        if [[ \$line_count_asym -ge 2 ]]; then
        echo -e "\nDemultiplexing with LIMA (dual asymmetric barcodes)"
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

      echo -e "\nPooling of symmetric and asymmetric barcodes"
      cd LIMA
      find ../LIMAd -name "*.fq.gz" | parallel -j1 "ln -s {} ."
      find ../LIMAs -name "*.fq.gz" | parallel -j1 "ln -s {} ."
      cd ..

    fi


    ## Rename barcode combinations into sample names
    ## Only user-provided combinations whould be kept (based on `lima --biosample-csv`)
    if [[ ${params.lima_barcodetype} == "dual_asymmetric" ]] || [[ ${params.lima_barcodetype} == "dual" ]]; then

      echo -e "\n..Renaming files from tag IDs to sample names"
      brename -p "(.+)" -r "{kv}" -k ${file_renaming} LIMA/

    fi

    if [[ ${params.lima_barcodetype} == "dual_symmetric" ]] || [[ ${params.lima_barcodetype} == "single" ]]; then

      echo -e "\n..Renaming demultiplexed files"
      rename --filename \
        's/^lima.//g; s/--.*\$/.fq.gz/' \
        \$(find LIMA -name "*.fq.gz")
    
    fi


    ## Combine summary stats for dual barcodes (two LIMA runs)
    if [[ ${params.lima_barcodetype} == "dual" ]]; then

      echo -e "\n..Combining dual-barcode log files"

      if [ -f "LIMAd/lima.lima.summary" ]; then
        echo -e "Asymmetric barcodes summary\n\n" >> LIMA/lima.lima.summary
        cat LIMAd/lima.lima.summary >> LIMA/lima.lima.summary

        echo -e "Asymmetric barcodes counts\n\n" >> LIMA/lima.lima.counts
        cat LIMAd/lima.lima.counts >> LIMA/lima.lima.counts

        echo -e "Asymmetric barcodes report\n\n" >> LIMA/lima.lima.report
        cat LIMAd/lima.lima.report >> LIMA/lima.lima.report
      fi

      if [ -f "LIMAs/lima.lima.summary" ]; then
        echo -e "\n\nSymmetric barcodes summary\n\n" >> LIMA/lima.lima.summary
        cat LIMAs/lima.lima.summary >> LIMA/lima.lima.summary

        echo -e "\n\nSymmetric barcodes counts\n\n" >> LIMA/lima.lima.counts
        cat LIMAs/lima.lima.counts >> LIMA/lima.lima.counts

        echo -e "\n\nSymmetric barcodes report\n\n" >> LIMA/lima.lima.report
        cat LIMAs/lima.lima.report >> LIMA/lima.lima.report
      fi

    fi  # end of dual logs pooling


    ## Compress logs
    echo -e "..Compressing log file"
    gzip -${params.gzip_compression} LIMA/lima.lima.report


    ## LIMA defaults:
    # SYMMETRIC  : --ccs --min-score 0 --min-end-score 80 --min-ref-span 0.75 --same --single-end
    # ASYMMETRIC : --ccs --min-score 80 --min-end-score 50 --min-ref-span 0.75 --different --min-scoring-regions 2

    echo -e "\nDemultiplexing finished"
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
    echo -e "Merging Illumina pair-end reads"

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

    echo -e "\n..done"

    ## Remove empty files (no sequences)
    echo -e "\nRemoving empty files"
    find . -type f -name "*.fq.gz" -size -29c -print -delete
    echo -e "..Done"

    echo -e "\nDemultiplexing finished"
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
    echo -e "\nDisambiguating reverse primer"
    disambiguate_primers.R \
      ${params.primer_reverse} \
      primer_R.fasta

    ## Reverse-complement primers
    echo -e "\nReverse-complementing primers"
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

    echo -e "\nCounting primers"
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
    echo -e "\nLooking for multiple primer occurrences"
    
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

      echo -e "\nNumber of artefacts found: " \$(wc -l < multiprimers.txt)

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

      echo -e "\nNo primer artefacts found"
      ln -s ${input} no_multiprimers.fq.gz
    
    fi
    echo -e "..Done"

    echo -e "\nReorienting sequences"

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

    echo -e "\nAll done"

    ## Clean up
    if [ -f no_multiprimers.fq.gz ]; then rm no_multiprimers.fq.gz; fi

    ## Remove empty file (no valid sequences)
    echo -e "\nRemoving empty files"
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
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.extraction.results", emit: itsx_details, optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.SSU.full_and_partial.fasta.gz",  emit: itsx_ssu_part,  optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.ITS1.full_and_partial.fasta.gz", emit: itsx_its1_part, optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.5_8S.full_and_partial.fasta.gz", emit: itsx_58s_part,  optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.ITS2.full_and_partial.fasta.gz", emit: itsx_its2_part, optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.LSU.full_and_partial.fasta.gz",  emit: itsx_lsu_part,  optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}_primertrimmed_sorted.fq.gz",     emit: trimmed_seqs,   optional: true
      path "parquet/*.parquet", emit: parquet, optional: true

    script:
    
    sampID="${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}"

    // Allow inclusion of sequences that only find a single domain, given that they meet the given E-value and score thresholds, on with parameters 1e-9,0 by default
    // singledomain = params.ITSx_singledomain ? "--allow_single_domain 1e-9,0" : ""

    """
    echo "Extraction of rRNA regions using ITSx\n"

    ## Trim primers    
    echo -e "Trimming primers\n"

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

    echo -e "..Done\n"

    ## Check if there are sequences in the output
    NUMSEQS=\$( seqkit stat --tabular --quiet ${sampID}_primertrimmed.fq.gz | awk -F'\t' 'NR==2 {print \$4}' )
    echo -e "Number of sequences after primer trimming: " \$NUMSEQS
    if [ \$NUMSEQS -lt 1 ]; then
      echo -e "\nIt looks like no reads remained after trimming the primers\n"
      exit 0
    fi
   
    ## Estimate sequence quality and sort sequences by quality
    echo -e "\nSorting by sequence quality"
    seqkit replace -p "\\s.+" ${sampID}_primertrimmed.fq.gz \
      | phredsort -i - -o - --metric meep --header avgphred,maxee,meep \
      | gzip -1 > ${sampID}_primertrimmed_sorted.fq.gz
    echo -e "..Done"

    ## Hash sequences, add sample ID to the header
    ## columns: Sample ID - Hash - PacBioID - AvgPhredScore - MaxEE - MEEP - Sequence - Quality - Length
    ## Convert to Parquet format
    echo -e "\nCreating hash table"
    seqhasher --hash sha1 --name ${sampID} ${sampID}_primertrimmed_sorted.fq.gz - \
      | seqkit fx2tab --length \
      | sed 's/;/\t/ ; s/;/\t/ ; s/ avgphred=/\t/ ; s/ maxee=/\t/ ; s/ meep=/\t/' \
      > ${sampID}_hash_table.txt
    echo -e "..Done"

    ## Check the number of fields per record (should be 9!)
    # awk '{print NF}' ${sampID}_hash_table.txt | sort | uniq -c
    # awk 'NF > 9 {print \$0 }' ${sampID}_hash_table.txt

    ## Dereplicate at sample level (use quality-sorted sequences to make sure that the representative sequence is with the highest quality)
    echo -e "\nDereplicating at sample level"
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
    echo -e "\nITSx extraction"
    ITSx \
      -i derep.fasta \
      --complement T \
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
      echo -e "Partial files found, removing empty sequences\n."

      find . -name "*.full_and_partial.fasta" \
        | parallel -j${task.cpus} "seqkit seq -m 1 -w 0 {} > {.}_tmp.fasta"

      rm *.full_and_partial.fasta
      brename -p "_tmp" -r "" -f "_tmp.fasta\$"

    fi


    ## Remove empty files (no sequences)
    echo -e "\nRemoving empty files"
    find . -type f -name "*.fasta" -empty -print -delete
    echo -e "..Done"

    ## Remove temporary file
    rm derep.fasta
    rm ${sampID}_primertrimmed.fq.gz

    ## Compress results
    echo -e "\nCompressing files"

    parallel -j${task.cpus} "gzip -${params.gzip_compression} {}" ::: \
      ${sampID}_hash_table.txt \
      ${sampID}_uc.uc \
      *.fasta

    ## Convert ITSx output to Parquet
    if [ ${params.ITSx_to_parquet} == true ]; then

      echo -e "\nConverting ITSx output to Parquet"
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

      echo -e "Parquet files created\n"

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

    echo -e "\n..Done"
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

    script:
    sampID="${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}"

    """
    echo -e "Trimming primers\n"
    echo -e "Input sample: "   ${sampID}
    echo -e "Forward primer: " ${params.primer_forward}
    echo -e "Reverse primer: " ${params.primer_reverse}

    ## Reverse-complement rev priver
    RR=\$(rc.sh ${params.primer_reverse})
    echo -e "Reverse primer RC: " "\$RR"

    echo -e "\nTrimming primers"
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
      echo -e "\nSorting by sequence quality"
      seqkit replace -p "\\s.+" ${sampID}_primertrimmed.fq.gz \
        | phredsort -i - -o - --metric meep --header avgphred,maxee,meep \
        | gzip -${params.gzip_compression} \
        > ${sampID}_primertrimmed_sorted.fq.gz
      echo -e "..Done"

      rm ${sampID}_primertrimmed.fq.gz

      ## Hash sequences, add sample ID to the header
      ## columns: Sample ID - Hash - PacBioID - AvgPhredScore - MaxEE - MEEP - Sequence - Quality - Length
      ## Convert to Parquet format
      echo -e "\nCreating hash table"
      seqhasher --hash sha1 --name ${sampID} ${sampID}_primertrimmed_sorted.fq.gz - \
        | seqkit fx2tab --length \
        | sed 's/;/\t/ ; s/;/\t/ ; s/ avgphred=/\t/ ; s/ maxee=/\t/ ; s/ meep=/\t/' \
        > ${sampID}_hash_table.txt
      echo -e "..Done"

      ## Compress results
      echo -e "Compressing result"
      gzip -${params.gzip_compression} ${sampID}_hash_table.txt

      ## Dereplicate at sample level
      echo -e "\nDereplicating at sample level"
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

      echo -e "\nNo sequences found after primer removal"
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

      echo -e "\n..All parts found"

      ## Prepare tables for ID matching
      echo -e "\n..Converting data to tabular format"
      seqkit fx2tab ${ITS1} | sed 's/\t\$//g' | csvtk add-header -t -n id,ITS1 > tmp_1_ITS1.txt
      seqkit fx2tab ${S58}  | sed 's/\t\$//g' | csvtk add-header -t -n id,58S  > tmp_1_s58.txt
      seqkit fx2tab ${ITS2} | sed 's/\t\$//g' | csvtk add-header -t -n id,ITS2 > tmp_1_ITS2.txt

      ## Join ITS fragments
      echo -e "\n..Joining ITS fragments"
      csvtk join -t -f   "id" tmp_1_ITS1.txt tmp_1_s58.txt tmp_1_ITS2.txt > tmp_2_ITS1_58S_ITS2.txt

      ## Check joining results
      NUMSEQS=\$(wc -l < tmp_2_ITS1_58S_ITS2.txt)
      echo "...Number of joined sequences: " \$((NUMSEQS - 1))

      if [ "\$NUMSEQS" -gt 1 ]; then

        ## Convert table back to fasta
        ## Remove leading and trailing Ns
        echo -e "\n..Preparing fasta"
        awk 'NR>1 { print \$1 "\t" \$2\$3\$4 }' tmp_2_ITS1_58S_ITS2.txt \
          | seqkit tab2fx -w 0  \
          | seqkit replace -p "^n+|n+\$" -r "" -is -w 0 \
          | gzip -${params.gzip_compression} > ${sampID}_ITS1_58S_ITS2.fasta.gz

      else
        echo "...There are no sequences with all ITS parts present\n"
        echo -e "\n..Skipping ITS assembly for this sample"
      fi

    else 
      echo -e "\n..Some or all parts are missing"
      echo -e "\n..Skipping ITS assembly for this sample"
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

    script:
    """
    echo -e "Aggregating sequence qualities"

    merge_hash_tables.sh \
      -i ./hash_tables \
      -o SeqQualities.parquet \
      -t ${task.cpus}

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
    echo -e "\nRe-clustering homopolymer-compressed data"
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
      echo -e "\nExtracting sequences"

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
      echo -e "Clustering homopolymer-compressed sequences returned to results\n"
      echo -e "(most likely, sequences were too short)\n"
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

    script:
    sampID="${input.getSimpleName()}"

    """

     vsearch \
          --derep_fulllength ${input} \
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

    """
}



// Reference-based chimera removal
process chimera_ref {

    label "main_container"

    publishDir "${out_5_chim}", mode: "${params.storagemode}"
    // cpus 1

    // Add sample ID to the log file
    tag "${input.getSimpleName().replaceAll(/_Homopolymer_compressed/, '')}"

    input:
      path input
      path DB

    output:
      path "${input.getSimpleName().replaceAll(/_Homopolymer_compressed/, '')}_NoChimera.fa.gz", emit: nonchimeric, optional: true
      path "${input.getSimpleName().replaceAll(/_Homopolymer_compressed/, '')}_Chimera.fa.gz", emit: chimeric, optional: true

    script:
    sampID="${input.getSimpleName().replaceAll(/_Homopolymer_compressed/, '')}"

    """

    ## Sample name will be added to the header of chimeric sequences
    # sampID="\$(basename ${input} _Homopolymer_compressed.fa.gz)"

    ## Reference database based chimera filtering
    echo -e "Reference-based chimera removal"
    vsearch \
      --uchime_ref ${input} \
      --db ${DB} \
      --selfid \
      --fasta_width 0 \
      --threads ${task.cpus} \
      --sizein --sizeout \
      --chimeras chimeras.fasta \
      --nonchimeras nonchimeras.fasta \
      --borderline borderline.fasta

    # --selfid = ignore reference sequences that are 100% identical to the query
    echo -e "..Done"


    ## Add borderline sequences to non-chimeric sequences
    if [ -e borderline.fasta ]
    then
        echo -e "\nBorderline sequences were added to non-chimeric sequences"
        cat borderline.fasta nonchimeras.fasta > nc_bo.fasta
        mv nc_bo.fasta nonchimeras.fasta
        rm borderline.fasta
    fi

    ## Chimeric sequences
    if [ -e chimeras.fasta ]
    then
      ## Add sample ID to the header and compress the file
      sed 's/>.*/&;sample='"${sampID}"';/' chimeras.fasta \
        | gzip -${params.gzip_compression} \
        > "${sampID}_Chimera.fa.gz"
      rm chimeras.fasta
    else
      echo -e "\nNo chimeras detected"
    fi

    ## Non-chimeric sequences
    if [ -e nonchimeras.fasta ]
    then
      gzip -c nonchimeras.fasta > "${sampID}_NoChimera.fa.gz"
      rm nonchimeras.fasta
    else
      echo "No non-chimeric sequences left"
    fi

    """
}


// Recovery of ref-based chimeric sequences with high occurrence
process chimera_rescue {

    label "main_container"

    publishDir "${out_5_chim}", mode: "${params.storagemode}"
    // cpus 1

    input:
      path input

    output:
      path "*_RescuedChimera.fa.gz", emit: rescuedchimeric, optional: true

    script:
    """

    ## Aggregate chimeric sequences from different samples
    echo -e "\nAggregating chimeric sequences"
    find . -name "*_Chimera.fa.gz" \
      | parallel -j1 "zcat {}" \
      | seqkit fx2tab \
      | sed -r 's:\t+:\t:g' | sed 's/\t\$//g' \
      | gzip -${params.gzip_compression} > All_chimeras.txt.gz
    echo -e "..Done"

    ### Inspect chimerae occurrence
    ## Rescue sequences that were annotated as chimeric, 
    ## but have high occurrence within sequenceing run (e.g., occurrence > 2)
    echo -e "\nInspecting occurrence of chimeric sequences"

    chimera_rescue.R \
      "All_chimeras.txt.gz" \
      ${params.chimera_rescueoccurrence} \
      "Rescued_Chimeric_sequences.fa.gz"

    echo -e "..Done"

    ## Split rescured sequences by sample
    if [ -e Rescued_Chimeric_sequences.fa.gz ]
    then
      echo -e "\n..Splitting rescued sequences by sample"
      seqkit split -i \
        --id-regexp ";sample=(.*);" \
        --threads ${task.cpus} \
        -w 0 \
        -O Rescued_by_sample \
        Rescued_Chimeric_sequences.fa.gz

      rename \
        --filename 's/^Rescued_Chimeric_sequences.id_//g; s/.fa.gz/_RescuedChimera.fa.gz/' \
        \$(find Rescued_by_sample -name "*.fa.gz")

      mv Rescued_by_sample/*_RescuedChimera.fa.gz .

      echo -e "..Done"
    fi

    ## Remove temporary files
    rm All_chimeras.txt.gz
    if [ -f Rescued_Chimeric_sequences.fa.gz ]; then rm Rescued_Chimeric_sequences.fa.gz; fi

    """
}



// De novo chimera identification
// NB! in uchime_denovo, sequences are compared on their plus strand only!
process chimera_denovo {

    label "main_container"

    publishDir "${out_5_chim}", mode: "${params.storagemode}"
    // cpus 1

    // Add sample ID to the log file
    tag "${input.getSimpleName().replaceAll(/_Homopolymer_compressed/, '')}"

    input:
      path input

    output:
      path "${input.getSimpleName().replaceAll(/_Homopolymer_compressed/, '')}_DeNovoChim.txt", emit: denovochim, optional: true

    script:
    sampID="${input.getSimpleName().replaceAll(/_Homopolymer_compressed/, '')}"

    """

    ## Input order matters for chimera detection,
    ## so sequences will be automatically sorted by decreasing abundance first

    echo -e "De novo chimera identification"

    vsearch \
      --uchime_denovo ${input} \
      --abskew ${params.chimeranov_abskew} \
      --dn ${params.chimeranov_dn} \
      --mindiffs ${params.chimeranov_mindiffs} \
      --mindiv ${params.chimeranov_mindiv} \
      --minh ${params.chimeranov_minh} \
      --xn ${params.chimeranov_xn} \
      --threads 1 \
      --qmask dust \
      --sizein --xsize \
      --fasta_width 0 \
      --fasta_score \
      --chimeras - \
    | seqkit seq --name \
    | sed 's/;+/;/g ; s/;/\t/g ; s/uchime_denovo=//' \
    | sed 's/\$/\t${sampID}/' \
    > ${sampID}_DeNovoChim.txt
      
      # --uchimeout uchimeout.txt

    ## Remove file, if empty
    find . -maxdepth 1 -name ${sampID}_DeNovoChim.txt -size 0 -print -delete

    echo -e "..Done"

    """
}

// Aggregate de novo chimeras into a single file
process chimera_denovo_agg {

    label "main_container"
    // cpus 1

    input:
      path input

    output:
      path "DeNovo_Chimera.txt", emit: alldenovochim, optional: true

    script:
    """
    echo -e "Aggregating de novo chimeric sequences"

    find . -name "*_DeNovoChim.txt" \
      | parallel -j1 "cat {}" \
      > DeNovo_Chimera.txt

    echo -e "..Done"

    """
}



// Merge and dereplicate sequences from all samples
process glob_derep {

    label "main_container"

    // cpus 1

    input:
      path input

    output:
      path "Derep_for_clust.fa.gz", emit: globderep
      path "Derep_for_clust.uc.gz", emit: globderep_uc

    script:
    """

    ## Concatenate non-chimeric sequences
    # zcat *_NoChimera.fa.gz *_RescuedChimera.fa.gz

    echo -e "\n ===== Files with non-chimeric sequences:"
    if test -n "\$(find . -maxdepth 1 -name '*_NoChimera.fa.gz' -print -quit)"
    then
        ls -l *_NoChimera.fa.gz;
    else
        echo "..No files"
    fi

    echo -e "\n ===== Files with rescued chimeric sequences:"
    if test -n "\$(find . -maxdepth 1 -name '*_RescuedChimera.fa.gz' -print -quit)"
    then
        ls -l *_RescuedChimera.fa.gz;
    else
        echo "..No files"
    fi

    echo -e "\nDereplicating sequences"
    find . -name "*Chimera.fa.gz" | parallel -j1 \
      "zcat {}" \
      | sed '/^>/ s/;sample=.*;/;/' \
      | vsearch \
        --derep_fulllength - \
        --output - \
        --strand both \
        --fasta_width 0 \
        --threads 1 \
        --sizein --sizeout \
        --uc Derep_for_clust.uc \
      | gzip -${params.gzip_compression} > Derep_for_clust.fa.gz
    echo -e "..Done"

    ## Compress UC file
    echo -e "\nCompressing UC file"
    gzip -${params.gzip_compression} Derep_for_clust.uc

    """
}


// De-novo clustering of sequences into OTUs (for tag-jump removal)
process otu_clust {

    label "main_container"

    publishDir "${out_6_tj}", mode: "${params.storagemode}"
    // cpus 10

    input:
      path input

    output:
      path "OTUs.fa.gz", emit: otus
      path "OTUs.uc.gz", emit: otus_uc

    script:
    """

    echo -e "\nClustering sequences"
    vsearch \
      --cluster_size ${input} \
      --id ${params.otu_id} \
      --iddef ${params.otu_iddef} \
      --sizein --sizeout \
      --qmask dust --strand both \
      --fasta_width 0 \
      --uc OTUs.uc \
      --threads ${task.cpus} \
      --centroids - \
    | gzip -${params.gzip_compression} > OTUs.fa.gz
    
    echo -e "..Done"

    ## Compress UC file
    echo -e "\nCompressing UC file"
    pigz -p ${task.cpus} -${params.gzip_compression} OTUs.uc

    """
}


// Pool sequences and add sample ID into header (for OTU and "ASV" table creation)
process pool_seqs {

    label "main_container"
    
    publishDir "${out_6_tj}", mode: "${params.storagemode}"
    // cpus 2

    input:
      path input

    output:
      path "Seq_tab_not_filtered.txt.gz", emit: seqtabnf
      path "Seq_not_filtered.fa.gz", emit: seqsnf

    script:
    """

    echo -e "\nPooling and renaming sequences"

    ## If there is a sample ID in the header already, remove it
    parallel -j 1 --group \
      "zcat {} \
        | sed -r '/^>/ s/;sample=[^;]*/;/g ; s/;;/;/g' \
        | sed 's/>.*/&;sample='{/.}';/ ; s/_NoChimera.fa//g ; s/_RescuedChimera.fa//g  ; s/_JoinedPE//g' \
        | sed 's/Rescued_Chimeric_sequences.part_//g' \
        | sed -r '/^>/ s/;;/;/g'" \
      ::: *.fa.gz \
      | gzip -${params.gzip_compression} \
      > Seq_not_filtered.fa.gz

    echo "..Done"

    echo -e "\nExtracting sequence count table"
    seqkit seq --name Seq_not_filtered.fa.gz \
      | sed 's/;/\t/g; s/size=//; s/sample=// ; s/\t*\$//' \
      | gzip -${params.gzip_compression} \
      > Seq_tab_not_filtered.txt.gz

    echo "..Done"

    """
}


// Create OTU table (for tag-jump removal)
process otu_tab {

    label "main_container"

    publishDir "${out_6_tj}", mode: "${params.storagemode}"
    // cpus 10

    input:
      path otus
      path sample_seqs

    output:
      path "OTU_tab_not_filtered.txt.gz", emit: otutab
      path "Sample_mapping.uc.gz", emit: samples_uc

    script:
    """
    
    ## Read mapping (for OTU table)
    echo -e "\nRead mapping"
    vsearch \
      --usearch_global "${sample_seqs}" \
      --db "${otus}" \
      --id ${params.otu_id} \
      --strand both \
      --qmask none \
      --dbmask none \
      --sizein --sizeout \
      --otutabout "OTU_tab_not_filtered.txt" \
      --uc Sample_mapping.uc \
      --threads ${task.cpus}

    echo -e "..Done"

    ## Compress UC file
    echo -e "\nCompressing results"
    pigz -p ${task.cpus} -${params.gzip_compression} OTU_tab_not_filtered.txt
    pigz -p ${task.cpus} -${params.gzip_compression} Sample_mapping.uc

    """
}


// Tag-jump removal
process tj {

    label "main_container"

    publishDir "${out_6_tj}", mode: "${params.storagemode}"
    // cpus 1

    input:
      path otutab

    output:
      path "OTU_tab_TagJumpFiltered.txt.gz", emit: otutabtj
      path "TagJump_OTUs.RData", emit: tjs
      path "TagJump_plot.pdf"

    script:
    """

    echo -e "Tag-jump removal"
    
    tag_jump_removal.R \
      ${otutab} \
      ${params.tj_f} \
      ${params.tj_p}

    echo "..Done"

    """
}


// Prepare a table with non-tag-jumped sequences
// Add quality estimate to singletons
// Add chimera-scores for putative de novo chimeras
process prep_seqtab {

    label "main_container"

    publishDir "${out_7_seq}", mode: "${params.storagemode}"
    // cpus 1

    input:
      path seqtabnf
      path seqsnf
      path mappings
      path tagjumps
      path denovos
      path quals

    output:
      path "Seqs.txt.gz",    emit: seq_tl
      path "Seq_tab.txt.gz", emit: seq_tw
      path "Seqs.fa.gz",     emit: seq_fa
      path "Seqs.RData",     emit: seq_rd
      tuple val("${task.process}"), val('R'), eval('Rscript -e "cat(R.version.string)" | sed "s/R version //"'),  topic: versions
      tuple val("${task.process}"), val('data.table'), eval('Rscript -e "cat(as.character(packageVersion(\'data.table\')))"'),  topic: versions
      tuple val("${task.process}"), val('arrow'), eval('Rscript -e "cat(as.character(packageVersion(\'arrow\')))"'),  topic: versions
      tuple val("${task.process}"), val('Biostrings'), eval('Rscript -e "cat(as.character(packageVersion(\'Biostrings\')))"'),  topic: versions

    script:
    """

    echo -e "Sequence table creation"
    
    seq_table_assembly.R \
      ${seqtabnf} \
      ${seqsnf}   \
      ${mappings} \
      ${tagjumps} \
      ${denovos}  \
      ${quals}

    echo "..Done"

    """
}



// Demultiplexing with cutadapt - for Illumina PE reads (only not merged)
// NB. it's possible to use anchored adapters (e.g., -g ^file:barcodes.fa),
//     but there could be a preceding nucleotides before the barcode,
//     therefore, modified barcodes would be used here (e.g., XN{30})
process demux_illumina_notmerged {

    label "main_container"

    publishDir "${out_1_demux}", mode: "${params.storagemode}"
    // cpus 20

    input:
      path input      // tuple of size 2
      path barcodes   // barcodes_modified.fa (e.g., XN{30})

    output:
      tuple path("Combined/*.R1.fastq.gz"), path("Combined/*.R2.fastq.gz"), emit: demux_pe, optional: true
      path "NonMerged_samples.txt", emit: samples_nonm_pe, optional: true

    script:
    """
    echo -e "\nDemultiplexing not-merged reads"

    echo -e "Input R1: " ${input[0]}
    echo -e "Input R2: " ${input[1]}
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
      ${input[0]} ${input[1]} \
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
    echo -e "\nRemoving unknowns"
    rm round1-unknown.R{1,2}.fastq.gz
    rm round2-unknown.R{1,2}.fastq.gz
    echo -e "..Done"

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
    echo -e "\nRemoving temporary files"
    find . -type f -name "round*.fastq.gz" -print -delete

    echo -e "\nDemultiplexing finished"
    """
}


// Trim primers of nonmerged PE reads
// + Estimate sequence qualities
process trim_primers_pe {

    label "main_container"

    publishDir "${out_3_trimPE}", mode: "${params.storagemode}"
    // cpus 2

    // Add sample ID to the log file
    tag "${input[0].getSimpleName()}"

    input:
      path input   // tuple of size 2

    output:
      tuple path("${input[0].getSimpleName()}_R1.fa.gz"), path("${input[0].getSimpleName()}_R2.fa.gz"), emit: nm_FA, optional: true
      tuple path("${input[0].getSimpleName()}_hash_table_R1.txt.gz"), path("${input[0].getSimpleName()}_hash_table_R2.txt.gz"), emit: nm_hashes, optional: true
      tuple path("${input[0].getSimpleName()}_R1.fq.gz"), path("${input[0].getSimpleName()}_R2.fq.gz"), emit: nm_FQ, optional: true
      tuple path("${input[0].getSimpleName()}_uc_R1.uc.gz"), path("${input[0].getSimpleName()}_uc_R2.uc.gz"), emit: nm_UC, optional: true

    script:
    sampID="${input[0].getSimpleName()}"
    
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
      zcat for_R1.fastq.gz | seqkit replace -p "\\s.+" | gzip -${params.gzip_compression} > OK_R1.fastq.gz
      zcat for_R2.fastq.gz | seqkit replace -p "\\s.+" | gzip -${params.gzip_compression} > OK_R2.fastq.gz
    fi

    if [ -s rev_R1.fastq.gz ]; then
      echo -e "..Adding sequences to the main pool"
      zcat rev_R1.fastq.gz | seqkit replace -p "\\s.+" | gzip -${params.gzip_compression} >> OK_R1.fastq.gz
      zcat rev_R2.fastq.gz | seqkit replace -p "\\s.+" | gzip -${params.gzip_compression} >> OK_R2.fastq.gz

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
      | gzip -${params.gzip_compression} \
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
      | gzip -${params.gzip_compression} \
      > ${sampID}_R2.fa.gz


      echo -e "..Done"

      ## Compress results
      echo -e "Compressing result"
      gzip -${params.gzip_compression} ${sampID}_hash_table_R1.txt
      gzip -${params.gzip_compression} ${sampID}_hash_table_R2.txt
      gzip -${params.gzip_compression} ${sampID}_uc_R1.uc
      gzip -${params.gzip_compression} ${sampID}_uc_R2.uc


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

    echo -e "\nJoining with N-pads"
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

      echo -e "\nJoining without N-pads (for quality estimation)"
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
      echo -e "\nCreating sequence hash table with average sequence quality"
        
      seqkit fx2tab --length --avg-qual tmp_for_qual.fq.gz \
        | hash_sequences.sh \
        | awk '{print \$1 "\t" \$6 "\t" \$4 "\t" \$5}' \
        > tmp_hash_table.txt
      
      echo -e "..Done"

      ## Estimating MaxEE
      echo -e "\nEstimating maximum number of expected errors per sequence"

      vsearch \
          --fastx_filter tmp_for_qual.fq.gz \
          --fastq_qmax 93 \
          --eeout \
          --fastaout - \
        | seqkit seq --name \
        | sed 's/;ee=/\t/g' \
        > tmp_ee.txt

      echo -e "..Done"

      echo -e "\nMerging quality estimates"

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
      echo -e "\nIt looks like there are no joined reads"
    fi

    ## Remove redundant symlinks
    find -L . -name "*.fastq.gz" | grep -v ${input} | parallel -j1 "rm {}"

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

    script:

    """
    echo -e "Summarizing run statistics\n"
    echo -e "Counting the number of reads in:\n"


    ## Count raw reads
    echo -e "\n..Raw data"
    seqkit stat --basename --tabular --threads ${task.cpus} --quiet \
      1_input/* > Counts_1.RawData.txt
    
    ## Count number of reads passed QC
    echo -e "\n..Sequenced passed QC"
    seqkit stat --basename --tabular --threads ${task.cpus} --quiet \
      2_qc/* > Counts_2.QC.txt
    
    ## Count demultiplexed reads
    echo -e "\n..Demultiplexed data"
    seqkit stat --basename --tabular --threads ${task.cpus} --quiet \
      3_demux/* > Counts_3.Demux.txt
    

    ## Count primer-checked reads
    echo -e "\n..Primer-checked data"
    if [ `find 4_primerch -name no_primerchecked 2>/dev/null` ]
    then
      echo -e "... No files found"
      touch Counts_4.PrimerCheck.txt
    else
      seqkit stat --basename --tabular --threads ${task.cpus} --quiet \
        4_primerch/* > Counts_4.PrimerCheck.txt
    fi


    ## Count primer-artefacts
    echo -e "\n..Primer-artefacts"
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
    echo -e "\n..ITSx- or primer-trimmed data"
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
    echo -e "\n..Counting homopolymer-corrected reads"
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
    echo -e "\n..Reference-based chimeras"
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
    echo -e "\n..De novo chimeras"
    if [ `find 7_chimdenov -name no_chimdenovo 2>/dev/null` ]
    then
      echo -e "... No files found"
      touch Counts_7.ChimDenov.txt
    else
      cat 7_chimdenov/* > Counts_7.ChimDenov.txt
    fi


    ## Rescued chimeras
    echo -e "\n..Rescued chimeric sequences"
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

    script:

    """
    echo -e "Summarizing run statistics\n"
    echo -e "Counting the number of reads in:\n"


    ## Count raw reads
    echo -e "\n..Raw data"
    seqkit stat --basename --tabular --threads ${task.cpus} \
      1_input/* > Counts_1.RawData.txt
    
    ## Count number of reads passed QC
    echo -e "\n..Sequenced passed QC"
    seqkit stat --basename --tabular --threads ${task.cpus} \
      2_qc/* > Counts_2.QC.txt
    
    ## Count demultiplexed reads
    echo -e "\n..Demultiplexed data"
    seqkit stat --basename --tabular --threads ${task.cpus} \
      3_demux/* > Counts_3.Demux.txt

    ## Count primer-checked reads
    echo -e "\n..Primer-checked data"
    if [ `find 4_primerch -name no_primerchecked 2>/dev/null` ]
    then
      echo -e "... No files found"
      touch Counts_4.PrimerCheck.txt
    else
      seqkit stat --basename --tabular --threads ${task.cpus} \
        4_primerch/* > Counts_4.PrimerCheck.txt
    fi

    ## Count primer-artefacts
    echo -e "\n..Primer-areifacts"
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

//  The default workflow
workflow {

  // Primer disambiguation
  disambiguate()


  // Demultiplex data
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

      // Demultiplexing with dual barcodes requires 3 additional files
      //  - "biosamples" with symmertic/asymmetirc tag combinations
      //  - and a table for assigning sample names to demuxed files
      // Create dummy files (for single or symmetic tags) if neccesary
      ch_biosamples_sym  = tag_validation.out.biosamples_sym.flatten().collect().ifEmpty(file("biosamples_sym"))
      ch_biosamples_asym = tag_validation.out.biosamples_asym.flatten().collect().ifEmpty(file("biosamples_asym"))
      ch_file_renaming   = tag_validation.out.file_renaming.flatten().collect().ifEmpty(file("file_renaming"))

      // Demultiplexing
      demux(
        qc_se.out.filtered,
        tag_validation.out.fasta,
        ch_biosamples_sym, 
        ch_biosamples_asym,
        ch_file_renaming)

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
        demux_illumina_notmerged(
          merge_pe.out.nm,
          prep_barcodes.out.barcodesm)

        // Channel of non-merged reads by sample (split into sample tuples)
        // ch_R1 = demux_illumina_notmerged.out.demux_pe....

        // Non-merged sample list
        ch_nonmerged = demux_illumina_notmerged.out.samples_nonm_pe.splitText().map{it -> it.trim()}

        // Trim primers of nonmerged PE reads
        // Estimate sequence qualities
        // Dereplicate R1 and R2 independently
        // trim_primers_pe(demux_illumina_notmerged.out.demux_pe.flatten())

        // Join nonmerged reads with poly-N pads
        join_pe(
          ch_nonmerged,
          demux_illumina_notmerged.out.demux_pe.flatten().collect()   // all non-merged R1 and R2 files
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
    
        // Reference-based chimera removal
        ch_chimerabd = Channel.value(params.chimera_db)
        chimera_ref(homopolymer.out.hc, ch_chimerabd)
    
        // De novo chimera search
        chimera_denovo(homopolymer.out.hc)

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

        // Reference-based chimera removal
        ch_chimerabd = Channel.value(params.chimera_db)
        chimera_ref(just_derep.out.nhc, ch_chimerabd)

        // De novo chimera search
        chimera_denovo(just_derep.out.nhc)

      }  // end of ITS

      // --Primer-trimmed sequences are already dereplicated
      if(params.its_region == "none"){
          
        // just_derep(trim_primers.out.primertrimmed_fa)

        // Reference-based chimera removal
        ch_chimerabd = Channel.value(params.chimera_db)
        chimera_ref(trim_primers.out.primertrimmed_fa, ch_chimerabd)

        // De novo chimera search
        chimera_denovo(trim_primers.out.primertrimmed_fa)
      }

      // Assembled ITS is also primer trimmed and dereplicated
      if(params.its_region == "ITS1_5.8S_ITS2"){

        // Reference-based chimera removal
        ch_chimerabd = Channel.value(params.chimera_db)
        chimera_ref(assemble_its.out.itsnf, ch_chimerabd)

        // De novo chimera search
        chimera_denovo(assemble_its.out.itsnf)
      }

    } // end of homopolymer correction condition

    // Chimera rescue
    ch_chimerafiles = chimera_ref.out.chimeric.collect()
    chimera_rescue(ch_chimerafiles)

    // Aggregate de novo chimeras into a single file
    chimera_denovo_agg(chimera_denovo.out.denovochim.collect())

    // Create channel with filtered reads
    ch_filteredseqs = chimera_ref.out.nonchimeric
      .concat(chimera_rescue.out.rescuedchimeric)
      .collect()

    // Global dereplication
    glob_derep(ch_filteredseqs)

    // Pool sequences (for a final sequence table)
    pool_seqs(ch_filteredseqs)

    // OTU clustering
    otu_clust(glob_derep.out.globderep)

    // Create OTU table
    otu_tab(
      otu_clust.out.otus,
      pool_seqs.out.seqsnf)

    // Tag-jump removal
    tj(otu_tab.out.otutab)


    // Check optional channel with de novo chimera scores
    ch_denovoscores = chimera_denovo_agg.out.alldenovochim.ifEmpty(file('DeNovo_Chimera.txt'))

    // Create sequence table
    prep_seqtab(
      pool_seqs.out.seqtabnf,               // non-filtered sequence table
      pool_seqs.out.seqsnf,                 // Sequences in FASTA format
      otu_tab.out.samples_uc,               // sequence mapping to OTUs
      tj.out.tjs,                           // tag-jumped OTU list
      ch_denovoscores,                      // de novo chimera scores
      seq_qual.out.quals                    // sequence qualities
      )


    
    //// Read count summary
    
      // Initial data - Per-sample input channels
      if( params.demultiplexed == false ){

        if(params.seqplatform == "PacBio"){
          ch_all_demux = demux.out.samples_demux.flatten().collect()
        }

        if(params.seqplatform == "Illumina"){
          ch_all_demux = demux_illumina.out.samples_demux.flatten().collect()
        }

      } else {
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
      ch_chimref     = chimera_ref.out.chimeric.flatten().collect().ifEmpty(file("no_chimref"))
      ch_chimdenovo  = chimera_denovo.out.denovochim.flatten().collect().ifEmpty(file("no_chimdenovo"))
      ch_chimrescued = chimera_rescue.out.rescuedchimeric.flatten().collect().ifEmpty(file("no_chimrescued"))

      // Tag-jump filtering channel
      ch_tj = tj.out.tjs


      // Count reads and prepare summary stats for the run
      // Currently, implemented only for PacBio
      // For Illumina, need replace:
      //   `ch_input` -> `ch_inputR1` & `ch_inputR2`
      //   `qc_se`    -> `qc_pe`

      if(params.seqplatform == "PacBio"){

      read_counts(
          ch_input,                // input data
          qc_se.out.filtered,      // data that passed QC
          ch_all_demux,            // demultiplexed sequences per sample
          ch_all_primerchecked,    // primer-cheched sequences
          ch_all_primerartefacts,  // multiprimer artefacts
          ch_all_trim,             // ITSx-extracted or primer-trimmed sequences
          ch_homopolymers,         // Homopolymer stats
          ch_chimref,              // Reference-based chimeras
          ch_chimdenovo,           // De novo chimeras
          ch_chimrescued,          // Rescued chimeras
          ch_tj,                   // Tag-jump filtering stats
          prep_seqtab.out.seq_rd   // Final table with sequences
          )

      } // end of read_counts for PacBio


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

  }
  
  
  // Dump the software versions to a file
  software_versions_to_yaml(Channel.topic('versions'))
      .collectFile(
          storeDir: "${params.outdir}/pipeline_info",
          name:     'software_versions.yml',
          sort:     true,
          newLine:  true
      )

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

  // Demultiplexing
  demux(
    qc_se.out.filtered,
    tag_validation.out.fasta,
    ch_biosamples_sym, 
    ch_biosamples_asym,
    ch_file_renaming)

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




// On completion
workflow.onComplete {
    println "Pipeline completed at : $workflow.complete"
    println "Duration              : ${workflow.duration}"
    println "Execution status      : ${workflow.success ? 'All done!' : 'Failed' }"
}

// On error
workflow.onError {
    println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
