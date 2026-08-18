#!/usr/bin/env nextflow
/*
============================================================================
  NextITS: Pipeline to process eukaryotic ITS amplicons
============================================================================
  License: Apache-2.0
  Github : https://github.com/vmikk/NextITS
  Website: https://Next-ITS.github.io/
----------------------------------------------------------------------------
*/

// NB!!:
// - provide absolute paths to the input data (e.g. --input and --barcodes)
// - File names should not contain period (.) characters (except for extensions)

/*
 * Typed pipeline parameters (Nextflow 26+)
 * Config-time defaults (outdir, tracedir, primers, qc_twocolor, chunking_n, storagemode)
 * are in conf/params.config and must not be duplicated here
 */
params {
    // Config-time parameters (defaults in conf/params.config)
    outdir: String
    tracedir: String

    // Step selection
    step: String = "Step1"

    // ~~~~~~~~~~~~~~ Step-1 params

    // Inputs
    // (String, not Path: may be file or directory; chimera_db may be a bare filename)
    input: String? = null               // FASTQ file or directory
    input_R1: String? = null
    input_R2: String? = null
    barcodes: String? = null            // FASTA file

    seqplatform: String = "PacBio"      // "PacBio" or "Illumina"
    
    // ITS part selector
    its_region: String = "full"
    //   "full" = default (full-length ITS sequence, after trimming SSU and LSU regions by ITSx)
    //   "ITS1" or "ITS2"
    //   "none" = just trim primers
    //   "ITS1_5.8S_ITS2" = assemble near-full-length ITS from ITSx output (useful in the case if primers are too close to SSU or LSU, and ITSx is not able to detect full-length sequence)
    //   "SSU" or "LSU"

    // Quality control
    qc_maxee: Float? = null             // only for single-end reads
    qc_maxeerate: Float = 0.01          // only for single-end reads
    qc_maxhomopolymerlen: Integer = 25  // max len of homopolymer regions (if >=, sequence will be removed)
    qc_maxn: Integer = 4

    // QC for Illumina reads
    qc_phredmin: Integer? = null
    qc_phredperc: Integer? = null
    qc_polyglen: Integer? = null
    qc_avgphred: Integer? = null      // Only for PE reads
    qc_twocolor: Boolean              // reduced resolution Phred-scores (two-color Illumina chemistry)
    

    // Is data demultiplexed?
    // If false (default), input = 1 fastq file and 1 fasta file
    // If true,            input = multiple fastq files
    demultiplexed: Boolean

    // Demultiplexing - PacBio & LIMA
    lima_barcodetype: String = "dual_symmetric"  // "single", "dual", "dual_symmetric", "dual_asymmetric"
    lima_minscore: Integer = 93         // minimum barcode score (93 is for 12bp-long barcodes)
    lima_minendscore: Integer = 50      // only useful for asymmetric barcoding schemes with different barcodes in a pair
    lima_minrefspan: Float = 0.75       // min read span relative to the barcode length
    lima_minscoringregions: Integer = 2 // for dual barcodes only (2 = requires both barcodes)
    lima_windowsize: Integer = 70       // window size (in base pairs)
    lima_minlen: Integer = 40           // minimum sequence length after clipping
    lima_remove_unknown: Boolean        // remove unknown barcode combinations (in dual-barcoding modes)

    // Demultiplexing single-end reads with cutadapt
    barcode_window: Integer = 30
    barcode_errors: Integer = 1
    barcode_overlap: Integer = 11

    // Illumina pair-end read assembly
    pe_minoverlap: Integer = 20
    pe_difflimit: Integer = 5
    pe_diffperclimit: Integer = 20
    pe_nlimit: Integer = 10
    pe_minlen: Integer = 30

    // What to do with not merged reads (Illumina-only)
    illumina_keep_notmerged: Boolean = true
    illumina_joinpadgap: String = "NNNNNNNNNN"
    illumina_joinpadqual: String = "IIIIIIIIII"   // quality score of 40

    // Primer trimming
    primer_forward: String
    primer_reverse: String
    primer_foverlap: Integer   // if not specified, will be derived automatically
    primer_roverlap: Integer   // if not specified, will be derived automatically

    primer_mismatches: Integer = 2

    // ITSx
    ITSx_evalue: Float = 0.1
    ITSx_partial: Integer = 0           // off, otherwise specify min length cutoff for partial ITS sequences to keep
    ITSx_tax: String = "all"
    ITSx_complement: String = "F"       // "F" (check single strand) or "T" (check both DNA strands for matches to HMM-profiles)
    ITSx_heuristics: Boolean            // use HMM heuristics for ITSx extraction
    /// ITSx_singledomain = true ....  optional arguments
    ITSx_to_parquet: Boolean = true     // convert ITSx output (FASTA files) to Parquet
      ITSx_chunk_size: Integer = 10000  // chunk size (number of dereplicated sequences per sample) for distributed ITSx processing; set to 0 to disable chunking


    // Primer trimming (for Illumina)
    trim_minlen: Integer = 10

    // Homopolymer compression
    hp: Boolean = true
    hp_similarity: Float = 0.999
    hp_iddef: Float = 2

    // Which chimera removal methods to use
    chimera_methods: String = "ref,denovo"   // null or "none" also supported

    // Reference-based chimera removal
    chimera_db: String = "Eukaryome_1.9.3_241222_FullITS_100-800.udb"
    chimera_rescueoccurrence: Integer = 2

    // De novo chimera identification (UCHIME1)
    chimeranov_abskew: Float = 2.0
    chimeranov_dn: Float = 1.4
    chimeranov_mindiffs: Float = 3
    chimeranov_mindiv: Float = 0.8
    chimeranov_minh: Float = 0.28
    chimeranov_xn: Float = 8.0

    // Tag-jump removal
    tj: Boolean = true        // run tag-jump removal
    tj_id: Float = 1          // depreplicate or pre-cluster: 1 = just dereplicate, < 1 (e.g., 0.99) = cluster at 99% similarity
    tj_iddef: Integer = 2
    tj_f: Float = 0.01        // UNCROSS parameter f
    tj_p: Float = 1

    // Singleton removal
    // singleton_minrelabundance = 1   // % of sample abundance   // not implemented yet - TODO

    // Collapsing similar sequences
    // coverage



    // ~~~~~~~~~~~~~~ Step-2 params

    data_path: String = "${launchDir}/Step1_Results"
    
    // Pool sample replicates (e.g., re-sequenced samples) in the final OTU
    merge_replicates: Boolean
    
    // Filtering sequences (trimmed amplicons) by length
    ampliconlen_min: Integer? = null
    ampliconlen_max: Integer? = null
    
    // Singleton and de novo chimera removal 
    max_MEEP: Float = 0.5
    max_ChimeraScore: Float = 0.6
    recover_lowqsingletons: Boolean = true
    recover_denovochimeras: Boolean = true
    
    // Number of chunks to split the dataset into prior clustering
    chunking_n: Integer?       // number of chunks
    chunking_id: Float = 0.6   // minimum sequence identity for clustering

    // Sequence denoising or pre-clustering ("none", "unoise", "dada2", "swarm_d1", "homopolymer")
    preclustering: String = "none"
    
    // Denoising with UNOISE
    unoise_alpha: Float = 6.0
    unoise_minsize: Integer = 1

    // Denoising with DADA2
    dada2_engine: String = "papa2"             // "dada2" (original R-based implementation) or "papa2" (Python-based implementation)
    dada2_pooling: String = "global"           // "global" or "byrun" (not implemented yet)
    dada2_error_estimation: String = "shared"  // "per_bucket" or "shared" (independent model for each chunk or shared model for all chunks)
    dada2_nbases: Float = 1e8
    dada2_bandsize: Integer = 16
    dada2_detectsingletons: Boolean = true
    dada2_omegaA: Float = 1e-20                // algorithm sensitivity (sets the p-value threshold at which new ASVs are inferred; reduce to increase sensitivity at the cost of more spurious ASVs)
    dada2_omegaC: Float = 1e-40                // controls the balance between correcting errors and throwing away errors
    dada2_omegaP: Float = 1e-4
    dada2_maxconsist: Integer = 10
    dada2_match: Integer = 4
    dada2_mismatch: Integer = -5
    dada2_gappenalty: Integer = -8
    dada2_maxreadsperseq: Integer = 1000   // cap max number of reads per sequence for error rate estimation (set to 0 to disable)

    // Sequence clustering method ("none" / "vsearch" / "swarm" / "shmatching")
    clustering: String = "vsearch"
    
    // VSEARCH clustering
    otu_id: Float = 0.98
    otu_iddef: Integer = 2
    otu_qmask: String = "dust"

    // SWARM clustering
    swarm_d: Integer = 1
    swarm_fastidious: Boolean = true
    swarm_d1boundary: Integer = 3

    // Alternative dereplication as in UNITE:
    // Allow query sequences to vary 4% in length at 100% similarity
    // These are only used when alignment_penalties == "UNITE"
    alignment_penalties: String = "default"   // "default" or "UNITE"
    unite_querycov: Float = 0.96
    unite_targetcov: Float = 0.96

    // VSEARCH alignment penalties
    // "default" mode: vsearch_gapopen = "20I/2E"
    // "UNITE"   mode: vsearch_gapopen = "0I/0E"   (usearch equivalent: "0.0/0.0E")
    // Override with --vsearch_gapopen / --vsearch_gapext when using UNITE mode
    vsearch_gapopen: String = "20I/2E"   // penalties for gap opening
    vsearch_gapext: String = "2I/1E"     // penalties for gap extension

    // LULU
    lulu: Boolean = true
    lulu_match: Float = 95.0
    lulu_ratio: Float = 1.0
    lulu_ratiotype: String = "min"
    lulu_relcooc: Float = 0.95
    lulu_maxhits: Integer = 0

    // Species-hypothesis (SH) matching
    sh: Boolean = true
    sh_thresholds: String = ""
    sh_coveragevariation: Float = 0.96


    // ~~~~~~~~~~~~~~ Generic params

    gzip_compression: Integer = 7
    storagemode: String

    // help supports --help and --help <paramName>
    help: String? = null    // nf-schema-based help message
    helpFull: Boolean       // nf-schema-based full help message
    helpMsg: Boolean        // custom help message

    // nf-schema configuration
    validate_params: Boolean = true
    showHidden: Boolean     // wired to validation.help.showHiddenParameter
    
    version: Boolean
    email: String? = null
    email_on_fail: String? = null
    plaintext_email: Boolean

    // Max Job request parameters
    max_cpus: Integer = 80
    max_memory: String = "132.GB"
    max_time: String = "240.h"

    monochrome_logs: Boolean  // use ANSI color output in pipeline messages

}

// nf-schema functions for parameter validation and help
include { validateParameters; paramsHelp } from 'plugin/nf-schema'

// Include custom help message function
include { helpMsg } from './modules/help_message.nf'

// Include custom parameter summary function
include { paramSummary } from './modules/parameter_summary'

// Include color utilities
include { getColors; colorize; colorizeMultiple; errorMsg; warningMsg; infoMsg; successMsg } from './modules/colors'

// Include workflows
// NB! `include` statements are static, meaning they are resolved at compile time rather than at runtime!
include { S1 } from './workflows/STEP1.nf'
include { S2 } from './workflows/STEP2.nf'
include { seqstats } from './workflows/STEP1.nf'


// Entry workflow
workflow {

  // Schema-based help (--help, --helpFull, --help <param>)
  if (params.help) {
    if (params.help == 'true') {
      log.info paramsHelp(hideWarning: true)
    } else {
      log.info paramsHelp(params.help, hideWarning: true)
    }
    exit 0
  }

  if (params.helpFull) {
    log.info paramsHelp(fullHelp: true, hideWarning: true)
    exit 0
  }

  // Print the version and exit
  if (params.version) {
    def ver = "NextITS " + workflow.manifest.version
    if (workflow.commitId) { ver += " revision " + workflow.commitId.substring(0, 7) }
    println "${ver}\n"
    exit(0)
  }

  // Show a custom help message and exit
  if (params.helpMsg) {
    helpMsg()
    exit(0)
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  VALIDATE INPUTS

  // Print NextITS logo
  def logoColors = getColors(params.monochrome_logs)

  def workflow_version = workflow.manifest.version ?: "unknown"
  if (workflow.commitId) { workflow_version += " (${workflow.commitId.substring(0, 7)})" }

  def logo = """
${logoColors.dim}----------------------------------------------------${logoColors.reset}
                    ${colorizeMultiple("Next", ['green', 'bold'], params.monochrome_logs)}${colorizeMultiple("ITS", ['purple', 'bold'], params.monochrome_logs)} ${logoColors.cyan}${workflow_version}${logoColors.reset}
${logoColors.green}        SSU     ${logoColors.purple}ITS1    ${logoColors.green}5.8S   ${logoColors.purple}ITS2     ${logoColors.green}LSU      ${logoColors.reset}
${logoColors.green}     ▒▒▒▒▒▒▒▒▒${logoColors.purple}░░░░░░░░░${logoColors.green}▒▒▒▒▒${logoColors.purple}░░░░░░░░░░${logoColors.green}▒▒▒▒▒▒▒▒▒▒▒▒${logoColors.reset}
${logoColors.dim}----------------------------------------------------${logoColors.reset}
"""

  log.info logo

  // Print all parameters using nf-schema plugin
  // include { paramsSummaryLog } from 'plugin/nf-schema'
  // log.info paramsSummaryLog(workflow)  // will print params from Step-1 and Step-2 simultaneously

  // Print step-specific parameter summary
  paramSummary(workflow, params)
  validateParameters()

  def is_demultiplexed = false
  def run_hp = true
  def keep_notmerged = false

  // Additional runtime parameter validation
  // These checks are performed after schema validation and handle
  // conditional logic and file existence checks that cannot be expressed in JSON Schema

  // Additional parameter validation for Step-1
  if (params.step == "Step1" || params.step == "seqstats") {
    is_demultiplexed = params.demultiplexed

    if (!params.input && params.seqplatform == "PacBio") {
        println( errorMsg("Please provide the input file with sequences in FASTQ.gz or BAM format with `--input` parameter.", params.monochrome_logs))
        exit(1)
    }
    if (!params.input_R1 && !params.input_R2 && params.seqplatform == "Illumina") {
        println( errorMsg("Please provide input files with sequences in FASTQ.gz format with `--input_R1` and `--input_R2` parameters.", params.monochrome_logs))
        exit(1)
    }
    if (!params.barcodes && !is_demultiplexed) {
        println( errorMsg("Please provide the file with sample barcodes in FASTA format with `--barcodes` parameter.", params.monochrome_logs))
        exit(1)
    }
  }

  if (params.step == "Step1") {
    run_hp = params.hp
    keep_notmerged = params.illumina_keep_notmerged

    // Reference-based chimera removal
    if (params.chimera_methods && params.chimera_methods.toLowerCase().split(',').contains('ref')) {
      if (!params.chimera_db || !file(params.chimera_db).exists()) {
          println( errorMsg("For reference-based chimera removal, please provide the database in UDB format with `--chimera_db` parameter.", params.monochrome_logs))
          println( colorize("       See https://Next-ITS.github.io/installation/#databases for more information.", 'red', params.monochrome_logs))
          println( colorize("Alternatively, you can disable reference-based chimera removal with `--chimera_methods` parameter (set it to `none` or `denovo`).", 'red', params.monochrome_logs))
          exit(1)
      }
      if (!(params.chimera_db.toLowerCase().endsWith('.udb'))) {
          println( errorMsg("The reference database file specified with `--chimera_db` parameter must be in UDB format.", params.monochrome_logs))
          println( colorize("       See https://Next-ITS.github.io/installation/#databases for more information.", 'red', params.monochrome_logs))
          exit(1)
      }
    }

    if (run_hp && params.seqplatform == "Illumina" && keep_notmerged) {
        println( errorMsg("Homopolymer compression is not implemented for Illumina non-merged reads (add `--hp false` to your command).", params.monochrome_logs))
        exit(1)
    }
    if (params.seqplatform == "Illumina" && is_demultiplexed) {
        println( errorMsg("Handling demultiplexed data for Illumina is not implemented yet.", params.monochrome_logs))
        exit(1)
    }

    if (params.seqplatform == "Illumina" && keep_notmerged && params.its_region != "none") {
        println( warningMsg("Unmerged Illumina reads are not compatible with ITSx. Amplicons will be primer-trimmed.", params.monochrome_logs))
    }

    // ITSx profiles validation
    if (params.its_region != "none") {

      /*
      Currently, the following regex pattern is used to pre-validate the `ITSx_tax` parameter (in schema):
      "^(?:all|
      (?:alveolata|bryophyta|bacillariophyta|amoebozoa|euglenozoa|fungi|chlorophyta|rhodophyta|phaeophyceae|marchantiophyta|metazoa|oomycota|haptophyceae|raphidophyceae|rhizaria|synurophyceae|tracheophyta|eustigmatophyceae|apusozoa|parabasalia)
      (?:,\\s*(?:alveolata|bryophyta|bacillariophyta|amoebozoa|euglenozoa|fungi|chlorophyta|rhodophyta|phaeophyceae|marchantiophyta|metazoa|oomycota|haptophyceae|raphidophyceae|rhizaria|synurophyceae|tracheophyta|eustigmatophyceae|apusozoa|parabasalia))*)$"

      this forbids:
      - mixing `all` with other values
      - empty elements and trailing commas
      - invalid values
      */

      def itsx_profiles = params.ITSx_tax

      // `ITSx_tax` must be a non-empty string (if specifying `--ITSx_tax ""`, Nextflow may coerce empty/flag to boolean)
      if (itsx_profiles == null || itsx_profiles instanceof Boolean) {
        println( errorMsg("Parameter --ITSx_tax must have a value (e.g. 'all' or 'fungi,rhizaria').", params.monochrome_logs) )
        exit(1)
      }
      if (itsx_profiles.toString().trim().isEmpty()) {
        println( errorMsg("Parameter --ITSx_tax cannot be empty. Use 'all' or a comma-separated list of taxa.", params.monochrome_logs) )
        exit(1)
      }

      // Allowed profiles
      def ITSX_ALLOWED = [
          'alveolata','bryophyta','bacillariophyta','amoebozoa','euglenozoa','fungi',
          'chlorophyta','rhodophyta','phaeophyceae','marchantiophyta','metazoa','oomycota',
          'haptophyceae','raphidophyceae','rhizaria','synurophyceae','tracheophyta',
          'eustigmatophyceae','apusozoa','parabasalia'
      ] as Set

      // Parse the specified profile string
      def itsx_items = itsx_profiles.toString().split(',', -1) as List<String>

      // Empty-item validation (empty or whitespace-only tokens, incl. ",," and trailing commas)
      def emptyIdx = []
      itsx_items.eachWithIndex { s, i ->
        if (s == null || s.trim().isEmpty()) emptyIdx << i
      }
      if (emptyIdx) {
        println( errorMsg("Parameter --ITSx_tax: empty entries are not allowed (check commas at positions: ${emptyIdx.join(', ')}).", params.monochrome_logs) )
        exit(1)
      }

      // Disallow internal whitespaces
      def whitespaces = itsx_items.findAll { s ->
        def tr = s.toString().trim()
        !(tr ==~ /\S+/)   // after trimming, token must be all non-whitespace
      }
      if (whitespaces) {
          println( errorMsg("Parameter --ITSx_tax: whitespace is not allowed in profile names.", params.monochrome_logs) )
          exit(1)
      }

      // Detect duplicates
      itsx_items = itsx_items.collect { s -> s.trim() }
      def dups = itsx_items.countBy { s -> s }.findAll { k, v -> v > 1 }.keySet().toList()
      if (dups) {
        println( errorMsg("Parameter --ITSx_tax: duplicated profile names are not allowed: ${dups.join(', ')}", params.monochrome_logs) )
        exit(1)
      }

      // Disallow mixing 'all' with specific profile names
      if (itsx_items.size() > 1 && itsx_items.contains('all')) {
        println( errorMsg("Parameter --ITSx_tax: do not combine 'all' with taxon-specific profile names.", params.monochrome_logs))
        exit(1)
      }

      // Validate values against the allow-list (skip when it's exactly ['all'])
      if (!(itsx_items.size() == 1 && itsx_items[0] == 'all')) {
        def invalid_profiles = (itsx_items as Set) - ITSX_ALLOWED
        if (invalid_profiles) {
            println( errorMsg("Parameter --ITSx_tax: invalid profile names - ${invalid_profiles.join(', ')}", params.monochrome_logs) )
            println( colorize("       Supported profiles: `all` OR a comma-separated list of the following: ${ITSX_ALLOWED.join(', ')}", 'red', params.monochrome_logs))
            exit(1)
        }
      }

      // Currently, there is no X.hmm profile (Apusozoa)
      if (itsx_items.contains('apusozoa')) {
        println( errorMsg("Parameter --ITSx_tax: Apusozoa profile is not yet supported in ITSx.", params.monochrome_logs))
        exit(1)
      }

    } // end of ITSx profiles validation

  }  // end of Step-1 parameter validation


  // Additional parameter validation for Step-2
  if (params.step == "Step2") {
    if (params.preclustering == "none" && params.clustering == "none" && params.lulu) {
      println errorMsg("LULU can not be applied when pre-clustering and clustering are set to 'none'", params.monochrome_logs)
      exit(1)
    }

    if (params.preclustering == "dada2" && params.dada2_pooling == "byrun" &&
        params.chunking_n != null && params.chunking_n > 1) {
      println errorMsg("By-sequencing-run pooling in DADA2 is not compatible with chunking.", params.monochrome_logs)
      println( colorize("Set `--chunking_n` to 1 to disable chunking OR use `--dada2_pooling global`.", 'red', params.monochrome_logs))
      exit(1)
    }

  }  // end of Step-2 parameter validation


  // Run the selected step
  if (params.step == "Step1") {
    S1()
  }

  if (params.step == "Step2") {
    S2()
  }

  if (params.step == "seqstats") {
    seqstats()
  }

  workflow.onComplete = {
    def completion_time = workflow.complete.format(java.time.format.DateTimeFormatter.ofPattern("yyyy-MM-dd'T'HH:mm:ssXXX"))

    println "Pipeline completed at : ${completion_time}"
    println "Duration              : ${workflow.duration}"
    println "CPU hours             : ${workflow.stats.getComputeTimeFmt()}"
    println "Execution status      : ${workflow.success ? 'All done!' : 'Failed' }"
  }

  workflow.onError = {
    println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
  }

}
