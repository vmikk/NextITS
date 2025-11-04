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

// Databases:
//  - UDB for chimera identification


// Enable DSL2 syntax
nextflow.enable.dsl = 2

// Print the version and exit
if (params.version) {
  ver = "NextITS " + workflow.manifest.version
  if (workflow.commitId) { ver += " revision " + workflow.commitId.substring(0, 7) }
  println "${ver}\n"
  exit(0)
}

// Note: nf-schema plugin handles --help automatically via configuration in nextflow.config

// Show a custom help message and exit
if (params.helpMsg){
  include { helpMsg } from './modules/help_message.nf'
  helpMsg()
  exit(0)
}


// Enable topic channels
// nextflow.preview.topic = true   // Nextflow < 25.04.0



// nf-schema functions for parameter validation
include { validateParameters } from 'plugin/nf-schema'

// Include custom parameter summary function
include { paramSummary } from './modules/parameter_summary'

// Include color utilities
include { getColors; colorize; colorizeMultiple; errorMsg; warningMsg; infoMsg; successMsg } from './modules/colors'

// Include workflows
// NB! `include` statements are static, meaning they are resolved at compile time rather than at runtime!
include { S1 } from './workflows/STEP1.nf'
include { S2 } from './workflows/STEP2.nf'
include { seqstats } from './workflows/STEP1.nf'


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


// Additional runtime parameter validation 
// These checks are performed after schema validation and handle 
// conditional logic and file existence checks that cannot be expressed in JSON Schema

// Additional parameter validation for Step-1
if( params.step == "Step1" || params.step == "seqstats" ) {

  if (params.input == false && params.seqplatform == "PacBio") {
      println( errorMsg("Please provide the input file with sequences in FASTQ.gz or BAM format with `--input` parameter.", params.monochrome_logs))
      exit(1)
  }
  if (params.input_R1 == false && params.input_R2 == false && params.seqplatform == "Illumina") {
      println( errorMsg("Please provide input files with sequences in FASTQ.gz format with `--input_R1` and `--input_R2` parameters.", params.monochrome_logs))
      exit(1)
  }
  if (params.barcodes == false && params.demultiplexed == false) {
      println( errorMsg("Please provide the file with sample barcodes in FASTA format with `--barcodes` parameter.", params.monochrome_logs))
      exit(1)
  }
}

if( params.step == "Step1" ) {

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
        exit 1
    }
  }

  if (params.hp == true && params.seqplatform == "Illumina" && params.illumina_keep_notmerged == true) {
      println( errorMsg("Homopolymer compression is not implemented for Illumina non-merged reads.", params.monochrome_logs))
      exit(1)
  }
  if (params.seqplatform == "Illumina" && params.demultiplexed == true) {
      println( errorMsg("Handling demultiplexed data for Illumina is not implemented yet.", params.monochrome_logs))
      exit(1)
  }

  if (params.seqplatform == "Illumina" && params.illumina_keep_notmerged == true && params.its_region != "none") {
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
    itsx_items = itsx_items.collect { it.trim() }
    def dups = itsx_items.countBy { it }.findAll { k, v -> v > 1 }.keySet().toList()
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
          println( colorize("       Supported profiles: ${ITSX_ALLOWED.join(', ')}", 'red', params.monochrome_logs))
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
if( params.step == "Step2" ) {

  if (params.preclustering == "none" && params.clustering == "none"){
    println errorMsg("Pre-clustering and clustering could not be both set to 'none'", params.monochrome_logs)
    exit(1)
  }

  if (params.preclustering == "dada2" && params.dada2_pooling == "byrun" && 
      (params.chunking_n > 1 || params.chunking_n != null)){
    println errorMsg("By-sequencing-run pooling in DADA2 is not compatible with chunking.", params.monochrome_logs)
    println( colorize("Set `--chunking_n` to 1 to disable chunking OR use `--dada2_pooling global`.", 'red', params.monochrome_logs))
    exit(1)
  }


}  // end of Step-2 parameter validation


// Run the workflow
workflow {

  // Print step-specific parameter summary
  paramSummary(workflow, params)
  validateParameters()

  if (params.step == "Step1") {
    S1()
  }

  if (params.step == "Step2") {
    S2()
  }

  if (params.step == "seqstats") {
    seqstats()
  }

}


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
