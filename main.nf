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

// Print the version and exit
if (params.version) {
  ver = "NextITS " + workflow.manifest.version
  if (workflow.commitId) { ver += " revision " + workflow.commitId.substring(0, 7) }
  println "${ver}\n"
  exit(0)
}

// Show a custom help message and exit
if (params.helpMsg){
  include { helpMsg } from './modules/help_message.nf'
  helpMsg()
  exit(0)
}


// Enable topic channels
// nextflow.preview.topic = true   // Nextflow < 25.04.0
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
${logoColors.green}                    𝗡𝗲𝘅𝘁${logoColors.purple}𝗜𝗧𝗦 ${logoColors.cyan}${workflow_version}${logoColors.reset}
${logoColors.green}        SSU     ${logoColors.purple}ITS1    ${logoColors.green}5.8S   ${logoColors.purple}ITS2     ${logoColors.green}LSU      ${logoColors.reset}
${logoColors.green}     ▒▒▒▒▒▒▒▒▒${logoColors.purple}░░░░░░░░░${logoColors.green}▒▒▒▒▒${logoColors.purple}░░░░░░░░░░${logoColors.green}▒▒▒▒▒▒▒▒▒▒▒▒${logoColors.reset}
${logoColors.dim}----------------------------------------------------${logoColors.reset}
"""

log.info logo

// Print all parameters using nf-schema plugin
// include { paramsSummaryLog } from 'plugin/nf-schema'
// log.info paramsSummaryLog(workflow)  // will print params from Step-1 and Step-2 simultaneously


// Include the pipeline initialisation subworkflow
// requires newer nf-core template and schema
// include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_NextITS_pipeline'


// Show help msg
if (params.helpMsg){
    helpMsg()
    exit(0)
}


// Additional parameter validation for Step-1
if( params.step == "Step1" ) {

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
  if (!params.chimera_db || !file(params.chimera_db).exists()) {
      println( errorMsg("Please provide the UDB file with reference sequences for chimera removal with `--chimera_db` parameter.", params.monochrome_logs))
      println( colorize("See https://Next-ITS.github.io/installation/#databases for more information.", 'red', params.monochrome_logs))
      exit(1)
  }
  if (!(params.chimera_db.toLowerCase().endsWith('.udb'))) {
      println( errorMsg("The reference database file specified with `--chimera_db` parameter must be in UDB format.", params.monochrome_logs))
      exit 1
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

}  // end of Step-1 parameter validation


// Additional parameter validation for Step-2
if( params.step == "Step2" ) {

  if (params.preclustering == "none" && params.clustering == "none"){
    println errorMsg("Pre-clustering and clustering could not be both set to 'none'", params.monochrome_logs)
    exit(1)
  }

}  // end of Step-2 parameter validation


// Run the workflow
workflow {

  // Print step-specific parameter summary
  paramSummary(workflow, params)

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
