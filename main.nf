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


// Enable topic channels
nextflow.preview.topic = true

// Include functions
include { helpMsg } from './modules/help_message.nf'


// Include workflows
// NB! `include` statements are static, meaning they are resolved at compile time rather than at runtime!
include { S1 } from './workflows/STEP1.nf'
include { S2 } from './workflows/STEP2.nf'
include { seqstats } from './workflows/STEP1.nf'


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  VALIDATE INPUTS

// Validate & print parameter summary
// NB! works only with old schema (`everit-json-schema` library doesn't support JSON Schema draft-2020-12)
if(params.step == "Step1") {
  WorkflowMain.initialise(workflow, params, log)
}

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

}  // end of Step-1 parameter validation



// Run the workflow
workflow {

  if (params.step == "Step1") {
    S1()
  }

  if (params.step == "Step2") {
    S2()
  }

  if (params.step == "seqstats") {
    quickstats()
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
