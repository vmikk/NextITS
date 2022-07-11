#!/usr/bin/env nextflow
/*

============================================================================
  NextITS: Pipeline to process fungal ITS amplicons
============================================================================
  Version: v0.0.1
  License: Apache-2.0 
  Github : https://github.com/vmikk/NextITS
  Website: TBA
----------------------------------------------------------------------------
*/


// Enable DSL2 syntax
nextflow.enable.dsl = 2

// Pipeline version
version = '0.0.1'

// Initialize parameters, set default values
params.data_path = "${projectDir}/pipeline_data"

params.input = false
params.outdir = "${launchDir}/results"


// Demultiplexing
params.barcodes = false
params.lima_minscore = 93
params.lima_dualbarcode = true
params.lima_W = 70
params.lima_minlen = 40


// Print the parameters to the console and to the log
log.info """
    =======================================================================
    NextITS ${version}
    =======================================================================
    Input data path: ${params.input}
    Barcodes:        ${params.barcodes}
    Output path:     ${params.outdir}
    """
    .stripIndent()



// Demultiplexing
process demux {

    label "main_container"

    publishDir "${out_1_demux}", mode: 'symlink'
    // cpus 10

    input:
      path input_fastq
      path barcodes

    output:
      path "LIMA/*.fq.gz", emit: samples_demux
      path "LIMA/lima.lima.report.gz", emit: lima_report
      path "LIMA/lima.lima.counts", emit: lima_counts
      path "LIMA/lima.lima.summary", emit: lima_summary

    script:

    // By default, demultiplex with dual barcodes
    barcodetype = params.lima_dualbarcode ? "" : "--single-side"

    """
    mkdir -p LIMA
    echo -e "Input file: " ${input_fastq}
    echo -e "Barcodes: " ${barcodes}

    echo -e "\nDemultiplexing with LIMA:"

    ## Demultiplex with LIMA
    lima \
      --same --ccs \
      ${barcodetype} \
      -W ${params.lima_W} \
      --min-length ${params.lima_minlen} \
      --min-score ${params.lima_minscore} \
      --split-named \
      --num-threads ${task.cpus} \
      --log-level INFO --log-file LIMA/_log.txt \
      ${input_fastq} \
      ${barcodes} \
      LIMA/lima.fq.gz

    echo -e "..done"

    ## Rename files
    echo -e "\n..Renaming demultiplexed files"
    rename --filename \
      's/^lima.//g; s/--.*\$/.fq.gz/' \
      \$(find LIMA -name "*.fq.gz")

    ## Compress logs
    echo -e "..Compressing log file"
    gzip -8 LIMA/lima.lima.report

    echo -e "Demultiplexing finished"
    """
}

