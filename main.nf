#!/usr/bin/env nextflow


// Enable DSL2 syntax
nextflow.enable.dsl = 2

// Pipeline version
version = '0.0.1'

// Initialize parameters, set default values
params.data_path = "${projectDir}/pipeline_data"

params.input = false
params.outdir = "${launchDir}/results"
