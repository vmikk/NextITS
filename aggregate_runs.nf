#!/usr/bin/env nextflow
/*

============================================================================
  NextITS: Pipeline to process fungal ITS amplicons
  Step II: aggregate sequencing runs
============================================================================
  Version: v0.2.0
  License: Apache-2.0
  Github : https://github.com/vmikk/NextITS
  Website: https://next-its.github.io/
----------------------------------------------------------------------------
*/

// Enable DSL2 syntax
nextflow.enable.dsl = 2

// Pipeline version
version = '0.2.0'

params.outdir = "Step2"


// Initialize parameters, set default values
params.data_path = "${projectDir}/Output"

// Pool and dereplicate sequences from all sequencing runs
process dereplication {

    label "main_container"

    // cpus 1

    input:
      path(inputs, stageAs: "?/*")

    output:
      path "Dereplicated.fa.gz", emit: derep
      path "Dereplicated.uc.gz", emit: derep_uc

    script:
    """
    echo -e "\nDereplicating sequences"

    find . -name "*.fa.gz" | parallel -j1 \
      "zcat {}" \
      | sed '/^>/ s/;sample=.*;/;/' \
      | vsearch \
        --derep_fulllength - \
        --output - \
        --strand both \
        --fasta_width 0 \
        --threads 1 \
        --sizein --sizeout \
        --uc Dereplicated.uc \
      | gzip -7 > Dereplicated.fa.gz
    
    echo -e "..Done"

    ## Compress UC file
    echo -e "\nCompressing UC file"
    gzip -7 Dereplicated.uc

    """
}




//  The default workflow
workflow {

    // Input files = FASTA files from individual sequencing runs
    // e.g. "*/07_ASV_table/ASVs.fa.gz"
    
    ch_seqs = Channel.fromPath(
      params.data_path + "/**/07_ASV_table/ASVs.fa.gz",
      checkIfExists: true).collect()

    // Pool and dereplicate all sequences
    dereplication(ch_seqs)

}

