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

// Step-2 workflow:
// - Dereplicate sequences
// - (optionally) Denoize with UNOISE
// - Cluster:
//   * SWARM
//   * VSEARCH
// - LULU (via MUMU implementation)
// - Prepare OTU table (wide, aggregate sequence abundance by OTU/Swarm cluster)

// Enable DSL2 syntax
nextflow.enable.dsl = 2

// Pipeline version
version = '0.2.0'

params.outdir = "Step2"


// Path to the Step-1 results
params.data_path = "${projectDir}/Output"

// Denoising
params.unoise         = false
params.unoise_alpha   = 2.0
params.unoise_minsize = 8
// Pool and dereplicate sequences from all sequencing runs
process dereplication {

    label "main_container"

    publishDir "${params.outdir}/01.Dereplicated", mode: 'symlink'
    // cpus 1

    input:
      path(inputs, stageAs: "?/*")

    output:
      path "Dereplicated.fa.gz", emit: derep
      path "Dereplicated.uc.gz", emit: derep_uc

    script:
    """
    echo -e "Dereplicating sequences\n"

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
      > Dereplicated.fa
    
    echo -e "..Dereplication finished"

    ## Compress results
    echo -e "\nCompressing results"
    parallel -j ${task.cpus} "gzip -7 {}" \
      ::: "Dereplicated.uc" "Dereplicated.fa"

    """
}


// Denoize sequences with UNOISE
process unoise {

    label "main_container"

    publishDir "${params.outdir}/02.UNOISE", mode: 'symlink'
    // cpus 10

    input:
      path input

    output:
      path "UNOISE.fa.gz", emit: unoise
      path "UNOISE.uc.gz", emit: unoise_uc

    script:
    """
    echo -e "Denoizing sequences with UNOISE\n"

    vsearch \
      --cluster_unoise ${input} \
      --unoise_alpha   ${params.unoise_alpha} \
      --minsize ${params.unoise_minsize} \
      --iddef   ${params.otu_iddef} \
      --qmask   ${params.otu_qmask} \
      --threads ${task.cpus} \
      --fasta_width 0 \
      --sizein --sizeout \
      --centroids UNOISE.fa \
      --uc UNOISE.uc

    echo -e "..UNOISE done\n"

    ## Compress results
    echo -e "\nCompressing UNOISE results"
    parallel -j ${task.cpus} "gzip -7 {}" \
      ::: "UNOISE.fa" "UNOISE.uc"

    """
}

    
    echo -e "..Done"

    ## Compress UC file
    echo -e "\nCompressing UC file"
    gzip -7 Dereplicated.uc

    """
}




//  The default workflow
workflow {

    // Input files = FASTA files from individual sequencing runs
    // e.g. "*/07_SeqTable/Seqs.fa.gz"
    
    ch_seqs = Channel.fromPath(
      params.data_path + "/**/07_SeqTable/Seqs.fa.gz",
      checkIfExists: true).collect()

    // Pool and dereplicate all sequences
    dereplication(ch_seqs)

    // Denoizing
    if ( params.unoise == true ) {
      unoise(dereplication.out.derep)
      unoize_ch = unoise.out.unoise
    } else {
      unoize_ch = dereplication.out.derep
    }
}

