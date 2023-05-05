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
version = '0.4.0'

params.outdir = "Step2"


// Path to the Step-1 results
params.data_path = "${projectDir}/Output"

// Denoising
params.unoise         = false
params.unoise_alpha   = 2.0
params.unoise_minsize = 8

// Sequence clustering method (VSEARCH / SWARM)
params.clustering_method = "vsearch"

// VSEARCH clustering
params.otu_id    = 0.98
params.otu_iddef = 2        // also for UNOISE
params.otu_qmask = "dust"   // also for UNOISE
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



// Cluster sequences with VSEARCH (fixed similarity threshold)
process cluster_vsearch {

    label "main_container"

    publishDir "${params.outdir}/03.Clustered_VSEARCH", mode: 'symlink'
    // cpus 10

    input:
      path input

    output:
      path "Clustered.fa.gz", emit: clust
      path "Clustered.uc.gz", emit: clust_uc

    script:
    """
    echo -e "Clustering sequences with VSEARCH\n"

    vsearch \
      --cluster_size ${input} \
      --id      ${params.otu_id} \
      --iddef   ${params.otu_iddef} \
      --qmask   ${params.otu_qmask} \
      --threads ${task.cpus} \
      --sizein --sizeout \
      --strand both \
      --fasta_width 0 \
      --uc Clustered.uc \
      --centroids - \
    | gzip -7 > Clustered.fa.gz
    
    echo -e "..Done"

    ## Compress UC file
    echo -e "\nCompressing UC file"
    gzip -7 Clustered.uc

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
      cluster_vsearch(unoize_ch)
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
