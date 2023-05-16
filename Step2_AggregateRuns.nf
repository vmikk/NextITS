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

// Default thresholds for singleton and de novo chimera removal 
params.max_MEEP         = 0.5
params.max_ChimeraScore = 0.6
params.recover_lowqsingletons = true
params.recover_denovochimeras = true

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





// Summarize sequence abundance by OTU
process summarize {

    label "main_container"

    publishDir "${params.outdir}/04.PooledResults", mode: 'symlink'
    // cpus 4

    input:
      path(seqtabs, stageAs: "?/*")
      path(uc_derep)
      path(uc_clust)

    output:
      path "OTU_table_wide.txt.gz", emit: otutabwide
      path "OTU_table_long.txt.gz", emit: otutablong
      path "OTU_table_wide.RData",  emit: otutabwider
      path "OTU_table_long.RData",  emit: otutablongr
      path "OTUs.fa.gz",            emit: seqs

    script:
    """

    pool_seq_runs.R \
      --ucderep "Dereplicated.uc.gz" \
      --ucclust "Clustered.uc.gz" \
      --maxmeep ${params.max_MEEP} \
      --maxchim ${params.max_ChimeraScore} \
      --recoverdenovo  ${params.recover_lowqsingletons} \
      --recoversinglet ${params.recover_denovochimeras} \
      --threads ${task.cpus} \

    """
}


// Post-clustering curation
process lulu {

    label "main_container"

    publishDir "${params.outdir}/05.LULU", mode: 'symlink'
    // cpus 10

    input:
      path otu_table
      path sequences

    output:
      path "OTU_table_LULU.txt.gz",  emit: lulu
      path "LULU_match_list.txt.gz", emit: matches

    script:
    """
    echo -e "Post-clustering curation with MUMU (C++ implementation of LULU)\n"

    ## If Clustered.fa.gz used as input
    ## (but there are sequences excluded from the OTU table)
    # echo -e "Removing size annotations from sequence headers"
    # zcat ${sequences} \
    #   | sed -r '/^>/ s/;size=[0-9]+//g' \
    #   | gzip -4 > tmp_sequences.fa.gz


    ## MUMU similarity threshold is specified as % (e.g., 84.0)
    ## while VSEARCH requires a value in 0-1 range (e.g., 0.84)
    
    ## With bc
    # VSID=\$(echo "scale=4; x = ${params.lulu_match} / 100; if(x<1) print 0; x" | bc)
    
    ## With awk
    VSID=\$(awk -v a=${params.lulu_match} 'BEGIN { print(a/100) }')

    echo -e "VSEARCH similarity threshold: " "\$VSID"

    ## Prepare match list (+ remove size annotations)
    echo -e "Preparing match list\n"
    vsearch \
      --usearch_global ${sequences}  \
      --db ${sequences}  \
      --self  \
      --id "\$VSID" \
      --iddef 1 \
      --userfields query+target+id \
      --maxaccepts 0 \
      --query_cov  0.9 \
      --maxhits 10 \
      --threads ${task.cpus} \
      --userout LULU_match_list.txt


    # Input otu_table  = tab-separated, samples in columns
    # Input match_list = tab-separated, OTU pairwise similarity scores


    echo -e "\nUnpacking OTU table\n"
    gunzip --stdout ${otu_table} > tmp_OTU_table.txt

    echo -e "\nRunning MUMU\n"
    mumu \
      --otu_table     tmp_OTU_table.txt \
      --match_list    LULU_match_list.txt \
      --new_otu_table OTU_table_LULU.txt \
      --log lulu.log \
      --threads                      ${task.cpus} \
      --minimum_match                ${params.lulu_match} \
      --minimum_ratio                ${params.lulu_ratio} \
      --minimum_ratio_type           ${params.lulu_ratiotype} \
      --minimum_relative_cooccurence ${params.lulu_relcooc}

    echo -e "..Compressing LULU-curated OTU table\n"
    parallel -j ${task.cpus} "gzip -7 {}" \
      ::: "OTU_table_LULU.txt" "lulu.log" "LULU_match_list.txt"

    echo -e "..LULU done\n"

    ## Remove temporary files
    # rm tmp_sequences.fa.gz
    rm tmp_OTU_table.txt

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
      cluster_ch = cluster_vsearch.out.clust
      clustuc_ch = cluster_vsearch.out.clust_uc

   // Pool sequence tables and aggregate at OTU level
   ch_seqtabs = Channel.fromPath(
     params.data_path + "/**/07_SeqTable/Seqs.RData",
     checkIfExists: true).collect()

   // Summarize sequence abundances by OTU and sample
   summarize(
    ch_seqtabs,
    dereplication.out.derep_uc,
    clustuc_ch
   )

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
