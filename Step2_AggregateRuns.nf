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

// SWARM clustering
params.swarm_d          = 2
params.swarm_fastidious = false
params.swarm_d1boundary = 3       // min mass of large OTUs, only for Fastidious + d=1

// Alignment parameters
// NB. vsearch scores = 2 * usearch scores  !!
// E.g., "20I/2E" = penalty 20 for opening internal gaps, and 2 for opening terminal gaps (left or right)

params.alignment_penalties = "UNITE"   // alternatively, "default"

if(params.alignment_penalties == "UNITE"){

  // Alternative dereplication as in UNITE
  // Allow query sequences vary 4% in length at 100% similarity
  params.unite_querycov  = 0.96 
  params.unite_targetcov = 0.96 

  // VSEARCH
  params.vsearch_gapopen = "0I/0E"   // penalties for gap opening     (usearch, "0.0/0.0E")
  params.vsearch_gapext  = "2I/1E"   // penalties for gap extension   (usearch, "1.0/0.5E")

}
if(params.alignment_penalties == "default"){
  
  // VSEARCH
  params.vsearch_gapopen = "20I/2E"
  params.vsearch_gapext  = "2I/1E"
}

// LULU
params.lulu = true
params.lulu_match     = 90.0    // minimum similarity threshold (default, 84.0)
params.lulu_ratio     = 1.0
params.lulu_ratiotype = "min"   // "min" or "avg"
params.lulu_relcooc   = 0.95    // relative co-occurrence (default, 0.95)


// Default thresholds for singleton and de novo chimera removal 
params.max_MEEP         = 0.5
params.max_ChimeraScore = 0.6
params.recover_lowqsingletons = true
params.recover_denovochimeras = true

// Pool and dereplicate sequences from all sequencing runs
process dereplication {

    label "main_container"

    publishDir "${params.outdir}/01.Dereplicated", mode: 'symlink'
    cpus 1

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


// Pool sequences from all sequencing runs,
// Dereplicate allowing query sequences to vary in length at 100% similarity (by default, 4% length variation allowed)
process dereplication_unite {

    label "main_container"

    publishDir "${params.outdir}/01.Dereplicated", mode: 'symlink'
    cpus 8

    input:
      path(inputs, stageAs: "?/*")

    output:
      path "Dereplicated.fa.gz", emit: derep
      path "Dereplicated.uc.gz", emit: derep_uc

    script:
    """
    echo -e "Dereplicating sequences\n"

    ## NB. by default, UNITE uses `cluster_fast`, which sorts sequences by length
    ## Here, we use `cluster_size`, which sorts by abundance

    find . -name "*.fa.gz" | parallel -j1 \
      "zcat {}" \
      | sed '/^>/ s/;sample=.*;/;/' \
      | vsearch \
        --cluster_size - \
        --id         1 \
        --iddef      0 \
        --query_cov  ${params.unite_querycov} \
        --target_cov ${params.unite_targetcov} \
        --strand both \
        --sizein --sizeout \
        --threads ${task.cpus} \
        --uc Dereplicated.uc \
        --centroids Dereplicated.fa
    
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
      --gapopen ${params.vsearch_gapopen} \
      --gapext  ${params.vsearch_gapext } \
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
      --gapopen ${params.vsearch_gapopen} \
      --gapext  ${params.vsearch_gapext } \
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



// Cluster sequences with SWARM (dynamic similarity threshold)
process cluster_swarm {

    label "main_container"

    publishDir "${params.outdir}/03.Clustered_SWARM", mode: 'symlink'
    // cpus 10

    input:
      path input

    output:
      path "SWARM_representatives.fa.gz", emit: clust
      path "SWARM.uc.gz",                 emit: clust_uc
      path "SWARM.swarms.gz",             emit: swarms
      path "SWARM.struct.gz",             emit: struct
      path "SWARM.stats.gz",              emit: stats

    script:
    fastidious = params.swarm_fastidious ? "--fastidious" : ""
    """
    echo -e "Clustering sequences with SWARM\n"
    echo -e "Note: sequences with ambiguous nucleotides will be excluded!\n"

    ## Swarm works with ACGTU alphabet only
    ## 1. So check if there are any sequences with ambiguities
    ## 2. If any, remove them
    ## 3. Cluster

    ## Count number of sequences with ambiguities (will go through the entire file)
    # AMBIGS=\$(seqkit grep --count --by-seq --use-regexp --ignore-case --pattern "[RYSWKMBDHVN]" ${input})

    ## Remove sequences with ambiguities
    zcat ${input} \
    | awk '{if (/^>/) {a = \$0} else {if (/^[ACGT]*\$/) {printf "%s\n%s\n", a, \$0}}}' \
    | swarm \
      --differences ${params.swarm_d} \
      --boundary ${params.swarm_d1boundary} \
      ${fastidious} \
      --threads ${task.cpus} \
      --usearch-abundance \
      --statistics-file    SWARM.stats \
      --internal-structure SWARM.struct \
      --uclust-file        SWARM.uc \
      --output-file        SWARM.swarms \
      --seeds              SWARM_representatives.fa

     # -r, --mothur                        output using mothur-like format

    echo -e "\n..Swarm clustering finished\n"

    ## Compress results
    echo -e "..Compressing results\n"
    parallel -j ${task.cpus} "gzip -6 {}" \
      ::: "SWARM_representatives.fa" "SWARM.uc" "SWARM.swarms" "SWARM.struct" "SWARM.stats"

    echo -e "..Done\n"
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
      path(otus_fasta)

    output:
      path "OTU_table_wide.txt.gz", emit: otutabwide
      path "OTU_table_long.txt.gz", emit: otutablong
      path "OTU_table_wide.RData",  emit: otutabwider
      path "OTU_table_long.RData",  emit: otutablongr
      path "OTUs.fa.gz",            emit: seqs

    script:
    """

    pool_seq_runs.R \
      --ucderep ${uc_derep} \
      --ucclust ${uc_clust} \
      --otus    ${otus_fasta} \
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
      --usearch_global ${sequences} \
      --db ${sequences} \
      --self \
      --id         "\$VSID" \
      --iddef      1 \
      --gapopen    ${params.vsearch_gapopen} \
      --gapext     ${params.vsearch_gapext } \
      --query_cov  0.9 \
      --userfields query+target+id \
      --maxaccepts 0 \
      --maxhits    10 \
      --threads    ${task.cpus} \
      --userout    LULU_match_list.txt

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
    derep_ch   = dereplication.out.derep
    derepuc_ch = dereplication.out.derep_uc


    // // Pool and dereplicate all sequences
    // if(params.alignment_penalties == "UNITE"){
    //   // Clustering-based dereplication, allowing for a slight length variation of sequences 
    //   dereplication_unite(ch_seqs)
    //   derep_ch   = dereplication_unite.out.derep
    //   derepuc_ch = dereplication_unite.out.derep_uc
    // }
    // if(params.alignment_penalties == "default"){
    //   // Fast, hash-based dereiplication
    //   dereplication(ch_seqs)
    //   derep_ch   = dereplication.out.derep
    //   derepuc_ch = dereplication.out.derep_uc
    // }
    // 
    // NB. In case with large number of sequences, UNITE-style dereplication is extremly slow.
    // Probably, it is possible to improve the speed, by using two steps:
    //   hash-based dereplication first, then additional round of clustering-based derep.
    // But it would add extra complexity to manage and combine two UC files.


    // Denoizing
    if ( params.unoise == true ) {
      unoise(derep_ch)
      unoize_ch = unoise.out.unoise
    } else {
      unoize_ch = derep_ch
    }

    // Clustering
    if ( params.clustering_method == "swarm" ) {
      cluster_swarm(unoize_ch)
      cluster_ch = cluster_swarm.out.clust
      clustuc_ch = cluster_swarm.out.clust_uc
    }
    if ( params.clustering_method == "vsearch" ) {
      cluster_vsearch(unoize_ch)
      cluster_ch = cluster_vsearch.out.clust
      clustuc_ch = cluster_vsearch.out.clust_uc
    }


   // Pool sequence tables and aggregate at OTU level
   ch_seqtabs = Channel.fromPath(
     params.data_path + "/**/07_SeqTable/Seqs.RData",
     checkIfExists: true).collect()

   // Summarize sequence abundances by OTU and sample
   summarize(
    ch_seqtabs,     // Step-1 sequnece tables in long format
    derepuc_ch,     // UC file with dereplication info
    clustuc_ch,     // UC file with OTU clustering info
    cluster_ch      // FASTA file with OTUs
   )

   // Post-clustering curation with LULU
   if ( params.lulu == true ) {
     lulu(
       summarize.out.otutabwide,
       summarize.out.seqs
       // cluster_ch        // In the Clustered.fa.gz, there are seqs excluded from OTU table
     )
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
