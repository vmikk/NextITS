#!/usr/bin/env nextflow
/*

============================================================================
  NextITS: Pipeline to process fungal ITS amplicons
  Step II: aggregate sequencing runs
============================================================================
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

include { software_versions_to_yaml } from '../modules/version_parser.nf'



// Aggregate sequences from all sequencing runs, remove de novo chimeras
process aggregate_sequences {

    label "main_container"

    // cpus 6

    input:
      path(inputs, stageAs: "?/*")

    output:
      path "Seqs.fa.gz",   emit: seqs
      path "Seqs.parquet", emit: seqs_parquet

    script:
    """
    echo -e "Aggregating sequences\n"
    
    aggregate_sequences.R \
      --seqtabs       . \
      --maxchim       ${params.max_ChimeraScore} \
      --recoverdenovo ${params.recover_denovochimeras} \
      --output        Seqs \
      --threads       ${task.cpus}

    """
}

// Pool and dereplicate sequences from all sequencing runs
process dereplication {

    label "main_container"

    publishDir "${params.outdir}/01.Dereplicated", mode: 'symlink'
    // cpus 8

    input:
      path seqs

    output:
      path "Dereplicated.fa.gz", emit: derep
      path "Dereplicated.uc.gz", emit: derep_uc

    script:
    def minlen = params.ampliconlen_min ? "--minseqlength ${params.ampliconlen_min}" : ""
    def maxlen = params.ampliconlen_max ? "--maxseqlength ${params.ampliconlen_max}" : ""

    // Calculate optimal number of threads
    def maxPigzThreads = 8  // Maximum threads per pigz instance
    def totalCPUs = task.cpus
    
    // Try to maximize CPUs per pigz while ensuring full CPU utilization
    def pigzCPUs = Math.min(maxPigzThreads, Math.ceil(Math.sqrt(totalCPUs * 2)).intValue())
    def parallelJobs = Math.max(1, Math.floor(totalCPUs / pigzCPUs).intValue())
    
    // Recalculate pigzCPUs to use all available CPUs
    pigzCPUs = Math.min(maxPigzThreads, Math.floor(totalCPUs / parallelJobs).intValue())

    """
    echo -e "Dereplicating sequences\n"

    vsearch \
      --derep_fulllength ${seqs} \
      --output Dereplicated.fa \
      --strand both \
      ${minlen} ${maxlen} \
      --fasta_width 0 \
      --threads 1 \
      --sizein --sizeout \
      --uc Dereplicated.uc

    echo -e "..Dereplication finished\n"

    ## Compress results
    echo -e "\nCompressing results"
    parallel -j ${parallelJobs} \
      "pigz -p ${pigzCPUs} -${params.gzip_compression} {}" \
      ::: "Dereplicated.uc" "Dereplicated.fa"

    """
}


// Pool sequences from all sequencing runs,
// Dereplicate allowing query sequences to vary in length at 100% similarity (by default, 4% length variation allowed)
process dereplication_unite {

    label "main_container"

    publishDir "${params.outdir}/01.Dereplicated", mode: 'symlink'
    // cpus 8

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
        --iddef      2 \
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
    parallel -j ${task.cpus} "gzip -${params.gzip_compression} {}" \
      ::: "Dereplicated.uc" "Dereplicated.fa"

    """
}


// Homopolymer correction (global, for pooled and dereplicated data)
process homopolymer {

    label "main_container"

    publishDir "${params.outdir}/02.Homopolymer", mode: "${params.storagemode}"
    // cpus 1

    input:
      path input

    output:
      path "HomopolymerCompressed.fa.gz", emit: hp
      path "HomopolymerCompressed.uc.gz", emit: hp_uc

    script:
    """
    ## Run homopolyer correction globally

    echo -e "Running homopolymer correction"

    echo -e "\nCompressing repeats"
    zcat ${input} \
      | homopolymer_compression.sh \
      | gzip -2 \
      > homo_compressed.fa.gz

    echo -e "\nAdditional dereplication"
    vsearch \
      --derep_fulllength homo_compressed.fa.gz \
      --output - \
      --strand both \
      --fasta_width 0 \
      --threads 1 \
      --sizein --sizeout \
      --uc HomopolymerCompressed.uc \
    > homo_compressed_dereplicated.fa

    ## Substitute homopolymer-comressed sequences with uncompressed ones
    ## (update size annotaions)
    echo -e "\nExtracting representative sequences"

    seqkit fx2tab ${input} > inp_tab.txt
    seqkit fx2tab homo_compressed_dereplicated.fa > clust_tab.txt

    if [ -s inp_tab.txt ]; then
      substitute_compressed_seqs.R \
        inp_tab.txt clust_tab.txt \
        HomopolymerCompressed_tmp.fa

      echo -e "..Done"
    else
      echo -e "..Input data looks empty, nothing to proceed with"
    fi


    ## Sort by number of reads
    vsearch \
      --sortbysize HomopolymerCompressed_tmp.fa \
      --sizein --sizeout \
      --threads ${task.cpus} \
      --fasta_width 0 \
      --output - \
      | gzip -${params.gzip_compression} \
      > HomopolymerCompressed.fa.gz


    #### combine_derep_and_hpcorrection.R

    echo -e "\nHomopolymer correction finished\n"

    ## Compress results
    echo -e "\nCompressing results"
    gzip -${params.gzip_compression} HomopolymerCompressed.uc
    
    ## Remove temporary files
    echo -e "\nRemoving temporary files"
    rm homo_compressed.fa.gz
    rm homo_compressed_dereplicated.fa
    rm HomopolymerCompressed_tmp.fa
    rm inp_tab.txt
    rm clust_tab.txt
    """
}


// Denoize sequences with UNOISE
process unoise {

    label "main_container"

    publishDir "${params.outdir}/02.UNOISE", mode: 'symlink'
    // cpus 8

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
    parallel -j 1 \
      "pigz -p ${task.cpus} -${params.gzip_compression} {}" \
      ::: "UNOISE.fa" "UNOISE.uc"

    """
}




// Preclustering with SWARM and d1
process precluster_swarm {

    label "main_container"

    publishDir "${params.outdir}/02.Preclustered_SWARM_d1", mode: 'symlink'
    // cpus 8

    input:
      path input

    output:
      path "SWARM_representatives.fa.gz", emit: clust
      path "SWARM.uc.gz",                 emit: clust_uc
      path "SWARM.swarms.gz",             emit: swarms
      path "SWARM.struct.gz",             emit: struct
      path "SWARM.stats.gz",              emit: stats

    script:
    """
    echo -e "Pre-clustering sequences with SWARM d=1\n"
    echo -e "Note: sequences with ambiguous nucleotides will be excluded!\n"

    ## Remove sequences with ambiguities
    zcat ${input} \
    | awk '{if (/^>/) {a = \$0} else {if (/^[ACGT]*\$/) {printf "%s\\n%s\\n", a, \$0}}}' \
    | swarm \
      --differences 1 \
      --boundary    ${params.swarm_d1boundary} \
      --fastidious \
      --threads     ${task.cpus} \
      --usearch-abundance \
      --statistics-file    SWARM.stats \
      --internal-structure SWARM.struct \
      --uclust-file        SWARM.uc \
      --seeds              SWARM_representatives.fa \
      > SWARM.swarms

    echo -e "\n..Swarm pre-clustering finished\n"

    ## Compress results
    echo -e "..Compressing results\n"
    parallel -j 1 \
      "pigz -p ${task.cpus} -${params.gzip_compression} {}" \
      ::: "SWARM_representatives.fa" "SWARM.uc" "SWARM.swarms" "SWARM.struct" "SWARM.stats"

    echo -e "..Done\n"
    """
}



// Cluster sequences with VSEARCH (fixed similarity threshold)
process cluster_vsearch {

    label "main_container"

    publishDir "${params.outdir}/03.Clustered_VSEARCH", mode: 'symlink'
    // cpus 8

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
    | gzip -${params.gzip_compression} > Clustered.fa.gz
    
    echo -e "..Done"

    ## Compress UC file
    echo -e "\nCompressing UC file"
    pigz -p ${task.cpus} -${params.gzip_compression} Clustered.uc

    """
}



// Cluster sequences with SWARM (dynamic similarity threshold)
process cluster_swarm {

    label "main_container"

    publishDir "${params.outdir}/03.Clustered_SWARM", mode: 'symlink'
    // cpus 8

    input:
      path input

    output:
      path "SWARM_representatives.fa.gz", emit: clust
      path "SWARM.uc.gz",                 emit: clust_uc
      path "SWARM.swarms.gz",             emit: swarms
      path "SWARM.struct.gz",             emit: struct
      path "SWARM.stats.gz",              emit: stats

    exec:
      fastidious = (params.swarm_fastidious.toBoolean() == true & params.swarm_d.toInteger() == 1) ? "--fastidious --boundary ${params.swarm_d1boundary}" : ""
      println("swarm_fastidious: ${params.swarm_fastidious}, swarm_d: ${params.swarm_d}")
      println("fastid option:  ${fastidious}")

    script:
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
    | awk '{if (/^>/) {a = \$0} else {if (/^[ACGT]*\$/) {printf "%s\\n%s\\n", a, \$0}}}' \
    | swarm \
        --differences ${params.swarm_d} \
        ${fastidious} \
        --threads     ${task.cpus} \
        --usearch-abundance \
        --statistics-file    SWARM.stats \
        --internal-structure SWARM.struct \
        --uclust-file        SWARM.uc \
        --seeds              SWARM_representatives.fa \
        > SWARM.swarms

    # --output-file SWARM.swarms  # to avoid buffering, it's better to stream data into a file (with >)
    # -r, --mothur                # output using mothur-like format

    echo -e "\n..Swarm clustering finished\n"

    ## Compress results
    echo -e "..Compressing results\n"
    parallel -j 1 \
      "pigz -p ${task.cpus} -${params.gzip_compression} {}" \
      ::: "SWARM_representatives.fa" "SWARM.uc" "SWARM.swarms" "SWARM.struct" "SWARM.stats"

    echo -e "..Done\n"
    """
}



// Merge UC files
process merge_uc {

    label "main_container"
    
    publishDir "${params.outdir}/04.PooledResults", mode: 'symlink'
    // cpus 4

    input:
      path(uc_derep)
      path(uc_preclust)
      path(uc_clust)

    output:
      path "UC_Pooled.parquet", emit: uc

    script:
    """
    echo -e "merge_uc process\n"

    ## Parse UC files from different steps, convert to parquet format
    echo -e "..Parsing dereplicated UC file\n"
    ucs --input ${uc_derep} --output UC_derep.parquet

    if [ -f ${uc_preclust} ] && [ "${uc_preclust}" != "NoPrecluster" ]; then
      echo -e "..Parsing pre-clustered UC file\n"
      ucs --input ${uc_preclust} --output UC_preclust.parquet
      UCPRECLUST="UC_preclust.parquet"
    else
      UCPRECLUST="NoPrecluster"
    fi

    if [ -f ${uc_clust} ]; then
      echo -e "..Parsing clustered UC file\n"
      ucs --input ${uc_clust} --output UC_clust.parquet
    fi

    ## Merge UC files into a single file
    echo -e "..Merging UC files\n"
    merge_uc_files.R \
      --ucderep    UC_derep.parquet \
      --ucpreclust \${UCPRECLUST} \
      --ucclust    UC_clust.parquet \
      --output     UC_Pooled.parquet

    """
}


// Summarize sequence abundance by OTU
process summarize {

    label "main_container"

    publishDir "${params.outdir}/04.PooledResults", mode: 'symlink'
    // cpus 4

    input:
      path(seqtab)
      path(uc_parquet)
      path(otus_fasta)

    output:
      path "OTU_table_wide.txt.gz", emit: otutabwide
      path "OTU_table_long.txt.gz", emit: otutablong
      path "OTU_table_wide.RData",  emit: otutabwider
      path "OTU_table_long.RData",  emit: otutablongr
      path "OTUs.fa.gz",            emit: seqs

    script:
    """

    summarize_clustered_data.R \
      --seqtab         ${seqtab} \
      --uc             ${uc_parquet} \
      --otus           ${otus_fasta} \
      --maxmeep        ${params.max_MEEP} \
      --recoversinglet ${params.recover_lowqsingletons} \
      --mergesamples   ${params.merge_replicates} \
      --threads        ${task.cpus}

    """
}


// Post-clustering curation
process lulu {

    label "main_container"

    publishDir "${params.outdir}/05.LULU", mode: 'symlink'
    // cpus 8

    input:
      path otu_table
      path sequences

    output:
      path "OTU_table_LULU.txt.gz",          emit: lulu
      path "LULU_match_list.txt.gz",         emit: matches
      path "LULU_merging_statistics.txt.gz", emit: stats
      path "OTUs_LULU.fa.gz",                emit: fasta

    script:
    """
    echo -e "Post-clustering curation with MUMU (C++ implementation of LULU)\n"

    ## If Clustered.fa.gz used as input
    ## (but there are sequences excluded from the OTU table)
    # echo -e "Removing size annotations from sequence headers"
    # zcat ${sequences} \
    #   | sed -r '/^>/ s/;size=[0-9]+//g' \
    #   | gzip -${params.gzip_compression} > tmp_sequences.fa.gz


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
      --db             ${sequences} \
      --self \
      --id         "\$VSID" \
      --iddef      1 \
      --gapopen    ${params.vsearch_gapopen} \
      --gapext     ${params.vsearch_gapext } \
      --query_cov  0.9 \
      --userfields query+target+id \
      --maxaccepts 0 \
      --maxhits    ${params.lulu_maxhits} \
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
      --log LULU_merging_statistics.txt \
      --threads                      ${task.cpus} \
      --minimum_match                ${params.lulu_match} \
      --minimum_ratio                ${params.lulu_ratio} \
      --minimum_ratio_type           ${params.lulu_ratiotype} \
      --minimum_relative_cooccurence ${params.lulu_relcooc}

    echo -e "..Compressing LULU-curated OTU table\n"
    parallel -j 1 \
      "pigz -p ${task.cpus} -${params.gzip_compression} {}" \
      ::: "OTU_table_LULU.txt" "LULU_merging_statistics.txt" "LULU_match_list.txt"

    echo -e "..LULU done\n"

    echo -e "\nPreparing sequence subset\n"

    echo -e "..Extracting OTU IDs\n"
    zcat OTU_table_LULU.txt.gz \
      | awk 'NR > 1 {print \$1}' \
      > curated_OTU_ids.txt

    echo -e "..Extracting sequences\n"
    rg -z -A 1 \
      -f curated_OTU_ids.txt \
      --context-separator "" \
      --threads ${task.cpus} \
      ${sequences} \
    | sed '/^\$/d' \
    | gzip -${params.gzip_compression} \
    > OTUs_LULU.fa.gz

    ## Remove temporary files
    echo -e "\nAll done!\n"
    echo -e "Removing temporary files\n"
    # rm tmp_sequences.fa.gz
    rm tmp_OTU_table.txt curated_OTU_ids.txt

    """
}

// LULU merging statistics format:
// 1.  `query_otu_name` - name of query OTU
// 2.  `parent_otu_name` - name of potential parent OTU
// 3.  `similarity_pct` - percentage of similarity (0 to 100)
// 4.  `query_total_abundance` - total abundance of the query OTU (sum through all samples)
// 5.  `parent_total_abundance` - total abundance of the potential parent OTU (sum through all samples)
// 6.  `query_overlap_abundance` - overlap abundance of the query OTU (sum through all samples where the potential parent OTU is also present)
// 7.  `parent_overlap_abundance` - overlap abundance of the potential parent OTU (sum through all samples where the query OTU is also present)
// 8.  `query_incidence` - incidence of the query OTU (number of samples where the query OTU is present)
// 9.  `parent_incidence` - incidence of the potential parent OTU (number of samples where the potential parent OTU is present)
// 10. `both_incidence` - incidence of the potential parent OTU (number of samples where both the potential parent OTU and the query OTU are present)
// 11. `smallest_abundance_ratio` - smallest abundance ratio (for each sample, compute the abundance of the potential parent OTU divided by the abundance of the query OTU)
// 12. `sum_abundance_ratios` - sum of the abundance ratios
// 13. `avg_abundance_ratio` - average value of abundance ratios
// 14. `smallest_nonnull_ratio` - smallest non-null abundance ratio (exclude ratios for samples where the query OTU is present but not the potential parent OTU)
// 15. `avg_nonnull_ratio` - average value of non-null abundance ratios (exclude ratios for samples where the query OTU is present but not the potential parent OTU)
// 16. `largest_ratio` - largest ratio value
// 17. `relative_cooccurrence` - relative co-occurence value (number of samples where both the potential parent OTU and the query OTU are present divided by the number of samples where the query OTU is present)
// 18. `status` - status: 'accepted' or 'rejected'
//     The potential parent OTU is either accepted as a parent, or rejected

// Step-2 workflow
workflow S2 {

    // Find quality-filtered sequence tables
    ch_seqtabs = Channel.fromPath(
      params.data_path + "/**/07_SeqTable/Seqs.parquet",
      checkIfExists: true).collect()

    // Aggregate sequences, remove de novo chimeras
    aggregate_sequences(ch_seqtabs)

    // Pool and dereplicate all sequences
    dereplication(aggregate_sequences.out.seqs)
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


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Pre-clustering / denoising
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Denoizing
    if ( params.preclustering == "none" ) {
      denoise_ch    = derep_ch
      preclustuc_ch = file('NoPrecluster')
    
    // Denoise with UNOISE
    } else if ( params.preclustering == "unoise" ) {
      unoise(derep_ch)
      denoise_ch    = unoise.out.unoise
      preclustuc_ch = unoise.out.unoise_uc
    

    // Precluster with SWARM
    } else if ( params.preclustering == "swarm_d1" ){
      precluster_swarm(derep_ch)
      denoise_ch    = precluster_swarm.out.clust
      preclustuc_ch = precluster_swarm.out.clust_uc
    }


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Clustering
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */


    // Greedy clustering with VSEARCH
    if ( params.clustering == "vsearch" ) {
      cluster_vsearch(denoise_ch)
      cluster_ch = cluster_vsearch.out.clust
      clustuc_ch = cluster_vsearch.out.clust_uc
    
    // Clustering with SWARM
    } else if ( params.clustering == "swarm" ) {
      
      // If pre-clustering was already done with the same d, just take the previous results
      if(params.preclustering == "swarm_d1" & params.swarm_d == 1){
        cluster_ch = precluster_swarm.out.clust
        clustuc_ch = precluster_swarm.out.clust_uc
        preclustuc_ch = file('NoPrecluster')
      
      // Otherwise, run SWARM
      } else {
        cluster_swarm(denoise_ch)
        cluster_ch = cluster_swarm.out.clust
        clustuc_ch = cluster_swarm.out.clust_uc
      }

    // Do not cluster, use zOTUs from UNOISE
    } else if ( params.preclustering == "unoise" & params.clustering == "none" ) {
      cluster_ch = unoise.out.unoise
      clustuc_ch = unoise.out.unoise_uc
      preclustuc_ch = file('NoPrecluster')
    
    } else if ( params.preclustering == "none" & params.clustering == "none" ){
      println "No pre-clustering or clustering was done"

      // TODO: create table based on dereplicated sequences?
    }


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Result processing
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

 
    // Pool UC files
    merge_uc(
     derepuc_ch,       // UC file with dereplication info
     preclustuc_ch,    // UC file with pre-clustering or denoising (optional)
     clustuc_ch        // UC file with OTU clustering info
    )

    // Summarize sequence abundances by OTU and sample
    summarize(
     aggregate_sequences.out.seqs_parquet,  // Step-1 sequnece tables in long format with de novo chimeras removed
     merge_uc.out.uc,                       // Combined UC files with sequence membership info
     cluster_ch                             // FASTA file with OTUs
    )

    // Post-clustering curation with LULU
    if ( params.lulu == true ) {
      lulu(
        summarize.out.otutabwide,
        summarize.out.seqs
        // cluster_ch        // In the Clustered.fa.gz, there are seqs excluded from OTU table
      )
    }

    // Run statistics
    // run_summary()


  // Auto documentation of analysis procedures
  // document_analysis_step2


}
