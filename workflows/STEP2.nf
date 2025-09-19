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
// - (optionally) Denoize with UNOISE or DADA2
// - (optionally) Cluster:
//   * SWARM
//   * VSEARCH
// - LULU (via MUMU implementation)
// - Prepare OTU table (wide, aggregate sequence abundance by ASV/OTU/Swarm cluster)



// Enable DSL2 syntax
nextflow.enable.dsl = 2

include { software_versions_to_yaml } from '../modules/version_parser.nf'
include { CLUSTERING } from '../subworkflows/clustering_subworkflow.nf'





// Aggregate sequences from all sequencing runs, remove de novo chimeras
process aggregate_sequences {

    label "main_container"

    // cpus 6

    input:
      path(inputs, stageAs: "?/*")

    output:
      path "Seqs.fa.gz",   emit: seqs
      path "Seqs.parquet", emit: seqs_parquet
      tuple val("${task.process}"), val('R'), eval('Rscript -e "cat(R.version.string)" | sed "s/R version //"'),  topic: versions
      tuple val("${task.process}"), val('data.table'), eval('Rscript -e "cat(as.character(packageVersion(\'data.table\')))"'),  topic: versions
      tuple val("${task.process}"), val('arrow'), eval('Rscript -e "cat(as.character(packageVersion(\'arrow\')))"'),  topic: versions
      tuple val("${task.process}"), val('Biostrings'), eval('Rscript -e "cat(as.character(packageVersion(\'Biostrings\')))"'),  topic: versions

    script:
    """
    echo -e "Aggregating sequences\\n"
    
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
      tuple val("${task.process}"), val('vsearch'), eval('vsearch --version 2>&1 | head -n 1 | sed "s/vsearch //g" | sed "s/,.*//g" | sed "s/^v//" | sed "s/_.*//"'), topic: versions

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
    echo -e "Dereplicating sequences\\n"

    vsearch \
      --derep_fulllength ${seqs} \
      --output Dereplicated.fa \
      --strand both \
      ${minlen} ${maxlen} \
      --fasta_width 0 \
      --threads 1 \
      --sizein --sizeout \
      --uc Dereplicated.uc

    echo -e "..Dereplication finished\\n"

    ## Compress results
    echo -e "\\nCompressing results"
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
    echo -e "Dereplicating sequences\\n"

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
    echo -e "\\nCompressing results"
    parallel -j ${task.cpus} "gzip -${params.gzip_compression} {}" \
      ::: "Dereplicated.uc" "Dereplicated.fa"

    """
}


// Fast pre-clustering of the dataset (to split into chunks prior processing)
process linclust {

    label "main_container"

    input:
      path input

    output:
      path "DB_clu.tsv", emit: db_clu
      tuple val("${task.process}"), val('mmseqs'), eval('mmseqs version'), topic: versions

    script:
    """

    ## Create DB
    echo -e "..DB creation\\n"
  
    mmseqs createdb \
      --dbtype 2 \
      --createdb-mode 0 \
      --shuffle 0 \
      ${input} \
      mmseqs_db


    ## Run (cascaded) clustering
    echo -e "..Lin-Clustering\\n"

    mmseqs linclust \
      mmseqs_db \
      linclusters_db \
      tmplc \
      --min-seq-id ${params.chunking_id} \
      --cluster-mode 0 \
      --similarity-type 2 \
      -c 0.7 --cov-mode 0 \
      -k 15 \
      --kmer-per-seq 100 \
      --kmer-per-seq-scale 0.3 \
      --spaced-kmer-mode 0 \
      --mask 0 \
      --split-memory-limit 100G \
      --remove-tmp-files 1 \
      --threads ${task.cpus}

    ## Generate a TSV-formatted output of clustering
    echo -e "..Generating TSV-formatted output of clustering\\n"

    mmseqs createtsv \
      mmseqs_db mmseqs_db \
      linclusters_db \
      DB_clu.tsv \
      --threads ${task.cpus}

    """
}

// Bucketize sequences into clusters
process bucketize {

    label "main_container"

    input:
      path sequences
      path clusters

    output:
      path "bucket_*.fa.gz", emit: buckets
      tuple val("${task.process}"), val('R'), eval('Rscript -e "cat(R.version.string)" | sed "s/R version //"'),  topic: versions
      tuple val("${task.process}"), val('data.table'), eval('Rscript -e "cat(as.character(packageVersion(\'data.table\')))"'),  topic: versions
      tuple val("${task.process}"), val('Biostrings'), eval('Rscript -e "cat(as.character(packageVersion(\'Biostrings\')))"'),  topic: versions

    script:
    numchunks = params.chunking_n ? "--numbuckets ${params.chunking_n}" : ""
    """
    echo -e "..Bucketizing sequences\\n"

    bucketize_db.R \
      --db      ${clusters} \
      --fasta   ${sequences} \
      ${numchunks} \
      --summary bucket_summary.txt \
      --threads ${task.cpus}
  
  """
}





// Merge processed buckets (e.g., clustered sequences)
process merge_buckets {

    label "main_container"

    // Conditional publishing to match non-chunked directory structure
    // (only enabled when chunking is used)
    
    // Pre-clustering results (if any)
    publishDir "${params.outdir}/02.Homopolymer", 
               mode: 'symlink',
               enabled: (params.chunking_n != null && params.chunking_n >= 2) && params.preclustering == "homopolymer",
               pattern: "PreClustered.uc.gz",
               saveAs: { filename -> filename == "PreClustered.uc.gz" ? "HomopolymerCompressed.uc.gz" : null }
    
    publishDir "${params.outdir}/02.UNOISE", 
               mode: 'symlink',
               enabled: (params.chunking_n != null && params.chunking_n >= 2) && params.preclustering == "unoise",
               pattern: "PreClustered.uc.gz",
               saveAs: { filename -> filename == "PreClustered.uc.gz" ? "UNOISE.uc.gz" : null }
               
    publishDir "${params.outdir}/02.DADA2", 
               mode: 'symlink',
               enabled: (params.chunking_n != null && params.chunking_n >= 2) && params.preclustering == "dada2",
               pattern: "PreClustered.uc.gz",
               saveAs: { filename -> filename == "PreClustered.uc.gz" ? "DADA2_denoised.uc.gz" : null }
               
    publishDir "${params.outdir}/02.Preclustered_SWARM_d1", 
               mode: 'symlink',
               enabled: (params.chunking_n != null && params.chunking_n >= 2) && params.preclustering == "swarm_d1",
               pattern: "PreClustered.uc.gz",
               saveAs: { filename -> filename == "PreClustered.uc.gz" ? "SWARM.uc.gz" : null }
    
    // Final clustering results - publish to clustering directory if clustering != "none"
    publishDir "${params.outdir}/03.Clustered_VSEARCH", 
               mode: 'symlink',
               enabled: (params.chunking_n != null && params.chunking_n >= 2) && params.clustering == "vsearch",
               pattern: "Clustered.{fa,uc}.gz"
               // No saveAs needed - files already have correct names for VSEARCH
               
    publishDir "${params.outdir}/03.Clustered_SWARM", 
               mode: 'symlink',
               enabled: (params.chunking_n != null && params.chunking_n >= 2) && params.clustering == "swarm",
               pattern: "Clustered.{fa,uc}.gz",
               saveAs: { filename -> 
                   switch(filename) {
                       case "Clustered.fa.gz": return "SWARM_representatives.fa.gz"
                       case "Clustered.uc.gz": return "SWARM.uc.gz"
                       default: return null
                   }
               }
    
    // Final results when no clustering is done (clustering == "none") - publish to preclustering directory
    publishDir "${params.outdir}/02.Homopolymer", 
               mode: 'symlink',
               enabled: (params.chunking_n != null && params.chunking_n >= 2) && params.clustering == "none" && params.preclustering == "homopolymer",
               pattern: "Clustered.{fa,uc}.gz",
               saveAs: { filename -> 
                   switch(filename) {
                       case "Clustered.fa.gz": return "HomopolymerCompressed.fa.gz"
                       case "Clustered.uc.gz": return "HomopolymerCompressed.uc.gz"
                       default: return null
                   }
               }
               
    publishDir "${params.outdir}/02.UNOISE", 
               mode: 'symlink',
               enabled: (params.chunking_n != null && params.chunking_n >= 2) && params.clustering == "none" && params.preclustering == "unoise",
               pattern: "Clustered.{fa,uc}.gz",
               saveAs: { filename -> 
                   switch(filename) {
                       case "Clustered.fa.gz": return "UNOISE.fa.gz"
                       case "Clustered.uc.gz": return "UNOISE.uc.gz"
                       default: return null
                   }
               }
               
    publishDir "${params.outdir}/02.DADA2", 
               mode: 'symlink',
               enabled: (params.chunking_n != null && params.chunking_n >= 2) && params.clustering == "none" && params.preclustering == "dada2",
               pattern: "Clustered.{fa,uc}.gz",
               saveAs: { filename -> 
                   switch(filename) {
                       case "Clustered.fa.gz": return "DADA2_denoised.fa.gz"
                       case "Clustered.uc.gz": return "DADA2_denoised.uc.gz"
                       default: return null
                   }
               }
               
    publishDir "${params.outdir}/02.Preclustered_SWARM_d1", 
               mode: 'symlink',
               enabled: (params.chunking_n != null && params.chunking_n >= 2) && params.clustering == "none" && params.preclustering == "swarm_d1",
               pattern: "Clustered.{fa,uc}.gz",
               saveAs: { filename -> 
                   switch(filename) {
                       case "Clustered.fa.gz": return "SWARM_representatives.fa.gz"
                       case "Clustered.uc.gz": return "SWARM.uc.gz"
                       default: return null
                   }
               }


    // Since there are name collisions, we need to stage files with unique names
    input:
      path(preclustuc_chunks, stageAs: "pre/?/*")  // UC files for pre-clustering (optional)
      path(cluster_chunks, stageAs:    "cls/?/*")  // Sequence representatives
      path(clustuc_chunks, stageAs:    "ucs/?/*")  // UC files for clustering

    output:
      path "PreClustered.uc.gz", emit: preclustuc_ch, optional: true
      path "Clustered.fa.gz",    emit: cluster_ch
      path "Clustered.uc.gz",    emit: clustuc_ch


    script:
    """
    echo -e "Merging buckets\\n"

    ## Pool sequence representatives
    echo -e "..Pooling sequence representatives\\n"
    find cls -name "*.fa.gz" \
      | parallel -j 1 "cat {}" \
      > Clustered.fa.gz

    ## Pool UC files
    echo -e "..Pooling UC files\\n"
    find ucs -name "*.uc.gz" \
      | parallel -j 1 "cat {}" \
      > Clustered.uc.gz

    ## Check if pre-clustering was performed
    if [[ -e pre/1/NoPrecluster || -L "pre/1/NoPrecluster" ]]; then
      echo -e "..Pre-clustering was not performed. Skipping pooling these data\\n"
    else
      echo -e "..Pre-clustering was performed\\n"
      find pre -name "*.uc.gz" \
        | parallel -j 1 "cat {}" \
        > PreClustered.uc.gz
    fi

    echo -e "..Done\\n"
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
      tuple val("${task.process}"), val('ucs'), eval('ucs --version | sed "s/ucs //"'), topic: versions
      tuple val("${task.process}"), val('R'), eval('Rscript -e "cat(R.version.string)" | sed "s/R version //"'),  topic: versions
      tuple val("${task.process}"), val('duckdb'), eval('Rscript -e "cat(as.character(packageVersion(\'duckdb\')))"'),  topic: versions

    script:
    """
    echo -e "Merging UC files\\n"

    ## Parse UC files from different steps, convert to parquet format
    echo -e "..Parsing dereplicated UC file\\n"
    ucs --input ${uc_derep} --output UC_derep.parquet

    if [ -f ${uc_preclust} ] && [ "${uc_preclust}" != "NoPrecluster" ]; then
      echo -e "..Parsing pre-clustered UC file\\n"
      ucs --input ${uc_preclust} --output UC_preclust.parquet
      UCPRECLUST="UC_preclust.parquet"
    else
      UCPRECLUST="NoPrecluster"
    fi

    if [ -f ${uc_clust} ]; then
      echo -e "..Parsing clustered UC file\\n"
      ucs --input ${uc_clust} --output UC_clust.parquet
    fi

    ## Merge UC files into a single file
    echo -e "..Merging UC files\\n"
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
      tuple val("${task.process}"), val('R'), eval('Rscript -e "cat(R.version.string)" | sed "s/R version //"'),  topic: versions
      tuple val("${task.process}"), val('data.table'), eval('Rscript -e "cat(as.character(packageVersion(\'data.table\')))"'),  topic: versions
      tuple val("${task.process}"), val('arrow'), eval('Rscript -e "cat(as.character(packageVersion(\'arrow\')))"'),  topic: versions
      tuple val("${task.process}"), val('Biostrings'), eval('Rscript -e "cat(as.character(packageVersion(\'Biostrings\')))"'),  topic: versions

    script:
    """
    echo -e "Summarizing clustered data\\n"

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
      tuple val("${task.process}"), val('mumu'), eval('mumu --version | head -n 1 | sed "s/mumu //"'), topic: versions
      tuple val("${task.process}"), val('vsearch'), eval('vsearch --version 2>&1 | head -n 1 | sed "s/vsearch //g" | sed "s/,.*//g" | sed "s/^v//" | sed "s/_.*//"'), topic: versions
      tuple val("${task.process}"), val('ripgrep'), eval('rg --version | head -1 | sed "s/ripgrep //"'), topic: versions

    script:
    """
    echo -e "Post-clustering curation with MUMU (C++ implementation of LULU)\\n"

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
    echo -e "Preparing match list\\n"
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


    echo -e "\\nUnpacking OTU table\\n"
    gunzip --stdout ${otu_table} > tmp_OTU_table.txt

    echo -e "\\nRunning MUMU\\n"
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

    echo -e "..Compressing LULU-curated OTU table\\n"
    parallel -j 1 \
      "pigz -p ${task.cpus} -${params.gzip_compression} {}" \
      ::: "OTU_table_LULU.txt" "LULU_merging_statistics.txt" "LULU_match_list.txt"

    echo -e "..LULU done\\n"

    echo -e "\\nPreparing sequence subset\\n"

    echo -e "..Extracting OTU IDs\\n"
    zcat OTU_table_LULU.txt.gz \
      | awk 'NR > 1 {print \$1}' \
      > curated_OTU_ids.txt

    echo -e "..Extracting sequences\\n"
    rg -z -A 1 \
      -f curated_OTU_ids.txt \
      --context-separator "" \
      --threads ${task.cpus} \
      ${sequences} \
    | sed '/^\$/d' \
    | gzip -${params.gzip_compression} \
    > OTUs_LULU.fa.gz

    ## Remove temporary files
    echo -e "\\nAll done!\\n"
    echo -e "Removing temporary files\\n"
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
        Clustering / pre-clustering / denoising   with optional chunking
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // No chunking (process all sequences at once)
    if(params.chunking_n == null || params.chunking_n < 2){

      CLUSTERING(derep_ch)

      preclustuc_ch = CLUSTERING.out.preclustuc_ch
      cluster_ch    = CLUSTERING.out.cluster_ch
      clustuc_ch    = CLUSTERING.out.clustuc_ch

    } else {
    // Chunking (process sequences in N chunks)

      // Groupd sequences into clusters
      linclust(derep_ch)

      // Bucketize sequence clusters into chunks
      bucketize(derep_ch, linclust.out.db_clu)
      buckets_ch = bucketize.out.buckets.flatten()

      // Run clustering/pre-clustering/denoising subworkflow
      CLUSTERING(buckets_ch)

      // collect UC and FASTA files from all chunks
      preclustuc_chunks = CLUSTERING.out.preclustuc_ch.collect()
      cluster_chunks    = CLUSTERING.out.cluster_ch.collect()
      clustuc_chunks    = CLUSTERING.out.clustuc_ch.collect()

      // Merge buckets into a single file
      merge_buckets(
        preclustuc_chunks,
        cluster_chunks,
        clustuc_chunks)

      cluster_ch    = merge_buckets.out.cluster_ch
      clustuc_ch    = merge_buckets.out.clustuc_ch
      preclustuc_ch = merge_buckets.out.preclustuc_ch.ifEmpty(file('NoPrecluster'))

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


  // Dump the software versions to a file
  software_versions_to_yaml(Channel.topic('versions'))
      .collectFile(
          storeDir: "${params.tracedir}",
          name:     'software_versions.yml',
          sort:     true,
          newLine:  true
      )

  // Auto documentation of analysis procedures
  // document_analysis_step2


}
