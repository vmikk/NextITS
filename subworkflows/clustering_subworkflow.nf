/*
============================================================================
  NextITS: Pipeline to process eukaryotic ITS amplicons
============================================================================
  License: Apache-2.0
  Github : https://github.com/vmikk/NextITS
  Website: https://Next-ITS.github.io/
----------------------------------------------------------------------------
*/

// Subworkflow for clustering sequences (with optional pre-clustering or denoising)


// Homopolymer correction (global, for pooled and dereplicated data)
process homopolymer {

    label "main_container"

    publishDir(
      "${params.outdir}/02.Homopolymer",
      mode: "${params.storagemode}", 
      enabled: params.chunking_n == null || params.chunking_n < 2
    )

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

    echo -e "\\nCompressing repeats"
    zcat ${input} \
      | homopolymer_compression.sh \
      | gzip -2 \
      > homo_compressed.fa.gz

    echo -e "\\nAdditional dereplication"
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
    echo -e "\\nExtracting representative sequences"

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

    echo -e "\\nHomopolymer correction finished\\n"

    ## Compress results
    echo -e "\\nCompressing results"
    gzip -${params.gzip_compression} HomopolymerCompressed.uc
    
    ## Remove temporary files
    echo -e "\\nRemoving temporary files"
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

    publishDir(
      "${params.outdir}/02.UNOISE",
      mode: "${params.storagemode}",
      enabled: params.chunking_n == null || params.chunking_n < 2
    )

    // cpus 8

    input:
      path input

    output:
      path "UNOISE.fa.gz", emit: unoise
      path "UNOISE.uc.gz", emit: unoise_uc
      tuple val("${task.process}"), val('vsearch'), eval('vsearch --version 2>&1 | head -n 1 | sed "s/vsearch //g" | sed "s/,.*//g" | sed "s/^v//" | sed "s/_.*//"'), topic: versions

    script:
    """
    echo -e "Denoizing sequences with UNOISE\\n"

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

    echo -e "..UNOISE done\\n"

    ## Compress results
    echo -e "\\nCompressing UNOISE results"
    parallel -j 1 \
      "pigz -p ${task.cpus} -${params.gzip_compression} {}" \
      ::: "UNOISE.fa" "UNOISE.uc"

    """
}



// Denoize sequences with DADA2
process dada2 {

    label "main_container"

    publishDir(
      "${params.outdir}/02.DADA2",
      mode: "${params.storagemode}",
      enabled: params.chunking_n == null || params.chunking_n < 2
    )

    // cpus 8

    input:
      path input

    output:
      path "DADA2_denoised.fa.gz",        emit: dada
      path "DADA2_denoised.uc.gz",        emit: dada_uc
      path "DADA2_UC.qs",                 emit: dada_ucr
      path "DADA2_denoising_summary.txt", emit: dada_summary
      // path "DADA2_ErrorRates_noqualErrfun.RData"
      // path "DADA2_InferedSeqs_noqualErrfun.RData"
      tuple val("${task.process}"), val('R'), eval('Rscript -e "cat(R.version.string)" | sed "s/R version //" | cut -d" " -f1'),  topic: versions
      tuple val("${task.process}"), val('dada2'), eval('Rscript -e "cat(as.character(packageVersion(\'dada2\')))"'),  topic: versions
      tuple val("${task.process}"), val('data.table'), eval('Rscript -e "cat(as.character(packageVersion(\'data.table\')))"'),  topic: versions

    script:
    """
    echo -e "Denoizing sequences with DADA2\\n"

    ## DADA2 works with ACGT alphabet only
    ## 1. So check if there are any sequences with ambiguities
    ## 2. If any, remove them
    ## 3. Sort by sequence abundance
    ## 4. Convert FASTA to pseudo-FASTQ
    ## 5. Denoise

    ## Remove sequences with ambiguities
    echo -e "..Preparing sequences\\n"
    zcat ${input} \
      | awk '{if (/^>/) {a = \$0} else {if (/^[ACGT]*\$/) {printf "%s\\n%s\\n", a, \$0}}}' \
      | vsearch --sortbysize - --output - --fasta_width 0 \
      | awk 'BEGIN {RS = ">" ; FS = "\\n"} NR > 1 {print "@"\$1"\\n"\$2"\\n+"\$1"\\n"gensub(/./, "I", "g", \$2)}' \
      | gzip -${params.gzip_compression} > no_ambigs.fq.gz

    echo -e "\\n\\n..Running DADA2\\n"
    dada2_no_quals.R \
      --input            no_ambigs.fq.gz \
      --nbases           ${params.dada2_nbases} \
      --bandsize         ${params.dada2_bandsize} \
      --detectsingletons ${params.dada2_detectsingletons} \
      --omegaA           ${params.dada2_omegaA} \
      --omegaC           ${params.dada2_omegaC} \
      --omegaP           ${params.dada2_omegaP} \
      --maxconsist       ${params.dada2_maxconsist} \
      --match            ${params.dada2_match} \
      --mismatch         ${params.dada2_mismatch} \
      --gappenalty       ${params.dada2_gappenalty} \
      --threads          ${task.cpus}

    echo -e "..Denoizing with DADA2 finished\\n"
    """
}




// Preclustering with SWARM and d1
process precluster_swarm {

    label "main_container"

    publishDir(
      "${params.outdir}/02.Preclustered_SWARM_d1",
      mode: "${params.storagemode}",
      enabled: params.chunking_n == null || params.chunking_n < 2
    )

    // cpus 8

    input:
      path input

    output:
      path "SWARM_representatives.fa.gz", emit: clust
      path "SWARM.uc.gz",                 emit: clust_uc
      path "SWARM.swarms.gz",             emit: swarms
      path "SWARM.struct.gz",             emit: struct
      path "SWARM.stats.gz",              emit: stats
      tuple val("${task.process}"), val('swarm'), eval('swarm --version 2>&1 | head -n 1 | sed "s/Swarm //"'), topic: versions

    script:
    """
    echo -e "Pre-clustering sequences with SWARM d=1\\n"
    echo -e "Note: sequences with ambiguous nucleotides will be excluded!\\n"

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

    echo -e "\\n..Swarm pre-clustering finished\\n"

    ## Compress results
    echo -e "\\n..Compressing results\\n"
    parallel -j 1 \
      "pigz -p ${task.cpus} -${params.gzip_compression} {}" \
      ::: "SWARM_representatives.fa" "SWARM.uc" "SWARM.swarms" "SWARM.struct" "SWARM.stats"

    echo -e "..Done\\n"
    """
}



// Cluster sequences with VSEARCH (fixed similarity threshold)
process cluster_vsearch {

    label "main_container"

    publishDir(
      "${params.outdir}/03.Clustered_VSEARCH",
      mode: "${params.storagemode}",
      enabled: params.chunking_n == null || params.chunking_n < 2
    )

    // cpus 8

    input:
      path input

    output:
      path "Clustered.fa.gz", emit: clust
      path "Clustered.uc.gz", emit: clust_uc
      tuple val("${task.process}"), val('vsearch'), eval('vsearch --version 2>&1 | head -n 1 | sed "s/vsearch //g" | sed "s/,.*//g" | sed "s/^v//" | sed "s/_.*//"'), topic: versions

    script:
    """
    echo -e "Clustering sequences with VSEARCH\\n"

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
    echo -e "\\nCompressing UC file"
    pigz -p ${task.cpus} -${params.gzip_compression} Clustered.uc

    """
}



// Cluster sequences with SWARM (dynamic similarity threshold)
process cluster_swarm {

    label "main_container"

    publishDir(
      "${params.outdir}/03.Clustered_SWARM",
      mode: "${params.storagemode}",
      enabled: params.chunking_n == null || params.chunking_n < 2
    )

    // cpus 8

    input:
      path input

    output:
      path "SWARM_representatives.fa.gz", emit: clust
      path "SWARM.uc.gz",                 emit: clust_uc
      path "SWARM.swarms.gz",             emit: swarms
      path "SWARM.struct.gz",             emit: struct
      path "SWARM.stats.gz",              emit: stats
      tuple val("${task.process}"), val('swarm'), eval('swarm --version 2>&1 | head -n 1 | sed "s/Swarm //"'), topic: versions

    exec:
      fastidious = (params.swarm_fastidious.toBoolean() == true & params.swarm_d.toInteger() == 1) ? "--fastidious --boundary ${params.swarm_d1boundary}" : ""
      // println("swarm_fastidious: ${params.swarm_fastidious}, swarm_d: ${params.swarm_d}")
      // println("fastid option:  ${fastidious}")

    script:
    """
    echo -e "Clustering sequences with SWARM\\n"
    echo -e "Note: sequences with ambiguous nucleotides will be excluded!\\n"

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

    echo -e "\\n..Swarm clustering finished\\n"

    ## Compress results
    echo -e "..Compressing results\\n"
    parallel -j 1 \
      "pigz -p ${task.cpus} -${params.gzip_compression} {}" \
      ::: "SWARM_representatives.fa" "SWARM.uc" "SWARM.swarms" "SWARM.struct" "SWARM.stats"

    echo -e "..Done\\n"
    """
}


// Pre-clustering / denoising / clustering subworkflow
workflow CLUSTERING {

    take:
    derep_ch

    main:

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Pre-clustering / denoising
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // No pre-clustering or denoizing
    if ( params.preclustering == "none" || params.preclustering == null ) {
      denoise_ch    = derep_ch
      preclustuc_ch = file('NoPrecluster')
      preclustaf_ch = file('NoPreclusterFASTA')
    
    // Denoise with UNOISE
    } else if ( params.preclustering == "unoise" ) {
      unoise(derep_ch)
      denoise_ch    = unoise.out.unoise
      preclustuc_ch = unoise.out.unoise_uc
      preclustaf_ch = denoise_ch
    
    // Denoise with DADA2
    } else if ( params.preclustering == "dada2" ) {

      // Denoise all dereplicated sequences
      if(params.dada2_pooling == "global"){
        dada2(derep_ch)
        denoise_ch    = dada2.out.dada
        preclustuc_ch = dada2.out.dada_uc
        preclustaf_ch = denoise_ch
      }

      // // Dereplicate and denoise by sequencing run
      // if(params.dada2_pooling == "byrun"){
      //   
      //   dereplication_byrun(ch_seqs)
      //   dada2(dereplication_byrun.out.dereps.flatten())
      //   
      //   dada2pool(
      //     dereplication.out.derep_uc,
      //     dada2.out.dada_ucr.collect()
      //     )
      // 
      //   /*
      //   denoise_ch    = dada2pool.out.dada
      //   preclustuc_ch = dada2pool.out.dada_uc
      //   preclustaf_ch = file('NoPreclusterFASTA')
      // 
      //   */
      // }


    // Precluster with SWARM
    } else if ( params.preclustering == "swarm_d1" ){
      precluster_swarm(derep_ch)
      denoise_ch    = precluster_swarm.out.clust
      preclustuc_ch = precluster_swarm.out.clust_uc
      preclustaf_ch = denoise_ch

    // Global homopolymer correction
    } else if ( params.preclustering == "homopolymer" ){
      homopolymer(derep_ch)
      denoise_ch    = homopolymer.out.hp
      preclustuc_ch = homopolymer.out.hp_uc
      preclustaf_ch = denoise_ch
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
        preclustaf_ch = file('NoPreclusterFASTA')
      
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
      preclustaf_ch = file('NoPreclusterFASTA')

    // Do not cluster, use ASVs from DADA2
    } else if ( params.preclustering == "dada2" & params.clustering == "none" ){
      
      if(params.dada2_pooling == "global"){
        cluster_ch = dada2.out.dada
        clustuc_ch = dada2.out.dada_uc
      } 
    
      preclustuc_ch = file('NoPrecluster')
      preclustaf_ch = file('NoPreclusterFASTA')
    
    } else if ( params.preclustering == "none" & params.clustering == "none" ){
      println "No pre-clustering or clustering was done"

      // TODO: create table based on dereplicated sequences?
    }


    emit:
    preclustuc_ch = preclustuc_ch    // UC file for pre-clustering
    preclustaf_ch = preclustaf_ch    // FASTA file for pre-clustering
    cluster_ch    = cluster_ch       // FASTA file for clustering
    clustuc_ch    = clustuc_ch       // UC file for clustering

} // end of subworkflow
