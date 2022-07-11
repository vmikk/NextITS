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


// Check primers
process primer_check {

    label "main_container"

    publishDir "${out_2_primer}", mode: 'symlink'
    // cpus 1

    // Add sample ID to the log file
    tag "${input.getSimpleName()}"

    input:
      path input

    output:
      path "${input.getSimpleName()}_PrimerChecked.fq.gz", emit: fq_primer_checked

    script:

    """

    echo -e "Reorinting sequences"
    echo -e "Input file: " ${input}
    echo -e "Forward primer: " ${params.primer_forward}
    echo -e "Reverse primer: " ${params.primer_reverse}
    echo -e "Number of mismathces allowed: " ${params.primer_mismatches}

    ## Convert IUPAC codes
    echo -e "\nConverting IUPAC codes"
    fwd_primer=\$(convert_IUPAC.sh ${params.primer_forward})
    rev_primer=\$(convert_IUPAC.sh ${params.primer_reverse})
    echo -e "IUPAC-expanded forward primer: " "\$fwd_primer"
    echo -e "IUPAC-expanded reverse primer: " "\$rev_primer"

    ## Forward (5' - 3' orinetation)
    ## + remove LIMA information from the header
    echo -e "\nSearching for forward primer"
    seqkit replace -p "\\s.+" ${input} \
    | fqgrep \
      -m ${params.primer_mismatches} \
      -i ${params.primer_mismatches_insertions} \
      -d ${params.primer_mismatches_deletions} \
      -e \
      -p "\$fwd_primer" \
      - \
      > 5_3.fastq
    echo -e "..Done"

    ## Reverse (3'-5' orientation - needs to be reverse-complemented)
    ## + remove LIMA information from the header
    echo -e "\nSearching for reverse primer"
    seqkit replace -p "\\s.+" ${input} \
    | fqgrep \
      -m ${params.primer_mismatches} \
      -i ${params.primer_mismatches_insertions} \
      -d ${params.primer_mismatches_deletions} \
      -e \
      -p "\$rev_primer" \
      - \
      > 3_5.fastq
    echo -e "..Done"

    ## If rev primer found, then make reverse complementary and merge with 5_3.fastq file
    echo -e "\nReverse-complementig 3'-5' sequences"
    if [ -s 3_5.fastq ]; then

      seqkit seq -t dna --validate-seq -r -p 3_5.fastq \
        >> 5_3.fastq

    fi
    echo -e "..Done"

    ## Searching for multi-primer artefacts
    echo -e "\nSearching for multi-primer artefacts"
    seqkit rmdup --by-name --threads ${task.cpus} \
      --dup-num-file duplicates.temp \
      --out-file 5_3.fastx.temp \
      5_3.fastq
    echo -e "..Done"


    # fqgrep -p "CGCCTGCGCTTAATTAT" -r Barcod01__A2201.fq.gz \
    #  | awk -F "\t" ' NR > 1 { print \$1 }' | sort | uniq -c | sort -r



    ## Remove multi-primer artefacts
    echo -e "\nRemoving multi-primer artefacts"
    ## If any duplicated sequences found
    if [ -s duplicates.temp ]; then
      
      ## Extract duplicate names
      awk 'BEGIN{FS=","}{print \$2}' duplicates.temp \
        | sed -e 's/^ //' > duplicates.names
      
      ## Remove duplicate seqs
      seqkit grep --invert-match --by-name \
        --pattern-file duplicates.names \
        --out-file out.fq \
        5_3.fastx.temp

      ## Get multi-primer artefacts
      # seqkit grep --by-name \
      #   --pattern-file duplicates.names \
      #   --out-file multiprimer_artefacts.fq \
      #   5_3.fastx.temp

      ## Count number of artefacts
      multiprimer_count=\$(wc -l duplicates.names | awk '{print \$1}')
      printf "Number of 'multi-primer' chimeric sequences found: \$multiprimer_count\n"

    ## No duplicated sequences found
    else

      mv 5_3.fastx.temp out.fq
      printf "No 'multi-primer' chimeric sequences found\n"

    fi
    echo -e "..Done"

    ## Check if reoriented output is empty; if yes, then report WARNING
    if [ -s out.fq ]; then
      :
    else
      printf '%s\n' "WARNING: primers not found (output is empty)"
    fi

    ## Compress results
    gzip -7 --stdout out.fq > "${input.getSimpleName()}_PrimerChecked.fq.gz"

    ## Remove temporary file
    rm out.fq
    if [ -f 5_3.fastq ];       then rm 5_3.fastq;       fi
    if [ -f 3_5.fastq ];       then rm 3_5.fastq;       fi
    if [ -f 5_3.fastx.temp ];  then rm 5_3.fastx.temp;  fi
    if [ -f duplicates.temp ]; then rm duplicates.temp; fi

    """
}


// Extract ITS region with ITSx
process itsx {

    label "main_container"

    publishDir "${out_3_itsx}", mode: 'symlink'
    // cpus 2

    // Add sample ID to the log file
    tag "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}"

    input:
      path input

    output:
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}_hash_table.txt.gz", emit: hashes
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}_uc.uc.gz", emit: uc
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.full.fasta.gz", emit: itsx_full, optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.SSU.fasta.gz",  emit: itsx_ssu, optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.ITS1.fasta.gz", emit: itsx_its1, optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.5_8S.fasta.gz", emit: itsx_58s, optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.ITS2.fasta.gz", emit: itsx_its2, optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.LSU.fasta.gz",  emit: itsx_lsu, optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.positions.txt", optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.problematic.txt", optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}_no_detections.fasta.gz", emit: itsx_nondetects, optional: true
      path "${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}.summary.txt", emit: itsx_summary

    script:
    
    sampID="${input.getSimpleName().replaceAll(/_PrimerChecked/, '')}"
    
    // Allow inclusion of sequences that only find a single domain, given that they meet the given E-value and score thresholds, on with parameters 1e-9,0 by default
    // singledomain = params.ITSx_singledomain ? "--allow_single_domain 1e-9,0" : ""

    """

    ## Sequence ID - Hash - Length - Average Phred score
    echo -e "Creating sequence hash table with average sequence quality"
    seqkit fx2tab --length --avg-qual ${input} \
      | hash_sequences.sh \
      | awk '{print \$1 "\t" \$6 "\t" \$4 "\t" \$5}' \
      > ${sampID}_hash_table.txt
    echo -e "..Done"

    ## Dereplicate at sample level
    echo -e "\nDereplicating at sample level"
    seqkit fq2fa -w 0 ${input} \
      | vsearch \
        --derep_fulllength - \
        --output - \
        --strand both \
        --fasta_width 0 \
        --threads 1 \
        --relabel_sha1 \
        --sizein --sizeout \
        --uc ${sampID}_uc.uc \
        --quiet \
      > derep.fasta
    echo -e "..Done"

    ## ITSx extraction
    echo -e "\nITSx extraction"
    ITSx \
      -i derep.fasta \
      --complement T \
      --save_regions all \
      --graphical F \
      --positions T \
      --not_found T \
      -E ${params.ITSx_evalue} \
      -t ${params.ITSx_tax} \
      --partial ${params.ITSx_partial} \
      --cpu ${task.cpus} \
      --preserve T \
      -o ${sampID}
    
    echo -e "..Done"

      # ITSx.full.fasta
      # ITSx.SSU.fasta
      # ITSx.ITS1.fasta
      # ITSx.5_8S.fasta
      # ITSx.ITS2.fasta
      # ITSx.LSU.fasta
      # ITSx.positions.txt
      # ITSx.problematic.txt
      # ITSx_no_detections.fasta
      # ITSx_no_detections.txt
      # ITSx.summary.txt

    ## Remove empty files (no sequences)
    echo -e "\nRemoving empty files"
    find . -type f -name "*.fasta" -empty -print -delete
    echo -e "..Done"

    ## Remove temporary file
    rm derep.fasta

    ## Compress results
    echo -e "\nCompressing files"
    gzip -7 ${sampID}_hash_table.txt
    gzip -7 ${sampID}_uc.uc
    gzip -7 *.fasta
    echo -e "..Done"

    """
}


// Homopolymer compression
process homopolymer {

    label "main_container"

    publishDir "${out_4_homop}", mode: 'symlink'
    // cpus 1

    // Add sample ID to the log file
    tag "${input.getSimpleName()}"

    input:
      path input

    output:
      path "${input.getSimpleName()}_Homopolymer_compressed.fa.gz", emit: hc, optional: true
      path "${input.getSimpleName()}_uch.uc.gz", emit: uch, optional: true

    script:
    sampID="${input.getSimpleName()}"

    """

    ## Homopolymer compression
    echo -e "Homopolymer compression"

    zcat ${input} \
      | homopolymer_compression.sh - \
      > homo_compressed.fa
    
    echo -e "..Done"

    ## Re-cluster homopolymer-compressed data
    echo -e "\nRe-clustering homopolymer-compressed data"
    vsearch \
      --cluster_size homo_compressed.fa \
      --id ${params.hp_similarity} \
      --iddef ${params.hp_iddef} \
      --qmask "dust" \
      --strand "both" \
      --fasta_width 0 \
      --threads ${task.cpus} \
      --sizein --sizeout \
      --centroids homo_clustered.fa \
      --uc ${sampID}_uch.uc
    echo -e "..Done"

    ## Compress UC file
    gzip -7 ${sampID}_uch.uc

    ## Substitute homopolymer-comressed sequences with uncompressed ones
    ## (update size annotaions)
    echo -e "\nExtracting sequences"

    seqkit fx2tab ${input} > inp_tab.txt
    seqkit fx2tab homo_clustered.fa > clust_tab.txt

    if [ -s inp_tab.txt ]; then
      substitute_compressed_seqs.R \
        inp_tab.txt clust_tab.txt res.fa

      echo -e "..Done"
    else
      echo -e "..Input data looks empty, nothin to proceed with"
    fi

    if [ -s res.fa ]; then
      gzip -c res.fa > ${sampID}_Homopolymer_compressed.fa.gz
    fi

    ## Remove temporary files
    rm homo_compressed.fa
    rm homo_clustered.fa
    rm inp_tab.txt
    rm clust_tab.txt
    rm res.fa

    """
}


