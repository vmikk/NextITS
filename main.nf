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

// Primer checks
params.primer_forward = "GTACACACCGCCCGTCG"    // ITS9mun
params.primer_reverse = "CCTSCSCTTANTDATATGC"  // ITS4ngsUni
params.primer_mismatches = 2
// params.primer_mismatches_insertions = 1
// params.primer_mismatches_deletions = 1
params.primer_foverlap = params.primer_forward.length() - 2
params.primer_roverlap = params.primer_reverse.length() - 2

// ITSx
params.ITSx_evalue = 1e-1
params.ITSx_partial = 0     // off, otherwise specify min length cutoff for partial ITS sequences to keep
params.ITSx_tax = "all"
/// params.ITSx_singledomain = true ....  optional arguments

// Homopolymer compression
params.hp_similarity = 0.999
params.hp_iddef = 2
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


// Primer disambiguation
process disambiguate {

    label "main_container"

    // publishDir "${out_2_primer}", mode: 'symlink'
    // cpus 1

    output:
      path "primer_F.fasta",  emit: F
      path "primer_R.fasta",  emit: R
      path "primer_Fr.fasta", emit: Fr
      path "primer_Rr.fasta", emit: Rr

    script:

    """

    ## Disambiguate forward primer
    echo -e "Disambiguating forward primer"
    disambiguate_primers.R \
      ${params.primer_forward} \
      primer_F.fasta

    ## Disambiguate reverse primer
    echo -e "\nDisambiguating reverse primer"
    disambiguate_primers.R \
      ${params.primer_reverse} \
      primer_R.fasta

    ## Reverse-complement primers
    echo -e "\nReverse-complementing primers"
    seqkit seq -r -p primer_F.fasta > primer_Fr.fasta
    seqkit seq -r -p primer_R.fasta > primer_Rr.fasta

    """
}


// Check primers + QC  //////////////////////////////////////////////////
// Count number of primer occurrences withnin a read,
// discard reads with > 1 primer occurrence
process primer_check {

    label "main_container"

    publishDir "${out_2_primer}", mode: 'symlink'

    // cpus 1

    // Add sample ID to the log file
    tag "${input.getSimpleName()}"

    input:
      path input
      path primer_F
      path primer_R
      path primer_Fr
      path primer_Rr

    output:
      path "${input.getSimpleName()}_PrimerChecked.fq.gz", emit: fq_primer_checked
      path "${input.getSimpleName()}_Mutiprimer.fq.gz", emit: mutiprimer, optional: true

    script:
    """
    echo -e "Input file: " ${input}

    ### Count number of pattern occurrences for each sequence
    count_primers (){
      # \$1 = file with primers

      seqkit locate \
        --max-mismatch ${params.primer_mismatches} \
        --only-positive-strand \
        --pattern-file "\$1" \
        --threads ${task.cpus} \
        ${input} \
        | awk -vOFS='\\t' 'NR > 1 { print \$1 , \$5 , \$6 }' \
        | bedtools merge -i stdin
    }

    echo -e "\nCounting primers"
    echo -e "..forward primer"
    count_primers ${primer_F}  >  PF.txt

    echo -e "..rc-forward primer"
    count_primers ${primer_Fr} >> PF.txt
    
    echo -e "..reverse primer"
    count_primers ${primer_R}  >  PR.txt

    echo -e "..rc-reverse primer"
    count_primers ${primer_Rr} >> PR.txt

    ## Sort by seqID and start position, remove overlapping regions,
    ## Find duplicated records
    echo -e "\nLooking for multiple primer occurrences"
    echo -e "..Processing forward primers"
    csvtk sort \
      -t -T -k 1:N -k 2:n \
      --num-cpus ${task.cpus} \
      PF.txt \
    | bedtools merge -i stdin \
    | awk '{ print \$1 }' \
    | runiq -i - \
    > multiprimer.txt

    echo -e "..Processing reverse primers"
    csvtk sort \
      -t -T -k 1:N -k 2:n \
      --num-cpus ${task.cpus} \
      PF.txt \
    | bedtools merge -i stdin \
    | awk '{ print \$1 }' \
    | runiq -i - \
    >> multiprimer.txt

    ## If some artifacts are found
    if [ -s multiprimer.txt ]; then
 
      ## Keep only uinque seqIDs
      runiq multiprimer.txt > multiprimers.txt
      rm multiprimer.txt
 
      echo -e "\nNumber of artefacts found: " \$(wc -l < multiprimers.txt)

      echo -e "..Removing artefacts"
      ## Remove multiprimer artefacts
      seqkit grep --invert-match --by-name \
        --threads ${task.cpus} \
        --pattern-file multiprimers.txt \
        --out-file no_multiprimers.fq.gz \
        ${input}

      ## Extract multiprimer artefacts
      echo -e "..Extracting artefacts"
      seqkit grep --by-name \
        --threads ${task.cpus} \
        --pattern-file multiprimers.txt \
        --out-file "${input.getSimpleName()}_Mutiprimer.fq.gz" \
        ${input}

      echo -e "..done"

    else

      echo -e "\nNo primer artefacts found"
      ln -s ${input} no_multiprimers.fq.gz
    
    fi
    echo -e "..Done"

    echo -e "\nReorienting sequences"

    cutadapt \
      -a ${params.primer_forward}";required;min_overlap=${params.primer_foverlap}"...${params.primer_reverse}";required;min_overlap=${params.primer_roverlap}" \
      --errors ${params.primer_mismatches} \
      --revcomp --rename "{header}" \
      --cores ${task.cpus} \
      --action none \
      --output ${input.getSimpleName()}_PrimerChecked.fq.gz \
      no_multiprimers.fq.gz \
      2> cutadapt.log

    echo -e "All done"

    ## Clean up
    if [ -f no_multiprimers.fq.gz ]; then rm no_multiprimers.fq.gz; fi

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




//  The default workflow
workflow {

    // Input file with multiplexed reads (FASTQ.gz)
    ch_input = Channel.value(params.input)

    // Input file with barcodes (FASTA)
    ch_barcodes = Channel.value(params.barcodes)

    // Demultiplexing
    demux(ch_input, ch_barcodes)

    // Primer disambiguation
    disambiguate()

    // Check primers
    primer_check(
      demux.out.samples_demux.flatten(),
      disambiguate.out.F,
      disambiguate.out.R,
      disambiguate.out.Fr,
      disambiguate.out.Rr
      )

    // Extract ITS
    itsx(primer_check.out.fq_primer_checked)

    // Homopolymer compression on full-length ITS sequences
    homopolymer(itsx.out.itsx_full)

}

