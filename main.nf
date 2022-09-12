#!/usr/bin/env nextflow
/*

============================================================================
  NextITS: Pipeline to process fungal ITS amplicons
============================================================================
  Version: v0.0.2
  License: Apache-2.0
  Github : https://github.com/vmikk/NextITS
  Website: TBA
----------------------------------------------------------------------------
*/


// Dependencies:
//  - lima >= 2.6.0
//  - fqgrep >= 0.4.4
//  - vsearch >= 2.21.1
//  - seqkit >= 2.2.0
//  - cutadapt >= 4.1
//  - fastp >= 0.23.2
//  - ITSx >= 1.1.3
//  - BLAST 2.12.0+
//  - R >= 4.1.0
//    -- data.table >= 1.14.0
//    -- Biostrings >= 2.60.0
//    -- plyr
//    -- DECIPHER >= 2.24.0
//    -- R.utils
//    -- ggplot2
//    -- openxlsx
//    -- optparse
//  - bioawk >= 20110810
//  - miller (mlr) >= 6.2.0
//  - bedtools >= 2.30.0
//  - GNU parallel
//  - rush >= 0.4.2
//  - csvtk >= 0.23.0
//  - runiq >= 1.2.1
//  - awk, sed, gzip, find, rename, cat, zcat, sha1sum
//
// Databases:
//  - UDB for chimera identification
//  - BlastDB for taxonomy annotation


// Enable DSL2 syntax
nextflow.enable.dsl = 2

// Pipeline version
version = '0.0.2'

// Initialize parameters, set default values
params.data_path = "${projectDir}/pipeline_data"

params.input = false
params.input_R1 = false
params.input_R2 = false
params.outdir = "${launchDir}/results"

// Help message flag
params.helpMsg = false

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

// Reference-based chimera removal
params.chimera_db = "/mnt/Dat2/DB/UNITE/Leho_Subset/UN95_chimera.udb"
params.chimera_rescueoccurrence = 2

// De novo chimera identification (UCHIME1)
params.chimeranov_abskew = 2.0
params.chimeranov_dn = 1.4
params.chimeranov_mindiffs = 3
params.chimeranov_mindiv = 0.8
params.chimeranov_minh = 0.28
params.chimeranov_xn = 8.0

// OTU clustering (for tag-jump removal)
params.otu_id = 0.98
params.otu_iddef = 2

// Tag-jump removal
params.tj_f = 0.01    // UNCROSS parameter f
params.tj_p = 1

// Taxonomy annotation
params.blast_taxdb = false
params.blast_task = "blastn"   // or "megablast" 
params.blast_chunksize = 100
params.blast_maxts = 10
params.blast_hsps = 1
// params.blast_wordsize


if(params.blast_taxdb){
  bastdb_name = file(params.blast_taxdb).name
  bastdb_dir = file(params.blast_taxdb).parent
}


// Pipeline help message
def helpMsg() {
    log.info"""
    =====================================================================
    NextITS ${version}
    =====================================================================
    
    Pipeline Usage:
    To run the pipeline, enter the following in the command line:
        nextflow run vmikk/nextits -r ${version} --input ... --outdir ...
    
    Options:
    REQUIRED:
        --input               Path to the directory with parquet files (GBIF occurrcence dump)
        --outdir              The output directory where the results will be saved

    OPTIONAL:
        --phylum              Phylum to analyze (multiple comma-separated values allowed); e.g., "Chordata"
        --class               Class to analyze (multiple comma-separated values allowed); e.g., "Mammalia"
    NEXTFLOW-SPECIFIC:
        -qs                   Queue size (max number of processes that can be executed in parallel); e.g., 8
    """.stripIndent()
}
// Show help msg
if (params.helpMsg){
    helpMsg()
    exit(0)
}

// Check if input path was provided
if (params.input == false) {
    println( "Please provide the input file with sequences in FASTQ.gz format with `--input` parameter.")
    exit(1)
}
if (params.barcodes == false) {
    println( "Please provide the file with sample barcodes in FASTA format with `--barcodes` parameter.")
    exit(1)
}
if (params.chimera_db == false) {
    println( "Please provide the UDB file with reference sequences for chimera removal with `--chimera_db` parameter.")
    exit(1)
}


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

log.info """
        Pipeline info:
          Pipeline profile:       ${workflow.profile}
          Config file used:       ${workflow.configFiles}
          Container engine:       ${workflow.containerEngine}
        """
        .stripIndent()

log.info """
        Core Nextflow options:
          launchDir:              ${workflow.launchDir}
          workDir:                ${workflow.workDir}
          projectDir:             ${workflow.projectDir}
        """
        .stripIndent()

log.info "======================================================================="
log.info "\n"



// Define output paths for different steps
out_0        = params.outdir
out_1_demux  = params.outdir + "/01_Demux"
out_2_primer = params.outdir + "/02_PrimerCheck"
out_3_itsx   = params.outdir + "/03_ITSx"
out_4_homop  = params.outdir + "/04_Homopolymer"
out_5_chim   = params.outdir + "/05_Chimera"
out_6_tj     = params.outdir + "/06_TagJumpFiltration"
out_7_asv    = params.outdir + "/07_ASV_table"
out_8_blast  = params.outdir + "/08_Taxonomy"




// Quality filtering for single-end reads
process qc_se {

    label "main_container"

    // cpus 10

    input:
      path input

    output:
      path "QC.fq.gz", emit: filtered

    script:
    filter_maxee      = params.qc_maxee      ? "--fastq_maxee ${params.qc_maxee}"          : ""
    filter_maxeerate  = params.qc_maxeerate  ? "--fastq_maxee_rate ${params.qc_maxeerate}" : ""
    """
    echo -e "QC"
    echo -e "Input file: " ${input}

    vsearch \
      --fastq_filter ${input} \
      --fastq_qmax 93 \
      ${filter_maxee} \
      ${filter_maxeerate} \
      --fastq_maxns ${params.qc_maxn} \
      --threads ${task.cpus} \
      --fastqout - \
    | gzip -7 \
    > QC.fq.gz

    echo -e "\nQC finished"
    """
}


// Quality filtering for pair-end reads
process qc_pe {

    label "main_container"

    // cpus 10

    input:
      path input_R1
      path input_R2

    output:
      path "QC_R1.fq.gz", emit: filtered_R1
      path "QC_R2.fq.gz", emit: filtered_R2

    script:
    filter_avgphred  = params.qc_avgphred  ? "--average_qual ${params.qc_avgphred}"                 : "--average_qual 0"
    filter_phredmin  = params.qc_phredmin  ? "--qualified_quality_phred ${params.qc_phredmin}"      : ""
    filter_phredperc = params.qc_phredperc ? "--unqualified_percent_limit ${params.qc_phredperc}"   : ""
    filter_polyglen  = params.qc_polyglen  ? "--trim_poly_g --poly_g_min_len ${params.qc_polyglen}" : ""
    """
    echo -e "QC"
    echo -e "Input R1: " ${input_R1}
    echo -e "Input R2: " ${input_R2}

    ## If `filter_phredmin` && `filter_phredperc` are specified,
    # Filtering based on percentage of unqualified bases
    # how many percents of bases are allowed to be unqualified (Q < 24)
        
    fastp \
      --in1 ${input_R1} \
      --in2 ${input_R2} \
      --disable_adapter_trimming \
      --n_base_limit ${params.qc_maxn} \
      ${filter_avgphred} \
      ${filter_phredmin} \
      ${filter_phredperc} \
      ${filter_polyglen} \
      --length_required 100 \
      --thread ${task.cpus} \
      --html qc.html \
      --json qc.json \
      --out1 QC_R1.fq.gz \
      --out2 QC_R2.fq.gz

    echo -e "\nQC finished"
    """
}




// Demultiplexing with LIMA - for PacBio reads
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


// Merge Illumina PE reads
process merge_pe {

    label "main_container"

    // publishDir "${out_1_demux}", mode: 'symlink'
    // cpus 10

    input:
      path input_R1
      path input_R2

    output:
      path "Merged.fq.gz", emit: r12
      tuple path("NotMerged_R1.fq.gz"), path("NotMerged_R2.fq.gz"), emit: nm, optional: true

    script:
    """
    echo -e "Merging Illumina pair-end reads"

    ## By default, fastp modifies sequences header
    ## e.g., `merged_150_15` means that 150bp are from read1, and 15bp are from read2
    ## But we'll preserve only sequence ID

    fastp \
      --in1 ${input_R1} \
      --in2 ${input_R2} \
      --merge --correction \
      --overlap_len_require ${params.pe_minoverlap} \
      --overlap_diff_limit ${params.pe_difflimit} \
      --overlap_diff_percent_limit ${params.pe_diffperclimit} \
      --length_required ${params.pe_minlen} \
      --disable_quality_filtering \
      --disable_adapter_trimming \
      --dont_eval_duplication \
      --compression 6 \
      --thread ${task.cpus} \
      --out1 NotMerged_R1.fq.gz \
      --out2 NotMerged_R2.fq.gz \
      --json log.json \
      --html log.html \
      --stdout \
    | seqkit seq --only-id \
    | gzip -7 \
    > Merged.fq.gz

    ##  --merged_out Merged.fq.gz \
    ##  --n_base_limit ${params.pe_nlimit} \

    # --overlap_len_require         the minimum length to detect overlapped region of PE reads
    # --overlap_diff_limit          the maximum number of mismatched bases to detect overlapped region of PE reads
    # --overlap_diff_percent_limit  the maximum percentage of mismatched bases to detect overlapped region of PE reads
    ## NB: reads should meet these three conditions simultaneously!

    echo -e "..done"
    """
}


// Modify barcodes for cutadapt (restrict the search window)
process prep_barcodes {

    label "main_container"

    // publishDir "${out_1_demux}", mode: 'symlink'
    // cpus 1

    input:
      path barcodes

    output:
      path "barcodes_modified.fa", emit: barcodesm

    script:
    """
    echo -e "Restricting the search window for barcode lookup"
    echo -e "Provided barcodes: " ${barcodes}

    ## Add `XN{30}` to the barcodes

    sed -e '/^>/! s/^/XN{${params.barcode_window}}/' \
      ${barcodes} \
      > barcodes_modified.fa

    echo -e "..Done"
    """
}


// Demultiplexing with cutadapt - for Illumina SE reads
process demux_illumina {

    label "main_container"

    publishDir "${out_1_demux}", mode: 'symlink'
    // cpus 10

    input:
      path input_fastq
      path barcodes

    output:
      path "*.fq.gz", emit: samples_demux

    script:
    """
    echo -e "Input file: " ${input_fastq}
    echo -e "Barcodes: " ${barcodes}

    echo -e "\nDemultiplexing with cutadapt:"

    ## Demultiplex with cutadapt
    cutadapt \
      -g file:${barcodes} \
      --revcomp --rename "{header}" \
      --errors ${params.barcode_errors} \
      --overlap ${params.barcode_overlap} \
      --no-indels \
      --cores ${task.cpus} \
      --discard-untrimmed \
      --action none \
      -o "{name}.fq.gz" \
      ${input_fastq} \
      > cutadapt.log

    echo -e "\n..done"

    ## Remove empty files (no sequences)
    echo -e "\nRemoving empty files"
    find . -type f -name "*.fq.gz" -size -29c -print -delete
    echo -e "..Done"

    echo -e "\nDemultiplexing finished"
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


// Merge tables with sequence qualities
process seq_qual {

    label "main_container"

    // cpus 1

    input:
      path input

    output:
      path "SeqQualities.txt.gz", emit: quals

    script:
    """
    echo -e "Aggregating sequence qualities"

    find . -maxdepth 1 -name "*_hash_table.txt.gz" \
      | parallel -j1 "merge_sequnce_qualities.sh {} {/.}" \
      | gzip -7 \
      > SeqQualities.txt.gz

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



// Reference-based chimera removal
process chimera_ref {

    label "main_container"

    publishDir "${out_5_chim}", mode: 'symlink'
    // cpus 1

    // Add sample ID to the log file
    tag "${input.getSimpleName().replaceAll(/_Homopolymer_compressed/, '')}"

    input:
      path input
      path DB

    output:
      path "${input.getSimpleName().replaceAll(/_Homopolymer_compressed/, '')}_NoChimera.fa.gz", emit: nonchimeric, optional: true
      path "${input.getSimpleName().replaceAll(/_Homopolymer_compressed/, '')}_Chimera.fa.gz", emit: chimeric, optional: true

    script:
    sampID="${input.getSimpleName().replaceAll(/_Homopolymer_compressed/, '')}"

    """

    ## Sample name will be added to the header of chimeric sequences
    # sampID="\$(basename ${input} _Homopolymer_compressed.fa.gz)"

    ## Reference database based chimera filtering
    echo -e "Reference-based chimera removal"
    vsearch \
      --uchime_ref ${input} \
      --db ${DB} \
      --selfid \
      --fasta_width 0 \
      --threads ${task.cpus} \
      --sizein --sizeout \
      --chimeras chimeras.fasta \
      --nonchimeras nonchimeras.fasta \
      --borderline borderline.fasta

    # --selfid = ignore reference sequences that are 100% identical to the query
    echo -e "..Done"


    ## Add borderline sequences to non-chimeric sequences
    if [ -e borderline.fasta ]
    then
        echo -e "\nBorderline sequences were added to non-chimeric sequences"
        cat borderline.fasta nonchimeras.fasta > nc_bo.fasta
        mv nc_bo.fasta nonchimeras.fasta
        rm borderline.fasta
    fi

    ## Chimeric sequences
    if [ -e chimeras.fasta ]
    then
      ## Add sample ID to the header and compress the file
      sed 's/>.*/&;sample='"${sampID}"';/' chimeras.fasta \
        | gzip -7 \
        > "${sampID}_Chimera.fa.gz"
      rm chimeras.fasta
    else
      echo -e "\nNo chimeras detected"
      touch "${sampID}_Chimera.fa.gz"
    fi

    ## Non-chimeric sequences
    if [ -e nonchimeras.fasta ]
    then
      gzip -c nonchimeras.fasta > "${sampID}_NoChimera.fa.gz"
      rm nonchimeras.fasta
    else
      echo "No non-chimeric sequences left"
      touch "${sampID}_NoChimera.fa.gz"
    fi

    """
}


// Recovery of ref-based chimeric sequences with high occurrence
process chimera_rescue {

    label "main_container"

    publishDir "${out_5_chim}", mode: 'symlink'
    // cpus 1

    input:
      path input

    output:
      path "*_RescuedChimera.fa.gz", emit: rescuedchimeric, optional: true

    script:
    """

    ## Aggregate chimeric sequences from different samples
    echo -e "\nAggregating chimeric sequences"
    find . -name "*_Chimera.fa.gz" \
      | parallel -j1 "zcat {}" \
      | seqkit fx2tab \
      | sed -r 's:\t+:\t:g' | sed 's/\t\$//g' \
      | gzip > All_chimeras.txt.gz
    echo -e "..Done"

    ### Inspect chimerae occurrence
    ## Rescue sequences that were annotated as chimeric, 
    ## but have high occurrence within sequenceing run (e.g., occurrence > 2)
    echo -e "\nInspecting occurrence of chimeric sequences"

    chimera_rescue.R \
      "All_chimeras.txt.gz" \
      ${params.chimera_rescueoccurrence} \
      "Rescued_Chimeric_sequences.fa.gz"

    echo -e "..Done"

    ## Split rescured sequences by sample
    if [ -e Rescued_Chimeric_sequences.fa.gz ]
    then
      echo -e "\n..Splitting rescued sequences by sample"
      seqkit split -i \
        --id-regexp ";sample=(.*);" \
        --threads ${task.cpus} \
        -w 0 \
        -O Rescued_by_sample \
        Rescued_Chimeric_sequences.fa.gz

      rename \
        --filename 's/^Rescued_Chimeric_sequences.id_//g; s/.fa.gz/_RescuedChimera.fa.gz/' \
        \$(find Rescued_by_sample -name "*.fa.gz")

      mv *_RescuedChimera.fa.gz .

      echo -e "..Done"
    fi

    ## Remove temporary files
    rm All_chimeras.txt.gz
    if [ -f Rescued_Chimeric_sequences.fa.gz ]; then rm Rescued_Chimeric_sequences.fa.gz; fi

    """
}



// De novo chimera identification
// NB! in uchime_denovo, sequences are compared on their plus strand only!
process chimera_denovo {

    label "main_container"

    publishDir "${out_5_chim}", mode: 'symlink'
    // cpus 1

    // Add sample ID to the log file
    tag "${input.getSimpleName().replaceAll(/_Homopolymer_compressed/, '')}"

    input:
      path input

    output:
      path "${input.getSimpleName().replaceAll(/_Homopolymer_compressed/, '')}_DeNovoChim.txt", emit: denovochim, optional: true

    script:
    sampID="${input.getSimpleName().replaceAll(/_Homopolymer_compressed/, '')}"

    """

    ## Input order matters for chimera detection,
    ## so sequences will be automatically sorted by decreasing abundance first

    echo -e "De novo chimera identification"

    vsearch \
      --uchime_denovo ${input} \
      --abskew ${params.chimeranov_abskew} \
      --dn ${params.chimeranov_dn} \
      --mindiffs ${params.chimeranov_mindiffs} \
      --mindiv ${params.chimeranov_mindiv} \
      --minh ${params.chimeranov_minh} \
      --xn ${params.chimeranov_xn} \
      --threads 1 \
      --qmask dust \
      --sizein --xsize \
      --fasta_width 0 \
      --fasta_score \
      --chimeras - \
    | seqkit seq --name \
    | sed 's/;+/;/g ; s/;/\t/g ; s/uchime_denovo=//' \
    | sed 's/\$/\t${sampID}/' \
    > ${sampID}_DeNovoChim.txt
      
      # --uchimeout uchimeout.txt

    ## Remove file, if empty
    find . -maxdepth 1 -name ${sampID}_DeNovoChim.txt -size 0 -print -delete

    echo -e "..Done"

    """
}

// Aggregate de novo chimeras into a single file
process chimera_denovo_agg {

    label "main_container"
    // cpus 1

    input:
      path input

    output:
      path "DeNovo_Chimera.txt", emit: alldenovochim

    script:
    """
    echo -e "Aggregating de novo chimeric sequences"

    find . -name "*_DeNovoChim.txt" \
      | parallel -j1 "cat {}" \
      > DeNovo_Chimera.txt

    echo -e "..Done"

    """
}



// Merge and dereplicate sequences from all samples
process glob_derep {

    label "main_container"

    // cpus 1

    input:
      path input

    output:
      path "Derep_for_clust.fa.gz", emit: globderep
      path "Derep_for_clust.uc.gz", emit: globderep_uc

    script:
    """

    ## Concatenate non-chimeric sequences
    # zcat *_NoChimera.fa.gz *_RescuedChimera.fa.gz

    echo -e "\n ===== Files with non-chimeric sequences:"
    if test -n "\$(find . -maxdepth 1 -name '*_NoChimera.fa.gz' -print -quit)"
    then
        ls -l *_NoChimera.fa.gz;
    else
        echo "..No files"
    fi

    echo -e "\n ===== Files with rescued chimeric sequences:"
    if test -n "\$(find . -maxdepth 1 -name '*_RescuedChimera.fa.gz' -print -quit)"
    then
        ls -l *_RescuedChimera.fa.gz;
    else
        echo "..No files"
    fi

    echo -e "\nDereplicating sequences"
    find . -name "*Chimera.fa.gz" | parallel -j1 \
      "zcat {}" \
      | sed '/^>/ s/;sample=.*;/;/' \
      | vsearch \
        --derep_fulllength - \
        --output - \
        --strand both \
        --fasta_width 0 \
        --threads 1 \
        --sizein --sizeout \
        --uc Derep_for_clust.uc \
      | gzip -7 > Derep_for_clust.fa.gz
    echo -e "..Done"

    ## Compress UC file
    echo -e "\nCompressing UC file"
    gzip -7 Derep_for_clust.uc

    """
}


// De-novo clustering of sequences into OTUs (for tag-jump removal)
process otu_clust {

    label "main_container"

    publishDir "${out_6_tj}", mode: 'symlink'
    // cpus 10

    input:
      path input

    output:
      path "OTUs.fa.gz", emit: otus
      path "OTUs.uc.gz", emit: otus_uc

    script:
    """

    echo -e "\nClustering sequences"
    vsearch \
      --cluster_size ${input} \
      --id ${params.otu_id} \
      --iddef ${params.otu_iddef} \
      --sizein --sizeout \
      --qmask dust --strand both \
      --fasta_width 0 \
      --uc OTUs.uc \
      --threads ${task.cpus} \
      --centroids - \
    | gzip -7 > OTUs.fa.gz
    
    echo -e "..Done"

    ## Compress UC file
    echo -e "\nCompressing UC file"
    gzip -7 OTUs.uc

    """
}


// Pool sequences and add sample ID into header (for OTU and ASV table creation)
process pool_seqs {

    label "main_container"
    
    publishDir "${out_6_tj}", mode: 'symlink'
    // cpus 3

    input:
      path input

    output:
      path "ASV_tab_not_filtered.txt.gz", emit: asvtabnf
      path "ASV_not_filtered.fa.gz", emit: asvsnf

    script:
    """

    echo -e "\nPooling and renaming sequences"

    parallel -j ${task.cpus} --group \
      "zcat {} \
        | sed 's/>.*/&;sample='{/.}';/ ; s/_NoChimera.fa//g ; s/_RescuedChimera//g '" \
      ::: *.fa.gz \
    | gzip -7 > ASV_not_filtered.fa.gz

    echo "..Done"


    echo -e "\nExtracting sequence count table"
    seqkit seq --name ASV_not_filtered.fa.gz \
      | sed 's/;/\t/g; s/size=//; s/sample=// ; s/\t*\$//' \
      | gzip -7 > ASV_tab_not_filtered.txt.gz

    echo "..Done"

    """
}


// Create OTU table (for tag-jump removal)
process otu_tab {

    label "main_container"

    publishDir "${out_6_tj}", mode: 'symlink'
    // cpus 10

    input:
      path otus
      path sample_seqs

    output:
      path "OTU_tab_not_filtered.txt.gz", emit: otutab
      path "Sample_mapping.uc.gz", emit: samples_uc

    script:
    """
    
    ## Read mapping (for OTU table)
    echo -e "\nRead mapping"
    vsearch \
      --usearch_global "${sample_seqs}" \
      --db "${otus}" \
      --id ${params.otu_id} \
      --strand both \
      --qmask none \
      --dbmask none \
      --sizein --sizeout \
      --otutabout "OTU_tab_not_filtered.txt" \
      --uc Sample_mapping.uc \
      --threads ${task.cpus}

    echo -e "..Done"

    ## Compress UC file
    echo -e "\nCompressing results"
    gzip -7 OTU_tab_not_filtered.txt
    gzip -7 Sample_mapping.uc

    """
}


// Tag-jump removal
process tj {

    label "main_container"

    publishDir "${out_6_tj}", mode: 'symlink'
    // cpus 1

    input:
      path otutab

    output:
      path "OTU_tab_TagJumpFiltered.txt.gz", emit: otutabtj
      path "TagJump_OTUs.RData", emit: tjs
      path "TagJump_plot.pdf"

    script:
    """

    echo -e "Tag-jump removal"
    
    tag_jump_removal.R \
      ${otutab} \
      ${params.tj_f} \
      ${params.tj_p}

    echo "..Done"

    """
}


// Prepare ASV table only with non-tag-jumped sequences
// Add quality estimate to singletons
// Add chimera-scores for putative de novo chimeras
process prep_asvtab {

    label "main_container"

    publishDir "${out_7_asv}", mode: 'symlink'
    // cpus 1

    input:
      path asvtabnf
      path asvsnf
      path mappings
      path tagjumps
      path denovos
      path quals

    output:
      path "ASVs.txt.gz", emit: asv_tl
      path "ASV_tab.txt.gz", emit: asv_tw
      path "ASVs.fa.gz", emit: asv_fa
      path "ASVs.RData", emit: asv_rd


    script:
    """

    echo -e "ASV table creation"
    
    asv_table_assembly.R \
      ${asvtabnf} \
      ${asvsnf}   \
      ${mappings} \
      ${tagjumps} \
      ${denovos}  \
      ${quals}

    echo "..Done"

    """
}




// // Count number of reads
// process read_counts {
// 
//     label "main_container"
// 
//     publishDir "${out_0}", mode: 'symlink'
//     // cpus 5
// 
//     input:
//       path input_fastq
//       path samples_demux
//       path samples_primerch
// 
//     output:
//       path "Counts_1.RawData.txt",     emit: counts_1
//       path "Counts_2.Demux.txt",       emit: counts_2
//       path "Counts_3.PrimerCheck.txt", emit: counts_3
// 
//     script:
// 
//     """
// 
//     ## Count raw reads
//     seqkit stat --basename --tabular --threads ${task.cpus} \
//       ${input_fastq} > Counts_1.RawData.txt
// 
//     ## Count demultiplexed reads
//     seqkit stat --basename --tabular --threads ${task.cpus} \
//       ${samples_demux} > Counts_2.Demux.txt
// 
//     ## Count primer-checked reads
//     seqkit stat --basename --tabular --threads ${task.cpus} \
//       ${samples_demux} > Counts_3.PrimerCheck.txt
// 
//     """
// }



// Taxonomy annotation
// use globally-dereplicated sequences
process blastn {

    label "main_container"

    // publishDir "${out_8_blast}", mode: 'symlink'
    // cpus 1

    input:
      path input
      path taxdb_dir

    output:
      path "${input.getBaseName()}.m8.gz", emit: blast

    script:
    """

    echo -e "Taxonomy annotation with BLASTn"
    echo -e "Input file: " ${input}

    blastn \
      -task ${params.blast_task} \
      -outfmt=6 -strand both \
      -query ${input} \
      -db ${taxdb_dir}/${bastdb_name} \
      -max_target_seqs ${params.blast_maxts} \
      -max_hsps ${params.blast_hsps} \
      -out ${input.getBaseName()}.m8 \
      -num_threads ${task.cpus}

      # -word_size
      # -evalue
      # -perc_identity

    ## Compress results
    gzip -7 ${input.getBaseName()}.m8

    echo "..Done"

    """
}


// Aggregate BLAST results
process blast_merge {

    label "main_container"

    publishDir "${out_8_blast}", mode: 'symlink'
    // cpus 1

    input:
      path input

    output:
      path "Blast_hits.m8.gz", emit: blast

    script:
    """

    echo -e "Aggregating BLASTn hits"
    
    cat *.m8.gz > Blast_hits.m8.gz

    echo "..Done"

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

    // Merge tables with sequence qualities
    seq_qual(itsx.out.hashes.collect())

    // Homopolymer compression on full-length ITS sequences
    homopolymer(itsx.out.itsx_full)

    // Reference-based chimera removal
    ch_chimerabd = Channel.value(params.chimera_db)
    chimera_ref(homopolymer.out.hc, ch_chimerabd)

    // Chimera rescue
    ch_chimerafiles = chimera_ref.out.chimeric.collect()
    chimera_rescue(ch_chimerafiles)

    // De novo chimera search
    chimera_denovo(homopolymer.out.hc)

    // Aggregate de novo chimeras into a single file
    chimera_denovo_agg(chimera_denovo.out.denovochim.collect())

    // Create channel with filtered reads
    ch_filteredseqs = chimera_ref.out.nonchimeric
      .concat(chimera_rescue.out.rescuedchimeric)
      .collect()

    // Global dereplication
    glob_derep(ch_filteredseqs)

    // Pool sequences (for ASV table)
    pool_seqs(ch_filteredseqs)

    // OTU clustering
    otu_clust(glob_derep.out.globderep)

    // Create OTU table
    otu_tab(
      otu_clust.out.otus,
      pool_seqs.out.asvsnf)

    // Tag-jump removal
    tj(otu_tab.out.otutab)

    // Create ASV table
    prep_asvtab(
      pool_seqs.out.asvtabnf,               // non-filtered ASV table
      pool_seqs.out.asvsnf,                 // ASV sequences in FASTA
      otu_tab.out.samples_uc,               // sequence mapping to OTUs
      tj.out.tjs,                           // tag-jumped OTU list
      chimera_denovo_agg.out.alldenovochim, // de novo chimera scores
      seq_qual.out.quals                    // sequence qualities
      )



    // BLAST (optional)
    if ( params.blast_taxdb ) {

      // Split FASTA sequences (globally dereplicated) into 
      ch_fasta = glob_derep.out.globderep
          .splitFasta(by: params.blast_chunksize, file:true)

      // Taxonomy annotation
      blastn(ch_fasta, bastdb_dir)

      // Aggregate BLAST results
      ch_blasthits = blastn.out.blast.collect()
      blast_merge(ch_blasthits)

      // Parse BLAST results
      // parse_blast()
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
