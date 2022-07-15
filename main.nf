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

}

