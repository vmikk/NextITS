


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
    echo -e "QC\\n"
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

    echo -e "\\nQC finished"
    """
}



// Merge Illumina PE reads
process merge_pe {

    label "main_container"

    // publishDir "${out_1_demux}", mode: "${params.storagemode}"
    // cpus 10

    input:
      path input_R1
      path input_R2

    output:
      path "Merged.fq.gz", emit: r12
      tuple path("NotMerged_R1.fq.gz"), path("NotMerged_R2.fq.gz"), emit: nm, optional: true

    script:
    """
    echo -e "Merging Illumina pair-end reads\\n"

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
    | gzip -${params.gzip_compression} \
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





// Demultiplexing with cutadapt - for Illumina PE reads (only not merged)
// NB. it's possible to use anchored adapters (e.g., -g ^file:barcodes.fa),
//     but there could be a preceding nucleotides before the barcode,
//     therefore, modified barcodes would be used here (e.g., XN{30})
process demux_illumina_notmerged {

    label "main_container"

    publishDir "${out_1_demux}", mode: 'symlink'
    // cpus 20

    input:
      path input_R1
      path input_R2
      path barcodes   // barcodes_modified.fa (e.g., XN{30})

    output:
      path "*.fq.gz", emit: samples_demux

    script:
    """
    echo -e "\nDemultiplexing not-merged reads"

    echo -e "Input R1: " ${input_R1}
    echo -e "Input R2: " ${input_R2}
    echo -e "Barcodes: " ${barcodes}

    ## First round
    echo -e "\nRound 1:"

    cutadapt -g file:${barcodes} \
      -o round1-{name}.R1.fastq.gz \
      -p round1-{name}.R2.fastq.gz \
      --errors ${params.barcode_errors} \
      --overlap ${params.barcode_overlap} \
      --no-indels \
      --cores ${task.cpus} \
      ${input_R1} ${input_R2} \
      > cutadapt_round_1.log

    echo -e ".. round 1 finished"
    
    ## Second round
    echo -e "\nRound 2:"

    cutadapt -g file:${barcodes} \
      -o round2-{name}.R2.fastq.gz \
      -p round2-{name}.R1.fastq.gz \
      --errors ${params.barcode_errors} \
      --overlap ${params.barcode_overlap} \
      --no-indels \
      --cores ${task.cpus} \
      round1-unknown.R2.fastq.gz round1-unknown.R1.fastq.gz \
      > cutadapt_round_2.log

    echo -e ".. round 2 finished"

    ## Remove empty files (no sequences)
    echo -e "\nRemoving empty files"
    find . -type f -name "round*.fastq.gz" -size -29c -print -delete
    echo -e "..Done"

    ## Remove unknowns
    echo -e "Removing unknowns"
    rm round1-unknown.R{1,2}.fastq.gz
    rm round2-unknown.R{1,2}.fastq.gz

    ## Combine sequences from round 1 and round 2 for each sample
    echo -e "\nCombining sequences from round 1 and round 2 for each sample"
    
    mkdir -p Combined
    
    find . -name "round*.R1.fastq.gz" | sort | parallel -j1 \
      "cat {} >> Combined/{= s/round1-//; s/round2-// =}"

    find . -name "round*.R2.fastq.gz" | sort | parallel -j1 \
      "cat {} >> Combined/{= s/round1-//; s/round2-// =}"

    echo -e "..Done"

    ## Clean up
    echo -e "..Removing temporary files"
    find . -type f -name "round*.fastq.gz" -print -delete



    echo -e "\nDemultiplexing finished"
    """
}


// Trim primers of nonmerged PE reads
// + Estimate sequence qualities
process trim_primers_pe {

    label "main_container"

    publishDir "${out_3_trimPE}", mode: 'symlink'
    // cpus 2

    // Add sample ID to the log file
    tag "${input.getSimpleName()}"

    input:
      path input   // tuple of size 2

    output:
      path "${input.getSimpleName()}_R1.fa.gz", emit: primertrimmed_fa_R1, optional: true
      path "${input.getSimpleName()}_R2.fa.gz", emit: primertrimmed_fa_R2, optional: true
      path "${input.getSimpleName()}_hash_table_R1.txt.gz", emit: hashes_R1, optional: true
      path "${input.getSimpleName()}_hash_table_R2.txt.gz", emit: hashes_R2, optional: true
      path "${input.getSimpleName()}_R1.fq.gz", emit: primertrimmed_fq_R1, optional: true
      path "${input.getSimpleName()}_R2.fq.gz", emit: primertrimmed_fq_R2, optional: true
      path "${input.getSimpleName()}_uc_R1.uc.gz", emit: ucR1, optional: true
      path "${input.getSimpleName()}_uc_R2.uc.gz", emit: ucR2, optional: true

    script:
    sampID="${input.getSimpleName()}"
    
    """
    echo -e "Input sample: " ${sampID}
    echo -e "Forward primer: " ${params.primer_forward}
    echo -e "Reverse primer: " ${params.primer_reverse}

    ## Reverse-complement primers
    FR=\$(rc.sh ${params.primer_forward})
    RR=\$(rc.sh ${params.primer_reverse})
    
    echo -e "Forward primer RC: " "\$FR"
    echo -e "Reverse primer RC: " "\$RR"

    ## Discard sequences without both primers
    echo -e "\nChecking primers"
    
    echo -e "..Forward strain"

    cutadapt \
      -a ${params.primer_forward}";required;min_overlap=${params.primer_foverlap}"..."\$RR"";required;min_overlap=${params.primer_roverlap}" \
      --errors ${params.primer_mismatches} \
      --cores ${task.cpus} \
      --action=none \
      --discard-untrimmed \
      -o for_R1.fastq.gz -p for_R2.fastq.gz \
      ${input[0]} ${input[1]} \
      > cutadapt_1.log

    
    echo -e "..Reverse strain"

    cutadapt \
      -a ${params.primer_reverse}";required;min_overlap=${params.primer_roverlap}"..."\$FR"";required;min_overlap=${params.primer_foverlap}" \
      --errors ${params.primer_mismatches} \
      --cores ${task.cpus} \
      --action=none \
      --discard-untrimmed \
      -p rev_R1.fastq.gz -o rev_R2.fastq.gz \
      ${input[0]} ${input[1]} \
      > cutadapt_2.log
    
    # cutadapt \
    #   -a FWDPRIMER...RCREVPRIMER \
    #   -A REVPRIMER...RCFWDPRIMER \
    #   --discard-untrimmed \
    #   -o out.1.fastq.gz -p out.2.fastq.gz \
    #   in.1.fastq.gz in.2.fastq.gz


    echo -e "\nReorienting"

    if [ -s for_R1.fastq.gz ]; then
      zcat for_R1.fastq.gz | seqkit replace -p "\\s.+" | gzip -7 > OK_R1.fastq.gz
      zcat for_R2.fastq.gz | seqkit replace -p "\\s.+" | gzip -7 > OK_R2.fastq.gz
    fi

    if [ -s rev_R1.fastq.gz ]; then
      echo -e "..Adding sequences to the main pool"
      zcat rev_R1.fastq.gz | seqkit replace -p "\\s.+" | gzip -7 >> OK_R1.fastq.gz
      zcat rev_R2.fastq.gz | seqkit replace -p "\\s.+" | gzip -7 >> OK_R2.fastq.gz

    else
      echo -e "..Probably all sequences are in forward orientation"
    fi


    echo -e "\nTrimming primers"
    if [ -s OK_R1.fastq.gz]; then

      cutadapt \
        -a ${params.primer_forward}";required;min_overlap=${params.primer_foverlap}"..."\$RR"";required;min_overlap=${params.primer_roverlap}" \
        --errors ${params.primer_mismatches} \
        --cores ${task.cpus} \
        --action=trim \
        --discard-untrimmed \
        --minimum-length ${params.trim_minlen} \
        --output ${sampID}_R1.fq.gz --paired-output ${sampID}_R2.fq.gz \
        OK_R1.fastq.gz OK_R2.fastq.gz

    fi


    ## Quality estimation and dereplication

    if [ -s ${sampID}_R1.fq.gz ]; then 

      ## Estimate sequence quality (for the extracted region)
      ## Sequence ID - Hash - Length - Average Phred score
      echo -e "\nCreating sequence hash table with average sequence quality"
      
      seqkit fx2tab --length --avg-qual ${sampID}_R1.fq.gz \
        | hash_sequences.sh \
        | awk '{print \$1 "\t" \$6 "\t" \$4 "\t" \$5}' \
        > tmp_hash_table_R1.txt

      seqkit fx2tab --length --avg-qual ${sampID}_R2.fq.gz \
        | hash_sequences.sh \
        | awk '{print \$1 "\t" \$6 "\t" \$4 "\t" \$5}' \
        > tmp_hash_table_R2.txt
      
      echo -e "..Done"


      ## Estimating MaxEE
      echo -e "\nEstimating maximum number of expected errors per sequence"

      vsearch \
          --fastx_filter ${sampID}_R1.fq.gz \
          --fastq_qmax 93 \
          --eeout \
          --fastaout - \
        | seqkit seq --name \
        | sed 's/;ee=/\t/g' \
        > tmp_ee_R1.txt

      vsearch \
          --fastx_filter ${sampID}_R2.fq.gz \
          --fastq_qmax 93 \
          --eeout \
          --fastaout - \
        | seqkit seq --name \
        | sed 's/;ee=/\t/g' \
        > tmp_ee_R2.txt

      echo -e "..Done"


      echo -e "\nMerging quality estimates"

      max_ee.R \
        tmp_hash_table_R1.txt \
        tmp_ee_R1.txt \
        ${sampID}_hash_table_R1.txt

      max_ee.R \
        tmp_hash_table_R2.txt \
        tmp_ee_R2.txt \
        ${sampID}_hash_table_R2.txt

      echo -e "..Done"


      ## Independent dereplication of pair-end reads
      echo -e "\nDereplicating R1 and R2 (independently)"
      
      seqkit fq2fa -w 0 ${sampID}_R1.fq.gz \
      | vsearch \
        --derep_fulllength - \
        --output - \
        --strand both \
        --fasta_width 0 \
        --threads 1 \
        --relabel_sha1 \
        --sizein --sizeout \
        --uc ${sampID}_uc_R1.uc \
        --quiet \
      | gzip -7 \
      > ${sampID}_R1.fa.gz

      seqkit fq2fa -w 0 ${sampID}_R2.fq.gz \
      | vsearch \
        --derep_fulllength - \
        --output - \
        --strand both \
        --fasta_width 0 \
        --threads 1 \
        --relabel_sha1 \
        --sizein --sizeout \
        --uc ${sampID}_uc_R2.uc \
        --quiet \
      | gzip -7 \
      > ${sampID}_R2.fa.gz


      echo -e "..Done"

      ## Compress results
      echo -e "Compressing result"
      gzip -7 ${sampID}_hash_table_R1.txt
      gzip -7 ${sampID}_hash_table_R2.txt
      gzip -7 ${sampID}_uc_R1.uc
      gzip -7 ${sampID}_uc_R2.uc


    else
      echo -e "\nNo sequences found after primer removal"
    fi

    ## Clean up
    if [ -f for_R1.fastq.gz ]; then rm for_R1.fastq.gz; fi
    if [ -f for_R2.fastq.gz ]; then rm for_R2.fastq.gz; fi
    if [ -f rev_R1.fastq.gz ]; then rm rev_R1.fastq.gz; fi
    if [ -f rev_R2.fastq.gz ]; then rm rev_R2.fastq.gz; fi
    if [ -f OK_R1.fastq.gz  ]; then rm OK_R1.fastq.gz; fi
    if [ -f OK_R2.fastq.gz  ]; then rm OK_R2.fastq.gz; fi
  
    echo -e "..Done"

    """
}




// Combine paired reads into single sequences 
// by reverse-complementing the reverse read and inserting poly-N padding
// + Estimate sequence qualities (without N pads!)
process join_pe {

    label "main_container"

    // publishDir "${out_1_joinPE}", mode: "${params.storagemode}"
    // cpus 2

    // Add sample ID to the log file
    tag "${input}"

    input:
      val input       // Sample name "(e.g., Barcode07_1__IS859)"
      path all_samples

    output:
      path "${input}_JoinedPE.fq.gz", emit: jj_FQ, optional: true
      path "${input}_JoinedPE_hash_table.txt.gz", emit: jj_hashes, optional: true

    script:
    sampID="${input}"
    
    """
    echo -e "Joining non-merged Illumina reads"
    echo -e "Input sample: " ${sampID}

    echo -e "\\nJoining with N-pads"
    vsearch \
      --fastq_join ${input}.R1.fastq.gz \
      --reverse ${input}.R2.fastq.gz \
      --join_padgap ${params.illumina_joinpadgap} \
      --join_padgapq ${params.illumina_joinpadqual} \
      --fastqout - \
    | seqkit replace -p "\\s.+" \
    | gzip -${params.gzip_compression} \
    > ${sampID}_JoinedPE.fq.gz

    ## Check if there are some sequences in the file
    if [ -n "\$(find . -name ${sampID}_JoinedPE.fq.gz -prune -size +29c)" ]; then

      echo -e "\\nJoining without N-pads (for quality estimation)"
      vsearch \
        --fastq_join ${input}.R1.fastq.gz \
        --reverse ${input}.R2.fastq.gz \
        --join_padgap "" \
        --join_padgapq "" \
        --fastqout - \
      | seqkit replace -p "\\s.+" \
      | gzip -${params.gzip_compression} \
      > tmp_for_qual.fq.gz


      ## Estimate sequence quality (without N pads!)
      ## Sequence ID - Hash - Length - Average Phred score
      echo -e "\\nCreating sequence hash table with average sequence quality"
        
      seqkit fx2tab --length --avg-qual tmp_for_qual.fq.gz \
        | hash_sequences.sh \
        | awk '{print \$1 "\t" \$6 "\t" \$4 "\t" \$5}' \
        > tmp_hash_table.txt
      
      echo -e "..Done"

      ## Estimating MaxEE
      echo -e "\\nEstimating maximum number of expected errors per sequence"

      vsearch \
          --fastx_filter tmp_for_qual.fq.gz \
          --fastq_qmax 93 \
          --eeout \
          --fastaout - \
        | seqkit seq --name \
        | sed 's/;ee=/\t/g' \
        > tmp_ee.txt

      echo -e "..Done"

      echo -e "\\nMerging quality estimates"

      max_ee.R \
        tmp_hash_table.txt \
        tmp_ee.txt \
        ${sampID}_JoinedPE_hash_table.txt

      echo -e "..Done"

      ## Compress results
      gzip -${params.gzip_compression} ${sampID}_JoinedPE_hash_table.txt

      ## Clean up
      rm tmp_for_qual.fq.gz
      rm tmp_hash_table.txt
      rm tmp_ee.txt
    
    else
      echo -e "\\nIt looks like there are no joined reads"
    fi

    ## Remove redundant symlinks
    find -L . -name "*.fastq.gz" | grep -v ${input} | parallel -j1 "rm {}"

    """
}

