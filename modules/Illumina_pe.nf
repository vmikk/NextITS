

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

