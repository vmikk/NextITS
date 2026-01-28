/*
============================================================================
  NextITS: Pipeline to process eukaryotic ITS amplicons
============================================================================
  License: Apache-2.0
  Github : https://github.com/vmikk/NextITS
  Website: https://Next-ITS.github.io/
----------------------------------------------------------------------------
*/

// Subworkflow for ITSx processing
// (which splits large dereplicated FASTAs into chunks).
// The workflow is as follows:
//   1. Trim primers and dereplicate at sample level
//   2. Split the dereplicated primer-trimmed sequences (at sample level) into chunks while preserving metadata
//   3. Run ITSx on each chunk
//   4. Group results back by original sample ID and concatenate + convert ITSx output to Parquet

// Path to the output results
out_3_itsx   = params.outdir + "/03_ITSx"
out_3_itsxp  = params.outdir + "/03_ITSx_PooledParts"


// Trim primers and dereplicate at sample level
process primer_trim {

    label "main_container"

    publishDir "${out_3_itsx}", mode: "${params.storagemode}"
    // cpus 2

    // Add sample ID to the log file
    tag "${meta.id}"

    input:
      tuple val(meta), path(fastq)

    output:
      tuple val(meta), path("${meta.id}_derep.fasta.gz"),    emit: derep,  optional: true
      tuple val(meta), path("${meta.id}_hash_table.txt.gz"), emit: hashes, optional: true
      tuple val(meta), path("${meta.id}_uc.uc.gz"),          emit: uc,     optional: true
      tuple val(meta), path("${meta.id}_primertrimmed_sorted.fq.gz"), emit: trimmed_seqs,   optional: true
      tuple val("${task.process}"), val('cutadapt'), eval('cutadapt --version'), topic: versions
      tuple val("${task.process}"), val('vsearch'), eval('vsearch --version 2>&1 | head -n 1 | sed "s/vsearch //g" | sed "s/,.*//g" | sed "s/^v//" | sed "s/_.*//"'), topic: versions
      tuple val("${task.process}"), val('seqkit'), eval('seqkit version | sed "s/seqkit v//"'), topic: versions
      tuple val("${task.process}"), val('phredsort'), eval('phredsort -v | sed "s/phredsort //"'), topic: versions
      tuple val("${task.process}"), val('seqhasher'), eval('seqhasher -v | sed "s/SeqHasher //"'), topic: versions
      tuple val("${task.process}"), val('parallel'), eval('parallel --version | head -n 1 | sed "s/GNU parallel //"'), topic: versions
      tuple val("${task.process}"), val('brename'), eval('brename --help | head -n 4 | tail -1 | sed "s/Version: //"'), topic: versions

    script:
    sampID="${meta.id}"
    """
    echo -e "Primer trimming and dereplication at sample level\\n"
    echo -e "Input sample: " ${sampID}

    ## Trim primers    
    echo -e "Trimming primers\\n"

    ## Reverse-complement rev primer
    RR=\$(rc.sh ${params.primer_reverse})

    cutadapt \
      -a ${params.primer_forward}";required;min_overlap=${params.primer_foverlap}"..."\$RR"";required;min_overlap=${params.primer_roverlap}" \
      --errors ${params.primer_mismatches} \
      --revcomp --rename "{id}" \
      --discard-untrimmed \
      --minimum-length ${params.trim_minlen} \
      --cores ${task.cpus} \
      --action trim \
      --output ${sampID}_primertrimmed.fq.gz \
      ${fastq}

    echo -e "..Done\\n"

    ## Check if there are sequences in the output
    NUMSEQS=\$( seqkit stat --tabular --quiet ${sampID}_primertrimmed.fq.gz | awk -F'\t' 'NR==2 {print \$4}' )
    echo -e "Number of sequences after primer trimming: " \$NUMSEQS
    if [ \$NUMSEQS -lt 1 ]; then
      echo -e "\\nIt looks like no reads remained after trimming the primers\\n"
      exit 0
    fi
   
    ## Estimate sequence quality and sort sequences by quality
    echo -e "\\nSorting by sequence quality"
    seqkit replace -p "\\s.+" ${sampID}_primertrimmed.fq.gz \
      | phredsort -i - -o - --metric meep --header avgphred,maxee,meep \
      | gzip -1 > ${sampID}_primertrimmed_sorted.fq.gz
    echo -e "..Done"

    ## Hash sequences, add sample ID to the header
    ## columns: Sample ID - Hash - PacBioID - AvgPhredScore - MaxEE - MEEP - Sequence - Quality - Length
    ## Convert to Parquet format
    echo -e "\\nCreating hash table"
    seqhasher --hash sha1 --name ${sampID} ${sampID}_primertrimmed_sorted.fq.gz - \
      | seqkit fx2tab --length \
      | sed 's/;/\t/ ; s/;/\t/ ; s/ avgphred=/\t/ ; s/ maxee=/\t/ ; s/ meep=/\t/' \
      > ${sampID}_hash_table.txt
    echo -e "..Done"

    ## Check the number of fields per record (should be 9!)
    # awk '{print NF}' ${sampID}_hash_table.txt | sort | uniq -c
    # awk 'NF > 9 {print \$0 }' ${sampID}_hash_table.txt

    ## Dereplicate at sample level (use quality-sorted sequences to make sure that the representative sequence is with the highest quality)
    echo -e "\\nDereplicating at sample level"
    seqkit fq2fa -w 0 ${sampID}_primertrimmed_sorted.fq.gz \
      | vsearch \
        --derep_fulllength - \
        --output - \
        --strand both \
        --fasta_width 0 \
        --threads 1 \
        --relabel_sha1 \
        --sizein --sizeout \
        --minseqlength ${params.trim_minlen} \
        --uc ${sampID}_uc.uc \
        --quiet \
      > ${sampID}_derep.fasta
    
    ## Remove temporary file
    rm ${sampID}_primertrimmed.fq.gz

    ## Compress results
    echo -e "\\nCompressing results"
    parallel -j${task.cpus} "gzip -${params.gzip_compression} {}" ::: \
      ${sampID}_hash_table.txt \
      ${sampID}_uc.uc \
      ${sampID}_derep.fasta

    echo -e "..Done"

    """
}
