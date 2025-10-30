
// Reference-based chimera removal
process chimera_ref {

    label "main_container"

    publishDir "${out_5_chim}", mode: "${params.storagemode}"
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
        echo -e "\\nBorderline sequences were added to non-chimeric sequences"
        cat borderline.fasta nonchimeras.fasta > nc_bo.fasta
        mv nc_bo.fasta nonchimeras.fasta
        rm borderline.fasta
    fi

    ## Chimeric sequences
    if [ -e chimeras.fasta ]
    then
      ## Add sample ID to the header and compress the file
      sed 's/>.*/&;sample='"${sampID}"';/' chimeras.fasta \
        | gzip -${params.gzip_compression} \
        > "${sampID}_Chimera.fa.gz"
      rm chimeras.fasta
    else
      echo -e "\\nNo chimeras detected"
    fi

    ## Non-chimeric sequences
    if [ -e nonchimeras.fasta ]
    then
      gzip -c nonchimeras.fasta > "${sampID}_NoChimera.fa.gz"
      rm nonchimeras.fasta
    else
      echo "No non-chimeric sequences left"
    fi

    """
}


// Recovery of ref-based chimeric sequences with high occurrence
process chimera_rescue {

    label "main_container"

    publishDir "${out_5_chim}", mode: "${params.storagemode}"
    // cpus 1

    input:
      path input

    output:
      path "*_RescuedChimera.fa.gz", emit: rescuedchimeric, optional: true

    script:
    """

    ## Aggregate chimeric sequences from different samples
    echo -e "\\nAggregating chimeric sequences"
    find . -name "*_Chimera.fa.gz" \
      | parallel -j1 "zcat {}" \
      | seqkit fx2tab \
      | sed -r 's:\t+:\t:g' | sed 's/\t\$//g' \
      | gzip -${params.gzip_compression} > All_chimeras.txt.gz
    echo -e "..Done"

    ### Inspect chimerae occurrence
    ## Rescue sequences that were annotated as chimeric, 
    ## but have high occurrence within sequenceing run (e.g., occurrence > 2)
    echo -e "\\nInspecting occurrence of chimeric sequences"

    chimera_rescue.R \
      "All_chimeras.txt.gz" \
      ${params.chimera_rescueoccurrence} \
      "Rescued_Chimeric_sequences.fa.gz"

    echo -e "..Done"

    ## Split rescured sequences by sample
    if [ -e Rescued_Chimeric_sequences.fa.gz ]
    then
      echo -e "\\n..Splitting rescued sequences by sample"
      seqkit split -i \
        --id-regexp ";sample=(.*);" \
        --threads ${task.cpus} \
        -w 0 \
        -O Rescued_by_sample \
        Rescued_Chimeric_sequences.fa.gz

      rename \
        --filename 's/^Rescued_Chimeric_sequences.id_//g; s/.fa.gz/_RescuedChimera.fa.gz/' \
        \$(find Rescued_by_sample -name "*.fa.gz")

      mv Rescued_by_sample/*_RescuedChimera.fa.gz .

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

    publishDir "${out_5_chim}", mode: "${params.storagemode}"
    // cpus 1

    // Add sample ID to the log file
    tag "${input.getSimpleName().replaceAll(/_Homopolymer_compressed/, '')}"

    input:
      path input

    output:
      path "${input.getSimpleName().replaceAll(/_Homopolymer_compressed/, '')}_DeNovoChim.txt", emit: denovochim, optional: true
      tuple val("${task.process}"), val('vsearch'), eval('vsearch --version 2>&1 | head -n 1 | sed "s/vsearch //g" | sed "s/,.*//g" | sed "s/^v//" | sed "s/_.*//"'), topic: versions

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
      path "DeNovo_Chimera.txt", emit: alldenovochim, optional: true

    script:
    """
    echo -e "Aggregating de novo chimeric sequences"

    find . -name "*_DeNovoChim.txt" \
      | parallel -j1 "cat {}" \
      > DeNovo_Chimera.txt

    echo -e "..Done"

    """
}

