process PARTITION_REFERENCE {
    label 'process_low'
    label 'optimotu'

    input:
    path training_taxonomy
    val opts

    output:
    path "partition/partition.tsv",     emit: partition
    path "partition/expected_rows.tsv", emit: expected_rows
    path "partition/partition.json",    emit: provenance
    path "partition/*.ids.txt",         emit: seq_ids

    script:
    """
    optimotu_partition_reference.R \\
        --training-taxonomy ${training_taxonomy} \\
        --output-dir partition \\
        --ranks ${opts.ranks} \\
        --budget-pairs ${opts.train_chunk_budget_pairs} \\
        --subsample-max ${opts.train_chunk_subsample_max} \\
        --subsample-min-per-taxon ${opts.train_chunk_min_per_taxon} \\
        --min-taxa ${opts.min_taxa} \\
        --min-refseq ${opts.min_refseq} \\
        --seed ${opts.train_seed} \\
        --overwrite
    """
}
