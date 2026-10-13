process TRAIN_THRESHOLDS {
    tag "${meta.id}"
    label 'optimotu'

    input:
    tuple val(meta), path(seq_ids)
    path training_fasta
    path training_taxonomy
    val opts

    output:
    tuple val(meta), path("${meta.id}.thresholds.tsv"),  emit: thresholds
    tuple val(meta), path("${meta.id}.thresholds.json"), emit: sidecar

    script:
    def penalties  = opts.dist_penalties ? "--dist-penalties ${opts.dist_penalties}" : ''
    def orphans    = opts.train_keep_orphan_ranks ? '--keep-orphan-ranks' : ''
    def label_arg  = opts.label ? "--label '" + opts.label.toString().replace("'", "'\\''") + "'" : ''
    // Chunks read the taxonomy REFERENCE_PREP already cleaned and homonym-qualified across the whole database;
    // re-qualifying inside a chunk would see only its own names and could resolve them differently
    def chunking = meta.own ? """--seq-ids ${seq_ids} \\
        --chunk-id ${meta.id} \\
        --own-superranks ${meta.own} \\
        --no-qualify-homonyms""" : ''
    """
    optimotu_train.R \\
        --training-fasta ${training_fasta} \\
        --training-taxonomy ${training_taxonomy} \\
        --output-prefix ${meta.id} \\
        --ranks ${meta.ranks} \\
        --root-taxon ${opts.root_taxon} \\
        --metric ${opts.train_metric} \\
        --dist-max ${opts.dist_max} \\
        --dist-step ${opts.dist_step} \\
        --min-taxa ${opts.min_taxa} \\
        --min-refseq ${opts.min_refseq} \\
        --dist-method ${opts.dist_method} \\
        --clust-method ${opts.clust_method} \\
        --threads ${task.cpus} \\
        ${chunking} \\
        ${penalties} ${orphans} ${label_arg} \\
        --overwrite
    """
}
