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
    def unpaired   = opts.train_drop_unpaired ? '--drop-unpaired' : ''
    def label_arg  = opts.label ? "--label '" + opts.label.toString().replace("'", "'\\''") + "'" : ''
    // Chunks read the taxonomy REFERENCE_PREP already cleaned and
    // homonym-qualified across the whole database; re-qualifying inside a chunk
    // would see only its own names and could resolve them differently.
    // A chunk's first rank is its own top rank; --drop-unpaired must only spare
    // the rows of the reference's root rank.
    def root_rank = opts.ranks.toString().tokenize(',')[0].trim()
    def chunking = meta.own ? """--seq-ids ${seq_ids} \\
        --chunk-id ${meta.id} \\
        --own-superranks ${meta.own} \\
        --reference-root-rank ${root_rank} \\
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
        --min-species-seqs ${opts.train_min_species_seqs ?: 1} \\
        --dist-method ${opts.dist_method} \\
        --clust-method ${opts.clust_method} \\
        --threads ${task.cpus} \\
        ${chunking} \\
        ${penalties} ${orphans} ${unpaired} ${label_arg} \\
        --overwrite
    """
}
