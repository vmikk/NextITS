process REFERENCE_PREP {
    label 'process_high_memory'
    label 'optimotu'

    input:
    path reference
    path taxonomy
    val opts

    output:
    path "reference.training.fasta",         emit: training_fasta
    path "reference.training_taxonomy.tsv",  emit: training_taxonomy
    path "reference.training.json",          emit: provenance
    path "reference.sintax.fasta", optional: true, emit: sintax_db

    script:
    def tax_arg  = taxonomy ? "--taxonomy-tsv ${taxonomy}" : ''
    def restrict = opts.train_restrict ? "--restrict ${opts.train_restrict}" : ''
    def sintax   = opts.train_emit_sintax ? '--emit-sintax' : ''
    def orphans  = opts.train_keep_orphan_ranks ? '--keep-orphan-ranks' : ''
    // With partitioned training the subsampling that keeps the cost tractable is
    // decided per chunk by optimotu_partition_reference.R, so this step must let
    // the whole reference through -- including the records identified only to
    // order or class, which the higher-rank chunks are trained on.  The default
    // of dropping everything unassigned at the tip rank would discard most of a
    // database like UNITE.
    def root_rank    = opts.ranks.toString().tokenize(',')[0].trim()
    def rank_default = opts.train_partition ? root_rank : null
    def required     = opts.train_require_rank ?: rank_default
    def require_rank = required ? "--require-rank ${required}" : ''
    """
    optimotu_reference_prep.R \\
        --reference-fasta ${reference} \\
        --output-prefix reference \\
        --taxonomy-format ${opts.taxonomy_format} \\
        --ranks ${opts.ranks} \\
        --root-taxon ${opts.root_taxon} \\
        --rank-letters ${opts.rank_letters} \\
        --min-per-taxon ${opts.train_min_per_taxon} \\
        --max-per-taxon ${opts.train_max_per_taxon} \\
        --max-seqs ${opts.train_max_seqs} \\
        --seed ${opts.train_seed} \\
        ${tax_arg} ${restrict} ${require_rank} ${sintax} ${orphans} \\
        --overwrite
    """
}
