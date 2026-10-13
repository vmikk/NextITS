process REFERENCE_SET {
    label 'process_high_memory'
    label 'optimotu'

    input:
    path thresholds
    path sidecar
    path training_taxonomy
    path classifier_fasta
    val opts

    output:
    path "reference_set/thresholds.tsv",             emit: thresholds
    path "reference_set/thresholds.json",            emit: sidecar
    path "reference_set/classifier.sintax.fasta.gz", emit: classifier_db
    path "reference_set/unreachable_thresholds.tsv", emit: unreachable
    path "reference_set/manifest.json",              emit: manifest

    script:
    def fields = [
        label            : opts.label,
        id               : opts.refset_id,
        title            : opts.refset_title,
        database         : opts.refset_database,
        'database-version': opts.refset_database_version,
        region           : opts.refset_region,
        citation         : opts.refset_citation,
        'index-name'     : opts.index_name
    ]
    // free-text fields (titles, citations) go through single quotes
    def metadata = fields.findAll { _key, value -> value }
        .collect { key, value -> "--${key} '" + value.toString().replace("'", "'\\''") + "'" }
        .join(' ')
    def orphans = opts.train_keep_orphan_ranks ? '--keep-orphan-ranks' : ''
    """
    optimotu_reference_set.R \\
        --fasta ${classifier_fasta} \\
        --training-taxonomy ${training_taxonomy} \\
        --thresholds ${thresholds} \\
        --taxonomy-format ${opts.taxonomy_format} \\
        --rank-letters ${opts.rank_letters} \\
        ${metadata} ${orphans} \\
        --output-dir reference_set \\
        --overwrite
    """
}
