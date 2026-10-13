process SINTAXER_INDEX {
    tag "${db.name}"
    label 'optimotu'

    input:
    path db
    // index file name without .stx; empty: derived from the database name
    val stem

    output:
    path "*.stx", emit: index
    tuple val("${task.process}"), val('sintaxer'), eval("sintaxer --version | sed 's/^sintaxer //'"), topic: versions

    script:
    // `simpleName` would cut at the first dot (General_EUK_ITS_v2.1.fasta ->
    // General_EUK_ITS_v2), so only the FASTA and compression suffixes are dropped
    def name = stem ?: db.name.replaceFirst(/\.(fa|fasta|fas|fna)(\.gz|\.zst)?$/, '')
    // sintaxer writes the index through a private (0600) temporary file; a
    // published index is a shared database, so make it readable like the FASTA
    """
    sintaxer index ${db} --output ${name}.stx
    chmod 0644 ${name}.stx
    """
}
