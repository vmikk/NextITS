process CLASSIFY {
    tag "${meta.id}"
    label 'process_medium'
    label 'optimotu'

    input:
    tuple val(meta), path(fasta)
    path db
    val opts

    output:
    tuple val(meta), path("${meta.id}.taxprob.tsv"),  emit: taxprob
    tuple val(meta), path("${meta.id}.taxprob.json"), emit: provenance
    tuple val("${task.process}"), val(opts.classifier), eval(
        opts.classifier == 'sintaxer' ? "sintaxer --version | sed 's/^sintaxer //'" :
        opts.classifier == 'none'     ? "echo NA" :
        "vsearch --version 2>&1 | sed -n '1s/^vsearch v\\([0-9.]*\\).*/\\1/p'"
    ), topic: versions

    script:
    def db_arg = db ? "--db ${db}" : ''
    def placeholders = opts.keep_placeholders ? '--keep-placeholders' : ''
    """
    optimotu_classify.R \\
        --input-fasta ${fasta} \\
        --output ${meta.id}.taxprob.tsv \\
        --classifier ${opts.classifier} \\
        ${db_arg} \\
        --min-identity ${opts.min_identity} \\
        --ranks ${opts.ranks} \\
        --root-taxon ${opts.root_taxon} \\
        --rank-letters ${opts.rank_letters} \\
        --taxonomy-format ${opts.taxonomy_format} \\
        --threads ${task.cpus} \\
        --seed ${opts.classifier_seed} \\
        ${placeholders} \\
        --overwrite
    """
}
