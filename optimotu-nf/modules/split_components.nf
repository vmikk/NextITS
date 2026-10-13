process SPLIT_COMPONENTS {
    tag "${meta.id}"
    label 'process_medium'
    label 'optimotu'

    input:
    tuple val(meta), path(fasta)
    path thresholds
    val opts

    output:
    tuple val(meta), path("components/*.fasta"),         emit: components
    tuple val(meta), path("${meta.id}.components.tsv"),  emit: manifest
    tuple val("${task.process}"), val('vsearch'), eval(
        "vsearch --version 2>&1 | sed -n '1s/^vsearch v\\([0-9.]*\\).*/\\1/p'"
    ), topic: versions

    script:
    def threshold_arg = thresholds ? "--thresholds ${thresholds}" : ''
    """
    optimotu_split_components.R \\
        --input-fasta ${fasta} \\
        --output-dir components \\
        --max-dist ${opts.component_max_dist} \\
        --prefix ${meta.id} \\
        --threads ${task.cpus} \\
        ${threshold_arg} \\
        --overwrite
    mv components/components.tsv ${meta.id}.components.tsv
    """
}
