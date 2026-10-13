process OPTIMOTU_CLUSTER {
    tag "${meta.id}"
    label 'process_high_memory'
    label 'optimotu'

    input:
    tuple val(meta), path(fasta), path(taxprob)
    path thresholds
    path sidecar
    val opts

    output:
    tuple val(meta), path("${meta.id}.membership.tsv"),        emit: membership
    tuple val(meta), path("${meta.id}.otus.fasta"),            emit: representatives
    tuple val(meta), path("${meta.id}.otu_taxonomy.tsv"),      emit: otu_taxonomy
    tuple val(meta), path("${meta.id}.sequence_taxonomy.tsv"), emit: sequence_taxonomy
    tuple val(meta), path("${meta.id}.uc"), optional: true,    emit: uc
    tuple val(meta), path("${meta.id}.run.json"),              emit: run_record
    tuple val("${task.process}"), val('optimotu'), eval(
        "Rscript -e 'cat(as.character(packageVersion(\"optimotu\")))'"
    ), topic: versions

    script:
    def penalties  = opts.dist_penalties ? "--dist-penalties ${opts.dist_penalties}" : ''
    def uc         = opts.emit_uc ? '' : '--no-uc'
    def closed_ref = opts.closed_ref_from ? "--closed-ref-from ${opts.closed_ref_from}" : ''
    def derep      = opts.dereplicate ? '--dereplicate' : ''
    """
    optimotu_cluster.R \\
        --input-fasta ${fasta} \\
        --thresholds ${thresholds} \\
        --tax-prob ${taxprob} \\
        --output-prefix ${meta.id} \\
        --chunk-id ${meta.id} \\
        --metric ${opts.metric} \\
        --confidence ${opts.confidence} \\
        --force-denovo ${opts.force_denovo} \\
        --ranks ${opts.ranks} \\
        --dist-method ${opts.dist_method} \\
        --clust-method ${opts.clust_method} \\
        --threads ${task.cpus} \\
        ${penalties} ${uc} ${closed_ref} ${derep} \\
        --overwrite
    """
}
