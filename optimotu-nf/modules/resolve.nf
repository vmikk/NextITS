// First pass of a taxonomic split: the rank loop up to opts.split_rank
//
// Writes the sequences it finishes (clustered de novo at or above the split rank)
// in the usual output tables, named <id>.p000.*, and every other sequence into
// per-taxon parts, <id>.parts/pNNN.{fasta,taxprob.tsv}, for OPTIMOTU_CLUSTER
process OPTIMOTU_RESOLVE {
    tag "${meta.id}"
    label 'process_high_memory'
    label 'optimotu'

    input:
    tuple val(meta), path(fasta), path(taxprob)
    path thresholds
    path sidecar
    val opts

    output:
    tuple val(meta), path("${meta.id}.p000.membership.tsv"),        emit: membership
    tuple val(meta), path("${meta.id}.p000.otus.fasta"),            emit: representatives
    tuple val(meta), path("${meta.id}.p000.otu_taxonomy.tsv"),      emit: otu_taxonomy
    tuple val(meta), path("${meta.id}.p000.sequence_taxonomy.tsv"), emit: sequence_taxonomy
    tuple val(meta), path("${meta.id}.p000.uc"), optional: true,    emit: uc
    tuple val(meta), path("${meta.id}.p000.run.json"),              emit: run_record
    tuple val(meta), path("${meta.id}.parts/p*.fasta", arity: '0..*'), path("${meta.id}.parts/p*.taxprob.tsv", arity: '0..*'), emit: parts
    tuple val("${task.process}"), val('optimotu'), eval(
        "Rscript -e 'cat(as.character(packageVersion(\"optimotu\")))'"
    ), topic: versions

    script:
    def penalties  = opts.dist_penalties ? "--dist-penalties ${opts.dist_penalties}" : ''
    def uc         = opts.emit_uc ? '' : '--no-uc'
    def closed_ref = opts.closed_ref_from ? "--closed-ref-from ${opts.closed_ref_from}" : ''
    """
    optimotu_cluster.R \\
        --input-fasta ${fasta} \\
        --thresholds ${thresholds} \\
        --tax-prob ${taxprob} \\
        --output-prefix ${meta.id}.p000 \\
        --chunk-id ${meta.id}.p000 \\
        --split-rank ${opts.split_rank} \\
        --split-target ${opts.split_target} \\
        --split-dir ${meta.id}.parts \\
        --metric ${opts.metric} \\
        --confidence ${opts.confidence} \\
        --force-denovo ${opts.force_denovo} \\
        --ranks ${opts.ranks} \\
        --dist-method ${opts.dist_method} \\
        --clust-method ${opts.clust_method} \\
        --threads ${task.cpus} \\
        ${penalties} ${uc} ${closed_ref} \\
        --overwrite
    """
}
