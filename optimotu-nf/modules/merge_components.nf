// A `path` input arrives as a List for several files but as a bare Path for
// one, so normalize before building the comma-separated argument.
def asCsv(files) {
    return (files instanceof List ? files : [files]).join(',')
}

process MERGE_COMPONENTS {
    tag "${meta.id}"
    label 'process_low'
    label 'optimotu'

    // Inputs are staged into their own directories, so a component can never
    // collide with the merged output names
    input:
    tuple val(meta), path(membership, stageAs: 'membership/*'), path(otu_taxonomy, stageAs: 'otu_taxonomy/*'),
          path(otus_fasta, stageAs: 'otus/*'), path(uc, stageAs: 'uc/*')

    output:
    tuple val(meta), path("${meta.id}.membership.tsv"),        emit: membership
    tuple val(meta), path("${meta.id}.otu_taxonomy.tsv"),      emit: otu_taxonomy
    tuple val(meta), path("${meta.id}.otus.fasta"),            emit: representatives
    tuple val(meta), path("${meta.id}.uc"), optional: true,    emit: uc

    script:
    def uc_arg = uc ? "--uc " + asCsv(uc) : ''
    """
    optimotu_merge.R \\
        --membership ${asCsv(membership)} \\
        --otu-taxonomy ${asCsv(otu_taxonomy)} \\
        --otus-fasta ${asCsv(otus_fasta)} \\
        ${uc_arg} \\
        --output-prefix ${meta.id} \\
        --overwrite
    """
}
