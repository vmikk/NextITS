process MERGE_THRESHOLDS {
    label 'process_low'
    label 'optimotu'

    input:
    path thresholds, stageAs: 'chunks/*'
    path sidecars,   stageAs: 'chunks/*'
    path expected_rows
    val opts

    output:
    path "reference.thresholds.tsv",            emit: thresholds
    path "reference.thresholds.json",           emit: sidecar
    path "reference.thresholds.provenance.tsv", emit: provenance

    script:
    def tsv_list  = thresholds instanceof List ? thresholds.join(',') : "${thresholds}"
    def json_list = sidecars   instanceof List ? sidecars.join(',')   : "${sidecars}"
    def keep      = opts.train_keep_underpowered ? '--keep-underpowered' : ''
    def label_arg = opts.label ? "--label '" + opts.label.toString().replace("'", "'\\''") + "'" : ''
    """
    optimotu_merge_thresholds.R \\
        --thresholds '${tsv_list}' \\
        --sidecars '${json_list}' \\
        --expected-rows ${expected_rows} \\
        ${keep} ${label_arg} \\
        --output-prefix reference \\
        --overwrite
    """
}
