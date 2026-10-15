process FILTER_THRESHOLDS {
    label 'process_low'
    label 'optimotu'

    input:
    path thresholds, stageAs: 'unfiltered/*'
    path sidecar,    stageAs: 'unfiltered/*'
    path training_taxonomy
    val opts

    output:
    path "reference.thresholds.tsv",             emit: thresholds
    path "reference.thresholds.json",            emit: sidecar
    path "reference.thresholds.diagnostics.tsv", emit: diagnostics

    script:
    // Every rule is off unless its option is set; with none, the table passes
    // through unchanged
    def values = [
        'min-taxa-ge2'      : opts.filter_min_taxa_ge2,
        'min-frac-taxa-ge2' : opts.filter_min_frac_taxa_ge2,
        'min-score'         : opts.filter_min_score,
        'max-dominance'     : opts.filter_max_dominance,
        'max-measure-gap'   : opts.filter_max_measure_gap,
        'measure-gap-ref'   : opts.filter_measure_gap_ref,
        'measure-gap-action': opts.filter_measure_gap_action,
        'max-threshold'     : opts.filter_max_threshold,
        'min-threshold'     : opts.filter_min_threshold,
        'above-rank'        : opts.filter_above_rank,
        'min-refseq-above'  : opts.filter_min_refseq_above,
        'min-taxa-above'    : opts.filter_min_taxa_above,
        'min-coverage'      : opts.filter_min_coverage,
        'monotonic'         : opts.filter_monotonic,
        'root-thresholds'   : opts.filter_root_thresholds,
        'root-threshold-values': opts.filter_root_threshold_values
    ]
    def rules = values.findAll { _key, value -> value != null && value.toString() != '' }
        .collect { key, value -> "--${key} '${value}'" }
        .join(' ')
    def edge = opts.filter_drop_grid_edge ? '--drop-grid-edge' : ''
    def zero = opts.filter_drop_grid_zero ? '--drop-grid-zero' : ''
    """
    optimotu_filter_thresholds.R \\
        --thresholds ${thresholds} \\
        --training-taxonomy ${training_taxonomy} \\
        ${rules} ${edge} ${zero} \\
        --output-prefix reference \\
        --overwrite
    """
}
