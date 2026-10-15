include { REFERENCE_PREP      } from '../modules/reference_prep'
include { PARTITION_REFERENCE } from '../modules/partition_reference'
include { TRAIN_THRESHOLDS    } from '../modules/train_thresholds'
include { MERGE_THRESHOLDS    } from '../modules/merge_thresholds'
include { FILTER_THRESHOLDS   } from '../modules/filter_thresholds'
include { SINTAXER_INDEX      } from '../modules/sintaxer_index'
include { SINTAXER_INDEX as REFERENCE_SET_INDEX } from '../modules/sintaxer_index'
include { REFERENCE_SET       } from '../modules/reference_set'
include { optimotuOptions     } from './optimotu_options'

/*
 * Build a training set from a reference database and optimize taxon-specific
 * clustering thresholds on it.
 *
 * Threshold optimization is all-against-all, so a single process is quadratic
 * in the size of the reference. With the train_partition option the reference
 * is instead split into taxonomic subtrees that are trained independently (and
 * concurrently) and merged, which is both faster and more accurate: every
 * taxon-specific threshold is trained on all of that taxon's sequences rather
 * than on a subsample of the whole database.
 *
 * The trained table is then filtered (FILTER_THRESHOLDS, option
 * filter_thresholds): rows whose training subset cannot support them are
 * dropped, so that clustering falls back to an ancestor's threshold. Every
 * filter rule is off unless its option is set, and the table then passes
 * through unchanged.
 *
 * The thresholds are then bundled with a classifier database whose taxon
 * names match them (REFERENCE_SET, option train_reference_set): the
 * train_classifier_reference FASTA, or the training reference, rewritten with
 * the training-time names. A reference set directory is what --reference_set
 * takes in cluster mode.
 *
 * Each job's meta carries `pairs`, the estimated number of pairwise alignments
 * (null for a single unpartitioned job), so the configuration can give large
 * chunks more threads than small ones.
 *
 * take:
 *   reference  reference FASTA (Path or String)
 *   taxonomy   taxonomy TSV (Path, String or null to parse the FASTA headers)
 *   options    map of option overrides, see optimotu_options.nf
 */
workflow OPTIMOTU_TRAIN {

    take:
    reference
    taxonomy
    options

    main:
    def opts = optimotuOptions(options)
    def reference_file = file(reference, checkIfExists: true)
    def taxonomy_file = taxonomy ? file(taxonomy, checkIfExists: true) : []
    // recorded in the threshold sidecars and the reference set manifest
    def run_label = opts.label ?:
        reference_file.name.replaceFirst(/(\.sintax)?\.(fa|fasta|fas|fna)(\.gz|\.zst)?$/, '')

    REFERENCE_PREP(
        reference_file,
        taxonomy_file,
        opts.subMap(['taxonomy_format', 'ranks', 'root_taxon', 'rank_letters', 'train_min_per_taxon',
                     'train_max_per_taxon', 'train_max_seqs', 'train_seed', 'train_restrict',
                     'train_require_rank', 'train_emit_sintax', 'train_keep_orphan_ranks',
                     'train_partition'])
    )

    // The SINTAX database serves --classifier sintax as it is, and sintaxer (the
    // default classifier) through an index built from it once, here
    if (opts.train_emit_sintax) {
        SINTAXER_INDEX(REFERENCE_PREP.out.sintax_db, '')
        ch_sintax_index = SINTAXER_INDEX.out.index
    }
    else {
        ch_sintax_index = channel.empty()
    }

    if (opts.train_partition) {
        PARTITION_REFERENCE(
            REFERENCE_PREP.out.training_taxonomy,
            opts.subMap(['ranks', 'train_chunk_budget_pairs', 'train_chunk_subsample_max',
                         'train_chunk_min_per_taxon', 'min_taxa', 'min_refseq',
                         'train_min_species_seqs', 'train_seed'])
        )

        // Pair every planned job with its own ID list.  The plan names the file;
        // resolving it against the collected paths keeps the two in step even
        // though Nextflow stages them as a flat list.
        ch_ids = PARTITION_REFERENCE.out.seq_ids
            .flatten()
            .map { ids -> [ids.name, ids] }

        ch_jobs = PARTITION_REFERENCE.out.partition
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                [
                    row.ids_file,
                    [
                        id:    row.chunk_id,
                        kind:  row.job_kind,
                        ranks: row.chunk_ranks,
                        own:   row.own_superranks,
                        pairs: row.est_pairs.toDouble().toLong()
                    ]
                ]
            }
            .join(ch_ids)
            .map { _name, meta, ids -> [meta, ids] }
    } else {
        // One job covering the whole reference; meta.own stays null, which is
        // what tells optimotu_train.R it is not a chunk.
        //
        // The subsetting defaults assume partitioning bounds the cost, so
        // turning it off without also capping the training set asks one
        // process to cluster the whole database all-against-all.
        def cap = opts.train_max_seqs.toString().toLong()
        if (cap > 20000) {
            // 1.0e-4 s per pair, measured on full-length ITS with WFA2
            def hours = Math.round(cap * (cap - 1) / 2 * 1.0e-4 / 3600)
            log.warn "--train_partition false trains in one process, which is " +
                "quadratic in the training set size, but --train_max_seqs allows up " +
                "to ${cap} sequences -- about ${hours} CPU-hours at that size. " +
                "Cap it (5000 was the old default) or leave partitioning on."
        }
        ch_jobs = channel.of([[id: 'reference', kind: 'monolithic',
                               ranks: opts.ranks, own: null, pairs: null], []])
    }

    TRAIN_THRESHOLDS(
        ch_jobs,
        REFERENCE_PREP.out.training_fasta,
        REFERENCE_PREP.out.training_taxonomy,
        opts.subMap(['ranks', 'root_taxon', 'train_metric', 'dist_max', 'dist_step', 'min_taxa', 'min_refseq',
                     'train_min_species_seqs', 'train_drop_unpaired', 'dist_method', 'dist_penalties', 'clust_method',
                     'train_keep_orphan_ranks']) +
            [label: run_label]
    )

    if (opts.train_partition) {
        // Two independent .collect() calls can order their lists differently,
        // which would pair each table with another chunk's sidecar.  Sorting
        // both by file name lines them up, since the stems match; the merge
        // script re-checks the pairing rather than trusting it.
        ch_tsv = TRAIN_THRESHOLDS.out.thresholds
            .map { _meta, tsv -> tsv }
            .toSortedList { a, b -> a.name <=> b.name }
        ch_json = TRAIN_THRESHOLDS.out.sidecar
            .map { _meta, json -> json }
            .toSortedList { a, b -> a.name <=> b.name }

        MERGE_THRESHOLDS(ch_tsv, ch_json, PARTITION_REFERENCE.out.expected_rows,
                         opts.subMap(['train_keep_underpowered']) + [label: run_label])
        ch_thresholds = MERGE_THRESHOLDS.out.thresholds
        ch_sidecar    = MERGE_THRESHOLDS.out.sidecar
    } else {
        ch_thresholds = TRAIN_THRESHOLDS.out.thresholds.map { _meta, tsv -> tsv }
        ch_sidecar    = TRAIN_THRESHOLDS.out.sidecar.map    { _meta, json -> json }
    }

    // Everything downstream (the reference set, the emitted thresholds) takes
    // the filtered table; the unfiltered one is published beside it
    if (opts.filter_thresholds) {
        FILTER_THRESHOLDS(
            ch_thresholds,
            ch_sidecar,
            REFERENCE_PREP.out.training_taxonomy,
            opts.subMap(['filter_min_taxa_ge2', 'filter_min_frac_taxa_ge2', 'filter_min_score',
                         'filter_max_dominance', 'filter_max_measure_gap', 'filter_measure_gap_ref',
                         'filter_measure_gap_action', 'filter_max_threshold', 'filter_min_threshold',
                         'filter_drop_grid_edge', 'filter_drop_grid_zero', 'filter_above_rank',
                         'filter_min_refseq_above', 'filter_min_taxa_above',
                         'filter_min_coverage', 'filter_monotonic',
                         'filter_root_thresholds', 'filter_root_threshold_values'])
        )
        ch_thresholds = FILTER_THRESHOLDS.out.thresholds
        ch_sidecar    = FILTER_THRESHOLDS.out.sidecar
    }

    // A classifier database rewritten from the training taxonomy, so that the
    // classifier returns the names the thresholds are keyed by
    def classifier_file = opts.train_classifier_reference
        ? file(opts.train_classifier_reference, checkIfExists: true)
        : reference_file
    def build_refset = opts.train_reference_set
    if (build_refset && taxonomy && !opts.train_classifier_reference) {
        log.warn "OPTIMOTU_TRAIN: the taxonomy comes from a separate table, so the " +
            "training reference headers cannot serve as a classifier database; " +
            "skipping the reference set (set train_classifier_reference to build one)"
        build_refset = false
    }
    if (build_refset) {
        def index_stem = 'classifier'
        REFERENCE_SET(
            ch_thresholds,
            ch_sidecar,
            REFERENCE_PREP.out.training_taxonomy,
            classifier_file,
            opts.subMap(['taxonomy_format', 'rank_letters', 'train_keep_orphan_ranks',
                         'refset_id', 'refset_title', 'refset_database',
                         'refset_database_version', 'refset_region', 'refset_citation']) +
                [label: run_label,
                 index_name: opts.train_reference_set_index ? "${index_stem}.stx" : null]
        )
        ch_refset_manifest = REFERENCE_SET.out.manifest
        ch_refset_db       = REFERENCE_SET.out.classifier_db
        if (opts.train_reference_set_index) {
            REFERENCE_SET_INDEX(REFERENCE_SET.out.classifier_db, index_stem)
            ch_refset_index = REFERENCE_SET_INDEX.out.index
        }
        else {
            ch_refset_index = channel.empty()
        }
    }
    else {
        ch_refset_manifest = channel.empty()
        ch_refset_db       = channel.empty()
        ch_refset_index    = channel.empty()
    }

    emit:
    thresholds        = ch_thresholds
    sidecar           = ch_sidecar
    training_fasta    = REFERENCE_PREP.out.training_fasta
    training_taxonomy = REFERENCE_PREP.out.training_taxonomy
    sintax_db         = REFERENCE_PREP.out.sintax_db
    sintax_index      = ch_sintax_index
    refset_manifest   = ch_refset_manifest
    refset_db         = ch_refset_db
    refset_index      = ch_refset_index
}
