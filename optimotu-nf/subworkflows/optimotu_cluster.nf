include { SPLIT_COMPONENTS } from '../modules/split_components'
include { SINTAXER_INDEX   } from '../modules/sintaxer_index'
include { CLASSIFY         } from '../modules/classify'
include { OPTIMOTU_CLUSTER } from '../modules/cluster'
include { OPTIMOTU_RESOLVE } from '../modules/resolve'
include { MERGE_COMPONENTS } from '../modules/merge_components'
include { optimotuOptions  } from './optimotu_options'

/*
 * Taxonomically guided OTU clustering
 *
 * take:
 *   ch_fasta       channel of tuple(meta, fasta); meta.id names the outputs, and
 *                  every input is clustered and reported on its own.
 *                  FASTA may be gzipped; `;size=N` header annotations are used
 *                  as abundances and stripped from the sequence IDs
 *   thresholds     thresholds TSV (Path or String) from the training entry point;
 *                  a sidecar <name>.json next to it is picked up automatically
 *   classifier_db  classifier database (Path, String or null): SINTAX-format
 *                  FASTA, or for sintaxer also a prebuilt .stx index
 *   taxprob        precomputed seq_id/rank/taxon/prob TSV (Path, String or null)
 *                  covering every input; skips the classifier
 *   options        map of option overrides, see optimotu_options.nf
 *
 * emit (all tuple(meta, file), with the input's meta):
 *   membership, representatives, otu_taxonomy, uc
 *   sequence_taxonomy, taxprob, run_records -- per clustered target; with
 *   split_components or split_rank these are per part, and meta.parent is the
 *   input's meta
 *
 * With split_rank, each target is clustered in two steps that together give the
 * same OTUs as one run (see modules/resources/usr/bin/optimotu_loop.R):
 * OPTIMOTU_RESOLVE runs the rank loop down to split_rank on the whole target and
 * finishes what it clusters de novo there (part p000); every taxon of that rank
 * is then clustered further on its own, packed into parts of up to split_target
 * sequences, in parallel OPTIMOTU_CLUSTER tasks. MERGE_COMPONENTS joins the parts
 *
 * The processes carry the label 'optimotu' (container selection) and need
 * `nextflow.enable.moduleBinaries = true`, which puts modules/resources/usr/bin
 * on the task PATH
 */
workflow OPTIMOTU {

    take:
    ch_fasta
    thresholds
    classifier_db
    taxprob
    options

    main:
    def opts = optimotuOptions(options)

    if (!thresholds) {
        error("OPTIMOTU: a thresholds table is required (train one with the TRAIN entry point)")
    }
    def thresholds_file = file(thresholds, checkIfExists: true)
    // The sidecar records the distance configuration the thresholds are valid
    // for; the clustering script refuses a mismatch, so staging it alongside is
    // what makes that check possible.
    def sidecar_file = file(thresholds_file.toString().replaceFirst(/\.tsv$/, '.json'))
    def sidecar = (sidecar_file != thresholds_file && sidecar_file.exists()) ? sidecar_file : []

    def db_file = classifier_db ? file(classifier_db, checkIfExists: true) : null
    def db_is_index = db_file && db_file.name.endsWith('.stx')
    if (!taxprob && opts.classifier != 'none') {
        if (!db_file) {
            error("OPTIMOTU: classifier '${opts.classifier}' needs a reference database, " +
                  "or use classifier 'none' for de novo clustering")
        }
        // only sintaxer reads its own index; VSEARCH needs the FASTA
        if (opts.classifier in ['sintax', 'vsearch'] && db_is_index) {
            error("OPTIMOTU: classifier '${opts.classifier}' needs a FASTA database, " +
                  "but got a sintaxer index: ${db_file}")
        }
    }

    // Optional pre-split into connected components. Each component becomes its
    // own task; this is only sound when the components really are complete at a
    // cutoff no stricter than the thresholds' max_clustering_distance, which
    // SPLIT_COMPONENTS asserts for us. A component keeps its input's meta as
    // `parent`, and the number of siblings, for the regrouping below
    if (opts.split_components) {
        SPLIT_COMPONENTS(ch_fasta, thresholds_file, opts.subMap(['component_max_dist']))
        ch_targets = SPLIT_COMPONENTS.out.components
            .flatMap { meta, components ->
                def parts = components instanceof List ? components : [components]
                parts.collect { part ->
                    def id = part.name.replaceFirst(/\.fasta$/, '')
                    [meta + [id: id, parent: meta, n_parts: parts.size()], part]
                }
            }
    }
    else {
        ch_targets = ch_fasta
    }

    // A precomputed classification short-circuits the classifier entirely, which
    // is how a different classifier gets plugged in: emit the same four-column
    // table and pass it as `taxprob`
    if (taxprob) {
        def taxprob_file = file(taxprob, checkIfExists: true)
        ch_taxprob = ch_targets.map { meta, _fasta -> [meta, taxprob_file] }
        ch_calls = ch_targets.map { meta, fasta -> [meta, fasta, taxprob_file] }
    }
    else {
        // sintaxer classifies against a prebuilt index: 
        // a FASTA database is indexed once here, and the index (a value channel) is shared by every target
        // A .stx database skips this step
        if (opts.classifier == 'sintaxer' && !db_is_index) {
            SINTAXER_INDEX(db_file, '')
            ch_classifier_db = SINTAXER_INDEX.out.index
        }
        else {
            ch_classifier_db = db_file ?: []
        }
        CLASSIFY(
            ch_targets,
            ch_classifier_db,
            opts.subMap(['classifier', 'min_identity', 'ranks', 'root_taxon', 'rank_letters',
                         'taxonomy_format', 'classifier_seed', 'keep_placeholders'])
        )
        ch_taxprob = CLASSIFY.out.taxprob
        ch_calls = ch_targets.join(CLASSIFY.out.taxprob)
    }

    def cluster_opts = opts.subMap(['metric', 'confidence', 'force_denovo', 'ranks', 'dist_method',
                                    'dist_penalties', 'clust_method', 'emit_uc', 'dereplicate',
                                    'closed_ref_from'])
    def split_by_taxon = opts.split_rank != 'none'

    if (split_by_taxon) {
        OPTIMOTU_RESOLVE(ch_calls, thresholds_file, sidecar,
                         opts.subMap(['split_rank', 'split_target']) + cluster_opts)

        // One task per part, paired with its own classification table by name
        ch_cluster_calls = OPTIMOTU_RESOLVE.out.parts
            .flatMap { meta, fastas, taxprobs ->
                def parts = fastas instanceof List ? fastas : [fastas]
                def calls = (taxprobs instanceof List ? taxprobs : [taxprobs])
                    .collectEntries { t -> [t.simpleName, t] }
                parts.collect { part ->
                    [meta + [id: "${meta.id}.${part.simpleName}", parent: meta, n_parts: parts.size() + 1],
                     part, calls[part.simpleName]]
                }
            }

        // The first pass's own results are part p000 of its input
        ch_n_parts = OPTIMOTU_RESOLVE.out.parts
            .map { meta, fastas, _taxprobs -> [meta, (fastas instanceof List ? fastas.size() : 1) + 1] }
        ch_first = OPTIMOTU_RESOLVE.out.membership
            .join(OPTIMOTU_RESOLVE.out.otu_taxonomy)
            .join(OPTIMOTU_RESOLVE.out.representatives)
        ch_first = opts.emit_uc
            ? ch_first.join(OPTIMOTU_RESOLVE.out.uc)
            : ch_first.map { meta, m, t, r -> [meta, m, t, r, null] }
        ch_first = ch_first
            .join(ch_n_parts)
            .map { meta, m, t, r, u, n -> [meta + [id: "${meta.id}.p000", parent: meta, n_parts: n], m, t, r, u] }
    }
    else {
        ch_cluster_calls = ch_calls
        ch_first = channel.empty()
    }

    OPTIMOTU_CLUSTER(ch_cluster_calls, thresholds_file, sidecar, cluster_opts)

    if (opts.split_components || split_by_taxon) {
        // OTU IDs are representative sequence IDs, so they stay unique across parts and merging is concatenation.
        // Each input is merged as soon as all of its parts are done; sorting by file name keeps the merged files deterministic
        ch_parts = OPTIMOTU_CLUSTER.out.membership
            .join(OPTIMOTU_CLUSTER.out.otu_taxonomy)
            .join(OPTIMOTU_CLUSTER.out.representatives)
        if (opts.emit_uc) {
            ch_parts = ch_parts.join(OPTIMOTU_CLUSTER.out.uc)
        }
        else {
            ch_parts = ch_parts.map { meta, m, t, r -> [meta, m, t, r, null] }
        }
        ch_merge = ch_parts
            .mix(ch_first)
            .map { meta, m, t, r, u -> [groupKey(meta.parent, meta.n_parts), [m, t, r, u]] }
            .groupTuple()
            .map { key, parts ->
                def sorted = parts.sort(false) { p -> p[0].name }
                [
                    key.getGroupTarget(),
                    sorted.collect { p -> p[0] },
                    sorted.collect { p -> p[1] },
                    sorted.collect { p -> p[2] },
                    opts.emit_uc ? sorted.collect { p -> p[3] } : []
                ]
            }
        MERGE_COMPONENTS(ch_merge)

        ch_membership      = MERGE_COMPONENTS.out.membership
        ch_representatives = MERGE_COMPONENTS.out.representatives
        ch_otu_taxonomy    = MERGE_COMPONENTS.out.otu_taxonomy
        ch_uc              = MERGE_COMPONENTS.out.uc
    }
    else {
        ch_membership      = OPTIMOTU_CLUSTER.out.membership
        ch_representatives = OPTIMOTU_CLUSTER.out.representatives
        ch_otu_taxonomy    = OPTIMOTU_CLUSTER.out.otu_taxonomy
        ch_uc              = OPTIMOTU_CLUSTER.out.uc
    }

    if (split_by_taxon) {
        def as_first_part = { meta, f -> [meta + [id: "${meta.id}.p000", parent: meta], f] }
        ch_sequence_taxonomy = OPTIMOTU_CLUSTER.out.sequence_taxonomy
            .mix(OPTIMOTU_RESOLVE.out.sequence_taxonomy.map(as_first_part))
        ch_run_records = OPTIMOTU_CLUSTER.out.run_record
            .mix(OPTIMOTU_RESOLVE.out.run_record.map(as_first_part))
    }
    else {
        ch_sequence_taxonomy = OPTIMOTU_CLUSTER.out.sequence_taxonomy
        ch_run_records       = OPTIMOTU_CLUSTER.out.run_record
    }

    emit:
    membership        = ch_membership
    representatives   = ch_representatives
    otu_taxonomy      = ch_otu_taxonomy
    uc                = ch_uc
    sequence_taxonomy = ch_sequence_taxonomy
    taxprob           = ch_taxprob
    run_records       = ch_run_records
}
