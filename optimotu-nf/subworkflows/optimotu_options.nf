/*
 * Options for the OPTIMOTU and OPTIMOTU_TRAIN subworkflows.
 *
 * The subworkflows and their processes never read `params`: a pipeline that
 * includes them has its own parameter namespace, and Nextflow loads `lib/` only
 * from the launched project. Options are passed in explicitly as one map, merged
 * over the defaults below by optimotuOptions().
 *
 * Standalone, main.nf builds that map from params (nextflow.config holds the
 * CLI defaults, which should match these). A host pipeline passes only what it changes, e.g.
 *
 *   OPTIMOTU(ch_fasta, thresholds, db, null, [confidence: params.optimotu_confidence])
 *
 * referenceSetFiles() resolves a reference set directory into those thresholds and db arguments
 */

def optimotuDefaults() {
    return [
        // --- shared --------------------------------------------------------
        ranks                     : 'root,kingdom,phylum,class,order,family,genus,species',
        root_taxon                : 'Eukaryota',
        rank_letters              : 'd=kingdom,k=kingdom,p=phylum,c=class,o=order,f=family,g=genus,s=species',
        taxonomy_format           : 'auto',
        dist_method               : 'wfa2',
        dist_penalties            : null,
        clust_method              : 'tree',

        // --- clustering ----------------------------------------------------
        classifier                : 'sintaxer',
        confidence                : 0.5,
        min_identity              : 0.8,
        classifier_seed           : 42,
        keep_placeholders         : false,
        force_denovo              : 'species',
        metric                    : 'FM',
        dereplicate               : false,
        emit_uc                   : true,
        split_components          : false,
        component_max_dist        : 0.4,
        split_rank                : 'none',
        closed_ref_from           : null,
        split_target              : 20000,

        // --- threshold training --------------------------------------------
        train_restrict            : null,
        train_require_rank        : null,
        train_min_per_taxon       : 1,
        train_max_per_taxon       : 100000000,
        train_max_seqs            : 100000000,
        train_seed                : 42,
        train_emit_sintax         : false,
        train_keep_orphan_ranks   : false,
        train_partition           : true,
        train_chunk_budget_pairs  : 25000000,
        train_chunk_subsample_max : 6000,
        train_chunk_min_per_taxon : 2,
        train_keep_underpowered   : false,
        dist_max                  : 0.4,
        dist_step                 : 0.001,
        min_taxa                  : 5,
        min_refseq                : 10,
        train_metric              : 'all',

        // --- reference set (thresholds + matching classifier database) -----
        train_reference_set       : true,
        train_reference_set_index : true,
        train_classifier_reference: null,   // null: the training reference
        label                     : null,   // null: the training reference file name stem
        refset_id                 : null,   // null: the label
        refset_title              : null,
        refset_database           : null,
        refset_database_version   : null,
        refset_region             : null,
        refset_citation           : null
    ]
}

/*
 * Boolean option values.
 *
 * With the strict parser, `--flag false` reaches the pipeline as the string
 * "false", which Groovy treats as true, so every boolean option goes through this.
 */
def asBool(value) {
    if (value instanceof Boolean) {
        return value
    }
    if (value == null) {
        return false
    }
    def text = value.toString().trim().toLowerCase()
    if (text in ['true', 'yes', '1', '']) {
        return true
    }
    if (text in ['false', 'no', '0', 'null']) {
        return false
    }
    error("Expected true or false, got '${value}'")
}

/*
 * Merge option overrides over optimotuDefaults().
 *
 * Unknown keys are an error, so a misspelled option fails here instead of being
 * ignored. Boolean options are normalized with asBool().
 */
def optimotuOptions(overrides) {
    def defaults = optimotuDefaults()
    def given = overrides ?: [:]
    def unknown = given.keySet().findAll { key -> !defaults.containsKey(key) }
    if (unknown) {
        error("Unknown OptimOTU option(s): ${unknown.join(', ')}. " +
              "Known options: ${defaults.keySet().join(', ')}")
    }
    def opts = defaults + given
    def flags = defaults.findAll { _key, value -> value instanceof Boolean }.keySet()
    flags.each { key -> opts[key] = asBool(opts[key]) }

    def classifiers = ['sintaxer', 'sintax', 'vsearch', 'none']
    if (opts.classifier !in classifiers) {
        error("Unknown classifier '${opts.classifier}'. Use one of: ${classifiers.join(', ')}")
    }

    // Closed-reference attachment from one rank down (null: at every rank)
    if (opts.closed_ref_from && opts.closed_ref_from.toString() != 'none') {
        def ranks = opts.ranks.toString().tokenize(',')*.trim()
        if (opts.closed_ref_from.toString() !in ranks.drop(1)) {
            error("Unknown closed_ref_from rank '${opts.closed_ref_from}'. Use one of: ${ranks.drop(1).join(', ')}")
        }
    }
    else {
        opts.closed_ref_from = null
    }

    // A taxonomic split stops after one rank strictly between the first and
    // the last, and it is not combined with a component split
    opts.split_rank = (opts.split_rank ?: 'none').toString()
    if (opts.split_rank != 'none') {
        def ranks = opts.ranks.toString().tokenize(',')*.trim()
        def allowed = ranks.size() > 2 ? ranks[1..-2] : []
        if (opts.split_rank !in allowed) {
            error("Unknown split rank '${opts.split_rank}'. Use 'none' or one of: ${allowed.join(', ')}")
        }
        if (opts.split_components) {
            error("split_rank and split_components cannot be combined")
        }
    }
    return opts
}

/*
 * Resolve the files of a reference set directory (written by --mode train,
 * see modules/reference_set.nf) for OPTIMOTU's thresholds and classifier_db.
 *
 * The thresholds and the classifier database are named in manifest.json. For
 * sintaxer, a prebuilt index is used when there is one: the manifest's
 * classifier_index, else <stem>.stx beside the database. Returns
 * [thresholds: Path, classifier_db: Path, manifest: Map]; the manifest's ranks
 * and root_taxon (absent in older manifests) are for the caller to compare
 * with its options.
 */
def referenceSetFiles(dir, classifier) {
    def base = file(dir.toString(), checkIfExists: true)
    def manifest_file = base.resolve('manifest.json')
    if (!manifest_file.exists()) {
        error("Reference set ${base} has no manifest.json")
    }
    def manifest = new groovy.json.JsonSlurper().parseText(manifest_file.text)
    if (!manifest.thresholds) {
        error("Reference set manifest ${manifest_file} names no thresholds")
    }
    def thresholds = base.resolve(manifest.thresholds.toString())
    if (!thresholds.exists()) {
        error("Reference set ${base} lacks its thresholds table ${manifest.thresholds}")
    }

    def db = manifest.classifier_db ? base.resolve(manifest.classifier_db.toString()) : null
    if (classifier == 'sintaxer') {
        // classifier.sintax.fasta.gz -> classifier.sintax.stx (SINTAXER_INDEX's
        // default name) or classifier.stx
        def stem = manifest.classifier_db
            ? manifest.classifier_db.toString().replaceFirst(/\.(fa|fasta|fas|fna)(\.gz|\.zst)?$/, '')
            : null
        def candidates = [
            manifest.classifier_index,
            manifest.files?.find { f -> f.role == 'classifier_index' }?.name,
            stem ? stem + '.stx' : null,
            stem ? stem.replaceFirst(/\.sintax$/, '') + '.stx' : null
        ].findAll { name -> name }
        def index = candidates.collect { name -> base.resolve(name.toString()) }.find { p -> p.exists() }
        if (index) {
            db = index
        }
    }
    if (classifier != 'none' && !(db && db.exists())) {
        error("Reference set ${base} lacks its classifier database ${manifest.classifier_db}")
    }
    return [thresholds: thresholds, classifier_db: db, manifest: manifest]
}
