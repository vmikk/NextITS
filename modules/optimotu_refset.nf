/*
 * OptimOTU reference sets (used with `--clustering optimotu`)
 *
 * A reference set is a directory `<db_dir>/optimotu/<id>/` with:
 *   manifest.json                   - id, title, region, classifier, file names, citation
 *   thresholds.tsv, thresholds.json - taxon-specific thresholds and the settings they were trained with
 *   classifier.sintax.fasta.gz      - classifier database (SINTAX-format FASTA)
 *   classifier.stx                  - optional prebuilt sintaxer index (manifest `classifier_index`). Used instead of indexing the FASTA in every run
 *
 * Thresholds are valid only with the taxonomy, region, and classifier they were trained with, so these are bundled together and checked here
 */

include { referenceSetFiles } from '../optimotu-nf/subworkflows/optimotu_options'

// Distance penalties of the default (edit) distance; other values are passed to OptimOTU explicitly
def optimotuDefaultPenalties() {
    return [match: 0, mismatch: 1, gap_open: 0, gap_extend: 1, gap_open2: 0, gap_extend2: 1]
}

/*
 * Locate and check the reference set given with `--optimotu_db` (an id or a directory path)
 * Returns [id, dir, manifest, thresholds, classifier_db, options],
 * or [error: message] if the set can not be used
 */
def resolveOptimotuRefset(params) {

    // Reference-set directory
    def dir = file(params.optimotu_db)
    if (!dir.isDirectory()) {
        def base = params.db_dir ?: System.getenv('NEXTITS_DB_DIR')
        if (!base) {
            return [error: "Cannot locate OptimOTU reference set `${params.optimotu_db}`: set `--db_dir` (or the NEXTITS_DB_DIR environment variable) to the folder with reference databases, or give the path to the reference-set folder with `--optimotu_db`."]
        }
        dir = file("${base}/optimotu/${params.optimotu_db}")
        if (!dir.isDirectory()) {
            def available = files("${base}/optimotu/*", type: 'dir').collect { d -> d.name }.sort()
            return [error: "OptimOTU reference set `${params.optimotu_db}` not found in ${base}/optimotu/ (available: ${available ? available.join(', ') : 'none'})."]
        }
    }

    // Manifest and the thresholds sidecar
    def manifest_file = dir.resolve('manifest.json')
    if (!manifest_file.exists()) {
        return [error: "OptimOTU reference set has no manifest.json: ${dir}"]
    }
    def manifest = new groovy.json.JsonSlurper().parseText(manifest_file.text)
    if (manifest.classifier != 'sintaxer') {
        return [error: "OptimOTU reference set ${dir} uses classifier `${manifest.classifier}`; only `sintaxer` is supported."]
    }
    if (manifest.region != 'ITS_full') {
        return [error: "OptimOTU reference set ${dir} is for region `${manifest.region}`; only full-length ITS (`ITS_full`) is supported."]
    }
    def thresholds    = dir.resolve(manifest.thresholds ?: 'thresholds.tsv')
    def sidecar_file  = file(thresholds.toString().replaceFirst(/\.tsv$/, '.json'))
    def classifier_db = dir.resolve(manifest.classifier_db ?: 'classifier.sintax.fasta.gz')
    def missing = [thresholds, sidecar_file, classifier_db].findAll { f -> !f.exists() }
    if (missing) {
        return [error: "OptimOTU reference set ${dir} is incomplete, missing: ${missing.join(', ')}"]
    }
    def sidecar = new groovy.json.JsonSlurper().parseText(sidecar_file.text)

    // Newer manifests repeat the ranks the thresholds were trained with
    if (manifest.ranks && sidecar.ranks && manifest.ranks != sidecar.ranks) {
        return [error: "OptimOTU reference set ${dir}: ranks in manifest.json (${manifest.ranks.join(',')}) differ from thresholds.json (${sidecar.ranks.join(',')})."]
    }

    // The prebuilt sintaxer index when the set has one, otherwise the FASTA
    def classifier_input = referenceSetFiles(dir, 'sintaxer').classifier_db

    def measures = sidecar.optimization?.measures
    if (measures && !(params.optimotu_metric in measures)) {
        return [error: "OptimOTU reference set ${dir} has no thresholds for metric `${params.optimotu_metric}` (available: ${measures.join(', ')})."]
    }

    // OptimOTU options that must match the training settings
    def dist = sidecar.distance ?: [:]
    def penalties = dist.findAll { k, _v -> k != 'method' }
    def defaults = optimotuDefaultPenalties()
    def custom_penalties = penalties.any { k, v -> !defaults.containsKey(k) || (v as Integer) != defaults[k] }
    def options = [
        ranks:          sidecar.ranks ? sidecar.ranks.join(',') : null,
        root_taxon:     sidecar.root_taxon,
        rank_letters:   sidecar.rank_letters ? sidecar.rank_letters.collect { k, v -> "${k}=${v}" }.join(',') : null,
        dist_method:    dist.method,
        dist_penalties: custom_penalties ? penalties.collect { k, v -> "${k}=${v}" }.join(',') : null,
        clust_method:   sidecar.clustering?.method
    ].findAll { _k, v -> v != null }

    return [
        id:            manifest.id ?: dir.name,
        dir:           dir.toString(),
        manifest:      manifest,
        thresholds:    thresholds.toString(),
        classifier_db: classifier_input.toString(),
        options:       options
    ]
}

// Same as resolveOptimotuRefset, but stops the run if the set can not be used
def optimotuRefset(params) {
    def refset = resolveOptimotuRefset(params)
    if (refset.error) {
        error(refset.error)
    }
    return refset
}
