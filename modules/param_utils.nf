/*
 * Helpers for parameters that can reach included scripts as strings when Nextflow uses the strict syntax parser
 */

def chunkingEnabled(value) {
    if (value == null) {
        return false
    }

    def normalized = value.toString().trim()

    if (!normalized.isInteger()) {
        error("Parameter --chunking_n must be an integer. Received: ${value}")
    }

    return normalized.toInteger() >= 2
}

// Poly-G trimming (`--qc_polyglen`), 0 or null = no trimming
// (in output declarations, a value from the command line is still a string, e.g. "0")
def polyGEnabled(value) {
    if (value == null) {
        return false
    }

    def normalized = value.toString().trim()

    if (!normalized.isInteger()) {
        error("Parameter --qc_polyglen must be an integer. Received: ${value}")
    }

    return normalized.toInteger() > 0
}

/*
 * Step-specific parameter schema
 * Parameters of the other step (or of the other sequencing platform-specific group) are kept in the schema
 * (so that their names are still recognized), but their values are not validated
 * Returns the path to a temporary schema file
 */
def stepSchema(schema_path, step) {

    def other_groups = [
        "Step1":    ["step_2_specific_parameters", "optimotu_clustering"],
        "seqstats": ["step_2_specific_parameters", "optimotu_clustering"],
        "Step2":    ["step_1_specific_parameters", "illumina_specific_parameters"]
    ][step] ?: []

    def schema = new groovy.json.JsonSlurper().parse(new File(schema_path.toString()))
    other_groups.each { grp ->
        def defn = schema['$defs'][grp]
        if (defn == null) {
            error("Parameter group `${grp}` not found in ${schema_path}")
        }
        defn.properties = defn.properties.collectEntries { k, _v -> [ (k): [:] ] }
        defn.remove('required')
        defn.remove('dependentRequired')
    }

    def tmp = java.nio.file.Files.createTempFile("NextITS_schema_${step}_", ".json").toFile()
    tmp.deleteOnExit()
    tmp.text = groovy.json.JsonOutput.toJson(schema)
    return tmp.absolutePath
}
