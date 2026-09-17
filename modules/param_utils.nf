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
