/*
 * Parameter helpers for values that may arrive as strings under the Nextflow strict syntax parser
 */

def parseBooleanParam(value, name) {
    if (value instanceof Boolean) {
        return value
    }

    if (value == null) {
        return false
    }

    def normalized = value.toString().trim().toLowerCase()

    if (['true', 't', '1', 'yes', 'y'].contains(normalized)) {
        return true
    }

    if (['false', 'f', '0', 'no', 'n'].contains(normalized)) {
        return false
    }

    error("Parameter --${name} must be a boolean value: true/false, yes/no, or 1/0. Received: ${value}")
}
