
// Custom function to parse software versions and return a YAML string

// Short commit the pipeline was run from
// `workflow.commitId` is only set when Nextflow pulled the pipeline itself (`nextflow run vmikk/NextITS`)
// for a local checkout, ask git directly
def nextits_revision() {
    if (workflow.commitId) {
        return workflow.commitId.substring(0, 7)
    }
    try {
        def proc = ["git", "-C", "${workflow.projectDir}", "rev-parse", "--short=7", "HEAD"].execute()
        proc.waitForOrKill(5000)
        if (proc.exitValue() == 0) {
            def rev = proc.text.trim()
            if (rev) {
                def dirty = ["git", "-C", "${workflow.projectDir}", "status", "--porcelain"].execute()
                dirty.waitForOrKill(5000)
                return (dirty.exitValue() == 0 && dirty.text.trim()) ? "${rev}-dirty" : rev
            }
        }
    } catch (Exception e) {
        // Not a git checkout, or no git available - the revision is simply omitted
    }
    return null
}

def software_versions_to_yaml(versions) {

    def revision = nextits_revision()

    def workflow_info = channel.of(
        "NextITS:\n" +
        "    version: ${workflow.manifest.version}\n" +
        (revision ? "    revision: ${revision}\n" : "") +
        "\nNextflow:\n" +
        "    version: ${nextflow.version}\n"
    )

    return workflow_info.mix(
        versions
            .unique()
            .map { name, tool, version ->
                [ name.tokenize(':')[-1], [ tool, version ] ]
            }
            .groupTuple()
            .map { processName, toolInfo ->
                def toolVersions = toolInfo.collect { tool, version -> "    ${tool}: ${version}" }.join('\n')
                "${processName}:\n${toolVersions}\n"
            }
            .map { it.trim() }
    )
}

