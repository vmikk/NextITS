
// Custom function to parse software versions and return a YAML string

def software_versions_to_yaml(versions) {
    
    def workflow_info = Channel.of(
        "NextITS:\n" +
        "    version: ${workflow.manifest.version}\n" +
        (workflow.commitId ? "    revision: ${workflow.commitId.substring(0,7)}\n" : "") +
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

