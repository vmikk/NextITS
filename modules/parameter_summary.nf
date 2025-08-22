/*
 * Custom parameter summary function for NextITS pipeline
 * Shows only step-specific parameters and relevant configuration
 */

def parameterSummary(workflow, params) {
    
    // Color definitions
    def colors = params.monochrome_logs ? [:] : [
        green: "\033[0;32m", blue: "\033[0;34m", yellow: "\033[0;33m", 
        cyan: "\033[0;36m", purple: "\033[0;35m", dim: "\033[2m", 
        bold: "\033[1m", reset: "\033[0m"
    ]
    
    def summary = ""
    summary += "${colors.bold}Parameters:${colors.reset}\n"
    summary += "  ${colors.green}step${colors.reset}             : ${colors.cyan}${params.step}${colors.reset}\n"
    
    // Step-specific parameters
    if (params.step == "Step1" || params.step == "seqstats") {
        
        if (params.seqplatform == "PacBio") {
            if (params.input) {
                summary += "  ${colors.green}input${colors.reset}            : ${colors.cyan}${params.input}${colors.reset}\n"
            }
        } else if (params.seqplatform == "Illumina") {
            if (params.input_R1) {
                summary += "  ${colors.green}input_R1${colors.reset}         : ${colors.cyan}${params.input_R1}${colors.reset}\n"
            }
            if (params.input_R2) {
                summary += "  ${colors.green}input_R2${colors.reset}         : ${colors.cyan}${params.input_R2}${colors.reset}\n"
            }
        }
        
        if (params.barcodes && !params.demultiplexed) {
            summary += "  ${colors.green}barcodes${colors.reset}         : ${colors.cyan}${params.barcodes}${colors.reset}\n"
        }
        
        if (params.demultiplexed) {
            summary += "  ${colors.green}demultiplexed${colors.reset}    : ${colors.cyan}${params.demultiplexed}${colors.reset}\n"
        }
        
        summary += "  ${colors.green}chimera_db${colors.reset}       : ${colors.cyan}${params.chimera_db}${colors.reset}\n"
        if(params.step == "Step1") {
            summary += "  ${colors.green}its_region${colors.reset}       : ${colors.cyan}${params.its_region}${colors.reset}\n"
        }

    } 
    else if (params.step == "Step2") {
        summary += "\n${colors.bold}Step 2 inputs:${colors.reset}\n"
        summary += "  ${colors.green}data_path${colors.reset}        : ${colors.cyan}${params.data_path}${colors.reset}\n"
        summary += "  ${colors.green}clustering${colors.reset}       : ${colors.cyan}${params.clustering}${colors.reset}\n"
        summary += "  ${colors.green}preclustering${colors.reset}    : ${colors.cyan}${params.preclustering}${colors.reset}\n"
    }

    
    // Common parameters
    summary += "\n${colors.bold}Output:${colors.reset}\n"
    summary += "  ${colors.green}outdir${colors.reset}           : ${colors.cyan}${params.outdir}${colors.reset}\n"
    summary += "  ${colors.green}workDir${colors.reset}          : ${colors.cyan}${workflow.workDir}${colors.reset}\n"
    
    // Nextflow configuration
    summary += "\n${colors.bold}Config:${colors.reset}\n"
    summary += "  ${colors.green}NextITS version${colors.reset}  : ${colors.cyan}${workflow.manifest.version}${colors.reset}\n"
    if(workflow.commitId) {
        summary += "  ${colors.green}NextITS revision${colors.reset}        : ${colors.cyan}${workflow.commitId.substring(0, 7)}${colors.reset}\n"
    }
    summary += "  ${colors.green}Profile${colors.reset}          : ${colors.cyan}${workflow.profile}${colors.reset}\n"
    summary += "  ${colors.green}Container engine${colors.reset} : ${colors.cyan}${workflow.containerEngine ?: 'none'}${colors.reset}\n"
    
    if (workflow.container) {
        summary += "  ${colors.green}Container${colors.reset}        : ${colors.cyan}${workflow.container}${colors.reset}\n"
    }
    
    return summary
}

// Export the function so it can be used in other files
workflow paramSummary {
    take:
    wf
    prms
    
    main:
    def summary = parameterSummary(wf, prms)
    log.info summary
}
