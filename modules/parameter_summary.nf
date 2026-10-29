/*
 * Custom parameter summary function for NextITS pipeline
 * Shows only step-specific parameters and relevant configuration
 */

include { resolveOptimotuRefset } from './optimotu_refset'
include { getColors } from './colors'

def parameterSummary(workflow, params) {
    
    // Color definitions
    def colors = getColors(params.monochrome_logs)
    
    def summary = ""
    summary += "${colors.bold}Parameters:${colors.reset}\n"
    summary += "  ${colors.green}step${colors.reset}             : ${colors.cyan}${params.step}${colors.reset}\n"
    
    // Step-specific parameters
    if (params.step == "Step1" || params.step == "seqstats") {
        
        if (params.seqplatform == "PacBio" || params.demultiplexed) {
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
        if (params.seqplatform == "Illumina") {
            summary += "  ${colors.green}seqplatform${colors.reset}      : ${colors.cyan}Illumina${colors.reset}\n"
            if (!params.demultiplexed) {
                summary += "  ${colors.green}tag layout${colors.reset}       : ${colors.cyan}${params.illumina_barcodetype}${colors.reset}\n"
            }
            summary += "  ${colors.green}pe_merger${colors.reset}        : ${colors.cyan}${params.pe_merger}${colors.reset}\n"
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
            if (params.long_chimera) {
                summary += "  ${colors.green}long_chimera${colors.reset}     : ${colors.cyan}true (minbits ${params.long_chimera_minbits}, mincov ${params.long_chimera_mincov})${colors.reset}\n"
            }
        }

    } 
    else if (params.step == "Step2") {
        summary += "\n${colors.bold}Step 2 inputs:${colors.reset}\n"
        summary += "  ${colors.green}data_path${colors.reset}        : ${colors.cyan}${params.data_path}${colors.reset}\n"
        summary += "  ${colors.green}clustering${colors.reset}       : ${colors.cyan}${params.clustering}${colors.reset}\n"
        summary += "  ${colors.green}preclustering${colors.reset}    : ${colors.cyan}${params.preclustering}${colors.reset}\n"
        if (params.clustering == "optimotu") {
            def refset = resolveOptimotuRefset(params)
            def refset_txt = refset.error ? "${params.optimotu_db} (not usable)" : "${refset.id} (${refset.dir})"
            summary += "  ${colors.green}optimotu_db${colors.reset}      : ${colors.cyan}${refset_txt}${colors.reset}\n"
            summary += "  ${colors.green}split_rank${colors.reset}       : ${colors.cyan}${params.optimotu_split_rank}${colors.reset}\n"
            summary += "  ${colors.green}closed_ref_from${colors.reset}  : ${colors.cyan}${params.optimotu_closed_ref_from}${colors.reset}\n"
        }
    }

    
    // Common parameters
    summary += "\n${colors.bold}Output:${colors.reset}\n"
    summary += "  ${colors.green}outdir${colors.reset}           : ${colors.cyan}${params.outdir}${colors.reset}\n"
    summary += "  ${colors.green}workDir${colors.reset}          : ${colors.cyan}${workflow.workDir}${colors.reset}\n"
    
    // Nextflow configuration
    summary += "\n${colors.bold}Config:${colors.reset}\n"
    summary += "  ${colors.green}NextITS version${colors.reset}  : ${colors.cyan}${workflow.manifest.version}${colors.reset}\n"
    if(workflow.commitId) {
        summary += "  ${colors.green}NextITS revision${colors.reset} : ${colors.cyan}${workflow.commitId.substring(0, 7)}${colors.reset}\n"
    }
    summary += "  ${colors.green}Profile${colors.reset}          : ${colors.cyan}${workflow.profile}${colors.reset}\n"
    summary += "  ${colors.green}Container engine${colors.reset} : ${colors.cyan}${workflow.containerEngine ?: 'none'}${colors.reset}\n"
    
    if (workflow.container) {
        def container_txt = workflow.container
        if (workflow.container instanceof Map) {
            // Main container first; OptimOTU container only when it is used (Step-2 with `--clustering optimotu`)
            def show_optimotu = params.step == "Step2" && params.clustering == "optimotu"
            def entries = workflow.container.findAll { k, v -> show_optimotu || !k.toString().contains('optimotu') }
            def main_ct = entries.remove('default')
            def other   = entries.collect { k, v -> "${v} (${k})" }
            container_txt = ([main_ct] + other).findAll { it }.join(', ')
        }
        summary += "  ${colors.green}Container${colors.reset}        : ${colors.cyan}${container_txt}${colors.reset}\n"
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
