#!/usr/bin/env nextflow
/*

============================================================================
  NextITS: Pipeline to process fungal ITS amplicons
============================================================================
  License: Apache-2.0
  Github : https://github.com/vmikk/NextITS
  Website: https://next-its.github.io/
----------------------------------------------------------------------------
*/

// NB!!:
// - provide absolute paths to the input data (e.g. --input and --barcodes)
// - File names should not contain period (.) characters (except for extensions)

// Databases:
//  - UDB for chimera identification


// Enable DSL2 syntax
nextflow.enable.dsl = 2

// Enable topic channels
nextflow.preview.topic = true

// Include functions
include { helpMsg }                   from './modules/help_message.nf'


// Include workflows
// NB! `include` statements are static, meaning they are resolved at compile time rather than at runtime!
include { S1 } from './workflows/STEP1.nf'
include { S2 } from './workflows/STEP2.nf'
include { seqstats } from './workflows/STEP1.nf'


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  VALIDATE INPUTS

// Validate & print parameter summary
// NB! works only with old schema (`everit-json-schema` library doesn't support JSON Schema draft-2020-12)
WorkflowMain.initialise(workflow, params, log)

// Include the pipeline initialisation subworkflow
// requires newer nf-core template and schema
// include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_NextITS_pipeline'


// Show help msg
if (params.helpMsg){
    helpMsg()
    exit(0)
}


// Additional parameter validation for Step-1
if( params.step == "Step1" ) {

  if (params.input == false && params.seqplatform == "PacBio") {
      println( "ERROR: Please provide the input file with sequences in FASTQ.gz or BAM format with `--input` parameter.")
      exit(1)
  }
  if (params.input_R1 == false && params.input_R2 == false && params.seqplatform == "Illumina") {
      println( "ERROR: Please provide input files with sequences in FASTQ.gz format with `--input_R1` and `--input_R2` parameters.")
      exit(1)
  }
  if (params.barcodes == false && params.demultiplexed == false) {
      println( "ERROR: Please provide the file with sample barcodes in FASTA format with `--barcodes` parameter.")
      exit(1)
  }
  if (!params.chimera_db || !file(params.chimera_db).exists()) {
      println( "ERROR: Please provide the UDB file with reference sequences for chimera removal with `--chimera_db` parameter.")
      exit(1)
  }
  if (!(params.chimera_db.toLowerCase().endsWith('.udb'))) {
      println( "ERROR: The reference database file specified with `--chimera_db` parameter must be in UDB format." )
      exit 1
  }
  if (params.hp == true && params.seqplatform == "Illumina" && params.illumina_keep_notmerged == true) {
      println( "ERROR: Homopolymer compression is not implemented for Illumina non-merged reads.")
      exit(1)
  }
  if (params.seqplatform == "Illumina" && params.demultiplexed == true) {
      println( "ERROR: Handling demultiplexed data for Illumina is not implemented yet.")
      exit(1)
  }

  if (params.seqplatform == "Illumina" && params.illumina_keep_notmerged == true && params.its_region != "none") {
      println( "WARNING: Unmerged Illumina reads are not compatible with ITSx. Amplicons will be primer-trimmed.")
  }

}  // end of Step-1 parameter validation



// Print the parameters to the console and to the log
/*
log.info """
    =======================================================================
    NextITS v.${workflow.manifest.version}
    =======================================================================
    Input data path: ${params.input}
    Barcodes:        ${params.barcodes}
    Output path:     ${params.outdir}
    ITS region:      ${params.its_region}
    """
    .stripIndent()

log.info """
        Pipeline info:
          Pipeline profile:       ${workflow.profile}
          Config file used:       ${workflow.configFiles}
          Container engine:       ${workflow.containerEngine}
        """
        .stripIndent()

log.info """
        Core Nextflow options:
          launchDir:              ${workflow.launchDir}
          workDir:                ${workflow.workDir}
          projectDir:             ${workflow.projectDir}
        """
        .stripIndent()

log.info "======================================================================="
log.info "\n"
*/


workflow {

  // Primer disambiguation
  disambiguate()


  // Demultiplex data
  if( params.demultiplexed == false ){
    
    // Input file with barcodes (FASTA)
    ch_barcodes = Channel.value(params.barcodes)

    // Validate tags
    tag_validation(ch_barcodes)

    // PacBio
    if ( params.seqplatform == "PacBio" ) {
      
      // Input file with multiplexed reads (FASTQ.gz or BAM)
      ch_input = Channel.value(params.input)

      // Check the extension of input
      input_type = file(params.input).getExtension() =~ /bam|BAM/ ? "bam" : "oth"
      // println("${input_type}")

      // If BAM is provided as input, convert it to FASTQ
      if ( input_type == 'bam'){

        // Add BAM index file
        ch_input_pbi = ch_input + ".pbi"

        bam2fastq(ch_input, ch_input_pbi)
        qc_se(bam2fastq.out.fastq)

      } else {

        // Initial QC
        qc_se(ch_input)

      }

      // Demultiplexing with dual barcodes requires 3 additional files
      //  - "biosamples" with symmertic/asymmetirc tag combinations
      //  - and a table for assigning sample names to demuxed files
      // Create dummy files (for single or symmetic tags) if neccesary
      ch_biosamples_sym  = tag_validation.out.biosamples_sym.flatten().collect().ifEmpty(file("biosamples_sym"))
      ch_biosamples_asym = tag_validation.out.biosamples_asym.flatten().collect().ifEmpty(file("biosamples_asym"))
      ch_file_renaming   = tag_validation.out.file_renaming.flatten().collect().ifEmpty(file("file_renaming"))

      // Demultiplexing
      demux(
        qc_se.out.filtered,
        tag_validation.out.fasta,
        ch_biosamples_sym, 
        ch_biosamples_asym,
        ch_file_renaming)

      // Check primers
      primer_check(
        demux.out.samples_demux.flatten(),
        disambiguate.out.F,
        disambiguate.out.R,
        disambiguate.out.Fr,
        disambiguate.out.Rr
        )

    } // end of PacBio-specific tasks

   // Illumina 
    if ( params.seqplatform == "Illumina" ) {
      
      // Input file with multiplexed pair-end reads (FASTQ.gz)
      ch_inputR1 = Channel.value(params.input_R1)
      ch_inputR2 = Channel.value(params.input_R2)

      // Initial QC
      qc_pe(ch_inputR1, ch_inputR2)

      // PE assembly
      merge_pe(
        qc_pe.out.filtered_R1,
        qc_pe.out.filtered_R2)

      // Modify barcodes (restict search window)
      prep_barcodes(tag_validation.out.fasta)

      // Demultiplexing
      demux_illumina(
        merge_pe.out.r12,
        prep_barcodes.out.barcodesm)

      ch_demux_merged = demux_illumina.out.samples_demux.flatten()

      // Illumina nonmerged PE reads sub-workflow (optional)
      if(params.illumina_keep_notmerged == true){

        // Demultiplexing non-merged reads
        demux_illumina_notmerged(
          merge_pe.out.nm,
          prep_barcodes.out.barcodesm)

        // Channel of non-merged reads by sample (split into sample tuples)
        // ch_R1 = demux_illumina_notmerged.out.demux_pe....

        // Non-merged sample list
        ch_nonmerged = demux_illumina_notmerged.out.samples_nonm_pe.splitText().map{it -> it.trim()}

        // Trim primers of nonmerged PE reads
        // Estimate sequence qualities
        // Dereplicate R1 and R2 independently
        // trim_primers_pe(demux_illumina_notmerged.out.demux_pe.flatten())

        // Join nonmerged reads with poly-N pads
        join_pe(
          ch_nonmerged,
          demux_illumina_notmerged.out.demux_pe.flatten().collect()   // all non-merged R1 and R2 files
          )

        // Add joined reads to the merged reads
        ch_joined = join_pe.out.jj_FQ.flatten()
        ch_demuxed = ch_demux_merged.concat(ch_joined)

      } else { // end of Illumina non-merged reads

        // Channel with demultiplexed reads
        ch_demuxed = ch_demux_merged

      }

      // Check primers
      primer_check(
        ch_demuxed,
        disambiguate.out.F,
        disambiguate.out.R,
        disambiguate.out.Fr,
        disambiguate.out.Rr
        )

    } // end of Illumina-specific tasks

  }   // end of demultiplexing



  // If samples were already demuliplexed
  if( params.demultiplexed == true ){

    // Input files with demultiplexed reads (FASTQ.gz)
    ch_input = Channel.fromPath( params.input + '/*.{fastq.gz,fastq,fq.gz,fq}' )

    // QC
    qc_se(ch_input)

    // Check primers
    primer_check(
      qc_se.out.filtered,
      disambiguate.out.F,
      disambiguate.out.R,
      disambiguate.out.Fr,
      disambiguate.out.Rr
      )

  }  // end of pre-demultiplexed branch



    // Extract ITS
    if(params.its_region == "full" || params.its_region == "ITS1" || params.its_region == "ITS2" || params.its_region == "SSU" || params.its_region == "LSU"){

      // Run ITSx
      itsx(primer_check.out.fq_primer_checked)

      // Merge tables with sequence qualities
      seq_qual(itsx.out.hashes.collect())
    }

    // Trim the primers (instead of ITS extraction)
    if(params.its_region == "none"){
      
      // Trim primers with cutadapt
      trim_primers(primer_check.out.fq_primer_checked)

      // Merge tables with sequence qualities
      seq_qual(trim_primers.out.hashes.collect())
    }

    // Trim the primers, run ITSx, and assemble near-full-length ITS
    if(params.its_region == "ITS1_5.8S_ITS2"){

      // Run ITSx
      itsx(primer_check.out.fq_primer_checked)

      // Assemble ITS1-5.8S-ITS2 from ITSx-extracted parts
      if (params.ITSx_partial == 0) {
        assemble_its(
          itsx.out.itsx_its1,
          itsx.out.itsx_58s,
          itsx.out.itsx_its2)
      } else {
        assemble_its(
          itsx.out.itsx_its1_part,
          itsx.out.itsx_58s,
          itsx.out.itsx_its2_part)
      }

      // Merge tables with sequence qualities
      seq_qual(itsx.out.hashes.collect())
    }


    // Homopolymer compression
    if(params.hp == true){

        // --Full-length ITS sequences
        if(params.its_region == "full"){
          homopolymer(itsx.out.itsx_full)
        }
        // --ITS1 sequences
        if(params.its_region == "ITS1"){
          if (params.ITSx_partial == 0) {
            homopolymer(itsx.out.itsx_its1)
          } else {
            homopolymer(itsx.out.itsx_its1_part)
          }
        }
        // --ITS2 sequences
        if(params.its_region == "ITS2"){
          if (params.ITSx_partial == 0) {
            homopolymer(itsx.out.itsx_its2)
          } else {
            homopolymer(itsx.out.itsx_its2_part)
          }
        }
        // --SSU sequences
        if(params.its_region == "SSU"){
          if (params.ITSx_partial == 0) {
            homopolymer(itsx.out.itsx_ssu)
          } else {
            homopolymer(itsx.out.itsx_ssu_part)
          }
        }
        // --LSU sequences
        if(params.its_region == "LSU"){
          if (params.ITSx_partial == 0) {
            homopolymer(itsx.out.itsx_lsu)
          } else {
            homopolymer(itsx.out.itsx_lsu_part)
          }
        }

        // --Primer-trimmed sequences
        if(params.its_region == "none"){
          homopolymer(trim_primers.out.primertrimmed_fa)
        }
        // Near-full-length ITS
        if(params.its_region == "ITS1_5.8S_ITS2"){
          homopolymer(assemble_its.out.itsnf)
        }
    
        // Reference-based chimera removal
        ch_chimerabd = Channel.value(params.chimera_db)
        chimera_ref(homopolymer.out.hc, ch_chimerabd)
    
        // De novo chimera search
        chimera_denovo(homopolymer.out.hc)

    } else {
      // No homopolymer comression is required,
      // Just dereplicate the data

      if(params.its_region == "full" || params.its_region == "ITS1" || params.its_region == "ITS2" || params.its_region == "SSU" || params.its_region == "LSU"){
        
        // --Full-length ITS sequences
        if(params.its_region == "full"){
          just_derep(itsx.out.itsx_full)
        }
        // --ITS1 sequences
        if(params.its_region == "ITS1"){
          if (params.ITSx_partial == 0) {
            just_derep(itsx.out.itsx_its1)
          } else {
            just_derep(itsx.out.itsx_its1_part)
          }
        }
        // --ITS2 sequences
        if(params.its_region == "ITS2"){
          if (params.ITSx_partial == 0) {
            just_derep(itsx.out.itsx_its2)
          } else {
            just_derep(itsx.out.itsx_its2_part)
          }
        }
        // --SSU sequences
        if(params.its_region == "SSU"){
          if (params.ITSx_partial == 0) {
            just_derep(itsx.out.itsx_ssu)
          } else {
            just_derep(itsx.out.itsx_ssu_part)
          }
        }
        // --LSU sequences
        if(params.its_region == "LSU"){
          if (params.ITSx_partial == 0) {
            just_derep(itsx.out.itsx_lsu)
          } else {
            just_derep(itsx.out.itsx_lsu_part)
          }
        }

        // Reference-based chimera removal
        ch_chimerabd = Channel.value(params.chimera_db)
        chimera_ref(just_derep.out.nhc, ch_chimerabd)

        // De novo chimera search
        chimera_denovo(just_derep.out.nhc)

      }  // end of ITS

      // --Primer-trimmed sequences are already dereplicated
      if(params.its_region == "none"){
          
        // just_derep(trim_primers.out.primertrimmed_fa)

        // Reference-based chimera removal
        ch_chimerabd = Channel.value(params.chimera_db)
        chimera_ref(trim_primers.out.primertrimmed_fa, ch_chimerabd)

        // De novo chimera search
        chimera_denovo(trim_primers.out.primertrimmed_fa)
      }

      // Assembled ITS is also primer trimmed and dereplicated
      if(params.its_region == "ITS1_5.8S_ITS2"){

        // Reference-based chimera removal
        ch_chimerabd = Channel.value(params.chimera_db)
        chimera_ref(assemble_its.out.itsnf, ch_chimerabd)

        // De novo chimera search
        chimera_denovo(assemble_its.out.itsnf)
      }

    } // end of homopolymer correction condition

    // Chimera rescue
    ch_chimerafiles = chimera_ref.out.chimeric.collect()
    chimera_rescue(ch_chimerafiles)

    // Aggregate de novo chimeras into a single file
    chimera_denovo_agg(chimera_denovo.out.denovochim.collect())

    // Create channel with filtered reads
    ch_filteredseqs = chimera_ref.out.nonchimeric
      .concat(chimera_rescue.out.rescuedchimeric)
      .collect()

    // Global dereplication
    glob_derep(ch_filteredseqs)

    // Pool sequences (for a final sequence table)
    pool_seqs(ch_filteredseqs)

    // OTU clustering
    otu_clust(glob_derep.out.globderep)

    // Create OTU table
    otu_tab(
      otu_clust.out.otus,
      pool_seqs.out.seqsnf)

    // Tag-jump removal
    tj(otu_tab.out.otutab)


    // Check optional channel with de novo chimera scores
    ch_denovoscores = chimera_denovo_agg.out.alldenovochim.ifEmpty(file('DeNovo_Chimera.txt'))

    // Create sequence table
    prep_seqtab(
      pool_seqs.out.seqtabnf,               // non-filtered sequence table
      pool_seqs.out.seqsnf,                 // Sequences in FASTA format
      otu_tab.out.samples_uc,               // sequence mapping to OTUs
      tj.out.tjs,                           // tag-jumped OTU list
      ch_denovoscores,                      // de novo chimera scores
      seq_qual.out.quals                    // sequence qualities
      )


    
    //// Read count summary
    
      // Initial data - Per-sample input channels
      if( params.demultiplexed == false ){

        if(params.seqplatform == "PacBio"){
          ch_all_demux = demux.out.samples_demux.flatten().collect()
        }

        if(params.seqplatform == "Illumina"){
          ch_all_demux = demux_illumina.out.samples_demux.flatten().collect()
        }

      } else {
        ch_all_demux = Channel.fromPath( params.input + '/*.{fastq.gz,fastq,fq.gz,fq}' ).flatten().collect()
      }
        
      // Primer-checked and multiprimer sequences
      ch_all_primerchecked = primer_check.out.fq_primer_checked.flatten().collect().ifEmpty(file("no_primerchecked"))
      ch_all_primerartefacts = primer_check.out.primerartefacts.flatten().collect().ifEmpty(file("no_multiprimer"))
      
      // ITSx and primer trimming channel
      if(params.its_region == "full"){
        ch_all_trim = itsx.out.itsx_full.flatten().collect().ifEmpty(file("no_itsx"))
      }
      if(params.its_region == "ITS1"){
        if (params.ITSx_partial == 0) {
          ch_all_trim = itsx.out.itsx_its1.flatten().collect().ifEmpty(file("no_itsx"))
        } else {
          ch_all_trim = itsx.out.itsx_its1_part.flatten().collect().ifEmpty(file("no_itsx"))
        }
      }
      if(params.its_region == "ITS2"){
        if (params.ITSx_partial == 0) {
          ch_all_trim = itsx.out.itsx_its2.flatten().collect().ifEmpty(file("no_itsx"))
        } else {
          ch_all_trim = itsx.out.itsx_its2_part.flatten().collect().ifEmpty(file("no_itsx"))
        }
      }
      if(params.its_region == "SSU"){
        if (params.ITSx_partial == 0) {
          ch_all_trim = itsx.out.itsx_ssu.flatten().collect().ifEmpty(file("no_itsx"))
        } else {
          ch_all_trim = itsx.out.itsx_ssu_part.flatten().collect().ifEmpty(file("no_itsx"))
        }
      }
      if(params.its_region == "LSU"){
        if (params.ITSx_partial == 0) {
          ch_all_trim = itsx.out.itsx_lsu.flatten().collect().ifEmpty(file("no_itsx"))
        } else {
          ch_all_trim = itsx.out.itsx_lsu_part.flatten().collect().ifEmpty(file("no_itsx"))
        }
      }
      if(params.its_region == "ITS1_5.8S_ITS2"){
        ch_all_trim = assemble_its.out.itsnf.flatten().collect().ifEmpty(file("no_itsx"))
      }
      if(params.its_region == "none"){
        ch_all_trim = trim_primers.out.primertrimmed_fq.flatten().collect().ifEmpty(file("no_primertrim"))
      }

      // Homopolymer-correction channel
      if(params.hp == true){
        ch_homopolymers = homopolymer.out.uch.flatten().collect().ifEmpty(file("no_homopolymer"))
      } else {
        ch_homopolymers = file("no_homopolymer")
      }

      // Chimeric channels
      ch_chimref     = chimera_ref.out.chimeric.flatten().collect().ifEmpty(file("no_chimref"))
      ch_chimdenovo  = chimera_denovo.out.denovochim.flatten().collect().ifEmpty(file("no_chimdenovo"))
      ch_chimrescued = chimera_rescue.out.rescuedchimeric.flatten().collect().ifEmpty(file("no_chimrescued"))

      // Tag-jump filtering channel
      ch_tj = tj.out.tjs


      // Count reads and prepare summary stats for the run
      // Currently, implemented only for PacBio
      // For Illumina, need replace:
      //   `ch_input` -> `ch_inputR1` & `ch_inputR2`
      //   `qc_se`    -> `qc_pe`

      if(params.seqplatform == "PacBio"){

      read_counts(
          ch_input,                // input data
          qc_se.out.filtered,      // data that passed QC
          ch_all_demux,            // demultiplexed sequences per sample
          ch_all_primerchecked,    // primer-cheched sequences
          ch_all_primerartefacts,  // multiprimer artefacts
          ch_all_trim,             // ITSx-extracted or primer-trimmed sequences
          ch_homopolymers,         // Homopolymer stats
          ch_chimref,              // Reference-based chimeras
          ch_chimdenovo,           // De novo chimeras
          ch_chimrescued,          // Rescued chimeras
          ch_tj,                   // Tag-jump filtering stats
          prep_seqtab.out.seq_rd   // Final table with sequences
          )

      } // end of read_counts for PacBio


  // Collect ITSx-extracted sequences
  if(params.its_region == "full" || params.its_region == "ITS1" || params.its_region == "ITS2" || params.its_region == "SSU" || params.its_region == "LSU" || params.its_region == "ITS1_5.8S_ITS2"){

    // Collect rRNA parts into separate channels
    ch_cc_full = itsx.out.itsx_full.flatten().collect().ifEmpty(file("NOFULL"))
    ch_cc_ssu  = itsx.out.itsx_ssu.flatten().collect().ifEmpty(file("NOSSU"))
    ch_cc_its1 = itsx.out.itsx_its1.flatten().collect().ifEmpty(file("NOITS1"))
    ch_cc_58s  = itsx.out.itsx_58s.flatten().collect().ifEmpty(file("NO58S"))
    ch_cc_its2 = itsx.out.itsx_its2.flatten().collect().ifEmpty(file("NOITS2"))
    ch_cc_lsu  = itsx.out.itsx_lsu.flatten().collect().ifEmpty(file("NOLSU"))
    
    ch_cc_ssu_part  = itsx.out.itsx_ssu_part.flatten().collect().ifEmpty(file("NOSSUPART"))
    ch_cc_its1_part = itsx.out.itsx_its1_part.flatten().collect().ifEmpty(file("NOITS1PART"))
    ch_cc_58s_part  = itsx.out.itsx_58s_part.flatten().collect().ifEmpty(file("NO58SPART"))
    ch_cc_its2_part = itsx.out.itsx_its2_part.flatten().collect().ifEmpty(file("NOITS2PART"))
    ch_cc_lsu_part  = itsx.out.itsx_lsu_part.flatten().collect().ifEmpty(file("NOLSUPART"))

    itsx_collect(
      ch_cc_full,
      ch_cc_ssu,
      ch_cc_its1,
      ch_cc_58s,
      ch_cc_its2,
      ch_cc_lsu,
      ch_cc_ssu_part,
      ch_cc_its1_part,
      ch_cc_58s_part,
      ch_cc_its2_part,
      ch_cc_lsu_part
      )

  }
  
  
  // Dump the software versions to a file
  software_versions_to_yaml(Channel.topic('versions'))
      .collectFile(
          storeDir: "${params.tracedir}",
          name:     'software_versions.yml',
          sort:     true,
          newLine:  true
      )

}


// Quick workflow for demultiplexing and estimation of the number of reads per sample
// Only PacBio non-demultiplexed reads are supported
workflow seqstats {

  // Primer disambiguation
  disambiguate()

  // Input file with barcodes (FASTA)
  ch_barcodes = Channel.value(params.barcodes)

  // Input file with multiplexed reads (FASTQ.gz)
  ch_input = Channel.value(params.input)

  // Initial QC
  qc_se(ch_input)

  // Validate tags
  tag_validation(ch_barcodes)

  // Tag-validation channels
  ch_biosamples_sym  = tag_validation.out.biosamples_sym.flatten().collect().ifEmpty(file("biosamples_sym"))
  ch_biosamples_asym = tag_validation.out.biosamples_asym.flatten().collect().ifEmpty(file("biosamples_asym"))
  ch_file_renaming   = tag_validation.out.file_renaming.flatten().collect().ifEmpty(file("file_renaming"))

  // Demultiplexing
  demux(
    qc_se.out.filtered,
    tag_validation.out.fasta,
    ch_biosamples_sym, 
    ch_biosamples_asym,
    ch_file_renaming)

  // Check primers
  primer_check(
    demux.out.samples_demux.flatten(),
    disambiguate.out.F,
    disambiguate.out.R,
    disambiguate.out.Fr,
    disambiguate.out.Rr
    )

  // Prepare input channels
  ch_all_demux = demux.out.samples_demux.flatten().collect()
  ch_all_primerchecked = primer_check.out.fq_primer_checked.flatten().collect().ifEmpty(file("no_primerchecked"))
  ch_all_primerartefacts = primer_check.out.primerartefacts.flatten().collect().ifEmpty(file("no_multiprimer"))

  // Count reads and prepare summary stats for the run
  quick_stats(
      ch_input,                // input data
      qc_se.out.filtered,      // data that passed QC
      ch_all_demux,            // demultiplexed sequences per sample
      ch_all_primerchecked,    // primer-cheched sequences
      ch_all_primerartefacts   // primer artefacts
      )

} // end of `seqstats` subworkflow




// On completion
workflow.onComplete {
    println "Pipeline completed at : $workflow.complete"
    println "Duration              : ${workflow.duration}"
    println "Execution status      : ${workflow.success ? 'All done!' : 'Failed' }"
}

// On error
workflow.onError {
    println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
