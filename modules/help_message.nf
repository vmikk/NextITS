

// Pipeline help message
def helpMsg() {
    log.info"""
    =====================================================================
    NextITS v.${workflow.manifest.version}
    =====================================================================
    
    Pipeline Usage:
    To run the pipeline, enter the following in the command line:
        nextflow run vmikk/nextits -r ${workflow.manifest.version} --input ... --outdir ...
    
    Options:
    REQUIRED:
        --input           File with single-end input sequences, PacBio (FASTQ or BAM) or a directory with pre-demultiplexed files
        --input_R1        Files with paired-end input sequences, Illumina (FASTQ)
        --input_R2
        --barcodes        Barcodes for demultiplexing (FASTA)
        --outdir          The output directory where the results will be saved

    OPTIONAL:
        --demultiplexed     Boolean, input is multiplexed (true, single FASTQ file) or pre-demultiplexed (multiple FASTQ files)
        --seqplatform       Sequencing platform type - "PacBio" (default) or "Illumina"
        --its_region        ITS part selector - "full" (defalut), "ITS1", "ITS2", "none" (trims primers only), or "ITS1_5.8S_ITS2"
        --primer_forward    Forward primer sequence (default, ITS9mun)
        --primer_reverse    Reverse primer sequence (default, ITS4ngsUni)
        --primer_mismatches
        --primer_foverlap   Min primer overlap (default, F primer length - 2)
        --primer_roverlap   Min primer overlap (default, R primer length - 2)
        --qc_maxn           Discard sequences with more than the specified number of N’s 
        --trim_minlen       Min sequence length after primer trimming (default, 10)
        --ITSx_tax          ITSx taxonomy profile (default, "all")
        --ITSx_evalue       ITSx E-value cutoff threshold (default, 1e-1)
        --ITSx_partial      Keep partial ITS sequences (defalt, off), otherwise specify min length cutoff
        --hp                Homopolymer compression (default, true)
        --hp_similarity     Allowed sequence similarity for homopolymer compression (default, 0.999)
        --hp_iddef          Sequence similarity definition for homopolymer compression (default, 2)
        
      # Chimera identification
        --chimera_db               Database for reference-based chimera removal
        --chimera_rescueoccurrence Min occurrence of chimeric sequences required to rescue them (default, 2)
        --chimeranov_abskew   De novo chimera identification `abskew` parameter (default, 2.0)
        --chimeranov_dn       De novo chimera identification `dn` parameter (default, 1.4)
        --chimeranov_mindiffs De novo chimera identification `mindiffs` parameter (default, 3)
        --chimeranov_mindiv   De novo chimera identification `mindiv` parameter (default, 0.8)
        --chimeranov_minh     De novo chimera identification `minh` parameter (default, 0.28)
        --chimeranov_xn       De novo chimera identification `xn` parameter (default, 8.0)

      # Tag-jump removal
        --tj_f      Tag-jump filtering, UNCROSS parameter `f` (default, 0.01)
        --tj_p      Tag-jump filtering parameter `p` (default, 1)
        --otu_id    Sequence similarity for OTU clustering (default, 0.98)
        --otu_iddef Sequence similarity definition for tag-jump removal step (default, 2)

      # PacBio-specific parameters
        --lima_barcodetype Tag type ("single", "dual", "dual_symmetric", "dual_asymmetric")
        --lima_minscore    Minimum barcode score for demultiplexing (default, 93)
        --lima_minendscore Minimum second barcode score (only for asymmetric and dual barcoding scheme; default, 50)
        --lima_minrefspan  Minimum read span relative to the barcode length (0-1; default, 0.75)
        --lima_minscoringregions Number of barcodes scored required for demultiplexing using dual barcodes (default, 2 = requires both barcodes)
        --lima_windowsize  Window size for barcode lookup (default, 70 bp)
        --lima_minlen      Minimum sequence length after clipping barcodes (default, 40)
        --qc_maxee         Maximum number of expected errors (default, false)
        --qc_maxeerate     Maximum number of expected errors per base (default, 0.01)
        --qc_maxhomopolymerlen  Threshold for a homopolymer region length in a sequence (default, 25)

      # Illumina-specific parameters
        --qc_avgphred      Average Phred score for QC (default, false)
        --qc_twocolor      Enable two-color chemistry mode, e.g. for Illumina NovaSeq (default, false)
        --qc_phredmin      Two-color mode: min Phred score of qualified bases (default, 24)
        --qc_phredperc     Two-color mode: Percentage of bases allowed to be unqualified (default, 30)
        --qc_polyglen      Two-color mode: minimum length of polyG tail (default, 8)
        --barcode_window   Window size for barcode lookup (default, 30 bp)
        --barcode_errors   Maximum allowed number of errors in barcodes (default, 1)
        --barcode_overlap  Min overlap between read and barcode (default, 11)
        --pe_minoverlap    Min length to detect overlapped region of PE reads (default, 20)
        --pe_difflimit     Max number of mismatched bases in PE overlap (default, 5)
        --pe_diffperclimit Max percentage of mismatched bases in PE overlap (default, 20)
        --pe_minlen        Min length of merged sequences (default, 30)
        --illumina_keep_notmerged  Keep not merged Illumina reads (default, true)
        --illumina_joinpadgap      Join not merged reads into one sequence using padding sequence string (default, NNNNNNNNNN)
        --illumina_joinpadqual     Join not merged reads into one sequence using padding quality string (default, IIIIIIIIII)

      # Miscellaneous parameters
        --gzip_compression Compression level for GZIP (default, 7; 1 = fastest (worst compression), 9 = slowest (best))

    NEXTFLOW-SPECIFIC:
        -profile       Configuration profile
        -resume        Execute the pipeline using the cached results (e.g., in case of )
        -work-dir      Path to the directory where intermediate result files are stored
        -qs            Queue size (max number of processes that can be executed in parallel); e.g., 8
        -r             Pipeline version to run (GitHub branch, tag, or SHA number)
    """.stripIndent()
}

