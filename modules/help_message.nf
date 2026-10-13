

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
        --qc_maxn           Discard sequences with more than the specified number of Ns 
        --trim_minlen       Min sequence length after primer trimming (default, 10)
        --itsx_tool         ITS extractor: "ITSx" (default, v1.x) or "ITSx2" (faster)
        --ITSx_chunk_size   Dereplicated sequences per chunk for parallel extraction (default, 9999)
        --ITSx_tax          ITSx taxonomy profile (default, "all"; ITSx 1.x only)
        --ITSx_evalue       ITSx E-value cutoff threshold (default, 1e-1; ITSx 1.x only)
        --ITSx_partial      Keep partial ITS sequences (default, off), otherwise specify min length cutoff (ITSx 1.x only)
        --hp                Homopolymer compression (default, true)
        --hp_similarity     Allowed sequence similarity for homopolymer compression (default, 0.999)
        --hp_iddef          Sequence similarity definition for homopolymer compression (default, 2)

      # Long-chimera removal (after primer trimming, before ITS extraction)
        --long_chimera         Remove sequences with rRNA in both orientations or with repeated rRNA regions (default, false)
        --long_chimera_minbits Min bit score of an rRNA hit (default, 50)
        --long_chimera_mincov  Min fraction of the rRNA model covered by a hit (default, 0.4)

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

      # Illumina-specific parameters (paired-end reads; also suitable for AVITI and DNBSEQ)
        Presets: `-profile miseq` (continuous Phred scores) or `-profile novaseq` (binned Phred scores, poly-G tails)
        --input_R1         File with multiplexed R1 reads (FASTQ.gz)
        --input_R2         File with multiplexed R2 reads (FASTQ.gz)
        --illumina_pe_pattern      File pattern of per-sample read pairs (with `--demultiplexed true`; default, "*_R{1,2}*.{fastq,fq}.gz")
        --illumina_barcodetype     Tag layout ("dual_symmetric" (default), "dual_asymmetric", "single")
        --illumina_demux_rescue    Rescue read pairs with a single readable tag (default, false)
        --illumina_revtag_orient   Orientation of reverse tags in `fwd...rev` format ("auto" (default), "forward", "revcomp")
        --barcode_window   Tag search window at the 5' end of reads (default, 30 bp)
        --barcode_errors   Maximum allowed number of mismatches in tags (default, 2)
        --barcode_overlap  Min overlap between read and tag (default, 11)
        --barcode_minlen   Min read length after tag removal (default, 50)
        --pe_merger        Read merging tool ("usearch" (default) or "vsearch")
        --pe_maxdiffs      Max number of mismatches in the overlap (default, 10)
        --pe_pctid         Min percent identity in the overlap (default, 80)
        --pe_minoverlap    Min overlap length (default, 16)
        --pe_minlen        Min length of merged reads (default, 50)
        --pe_maxlen        Max length of merged reads (default, none)
        --qc_polyglen      Trim 3' poly-G tails of at least this length prior to merging (default, 10)
        --illumina_keep_notmerged  Join not merged reads with an N-pad (default, false)
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

