/*
============================================================================
  NextITS: Pipeline to process eukaryotic ITS amplicons
============================================================================
  License: Apache-2.0
  Github : https://github.com/vmikk/NextITS
  Website: https://Next-ITS.github.io/
----------------------------------------------------------------------------
*/

// Subworkflow for Illumina paired-end reads:
//   quality-score check -> demultiplexing (cutadapt) -> reorientation -> read merging -> (optional) joining
// The resulting merged reads enter the common Step-1 workflow (QC, primer check, ITS extraction, ...)

include { illumina_qcheck; demux_illumina; reorient_pe; merge_pe; join_pe } from '../modules/Illumina_pe.nf'


// Report the detected type of Phred scores (binned or continuous),
// and warn if it does not match the selected preset (`-profile miseq` / `-profile novaseq`)
def report_quality_type(rows) {

    def binned = rows.every { r -> r.QualityType == "binned" }
    def qtype  = binned ? "binned" : "continuous"
    def qvals  = rows.collect { r -> "${r.Mate}: ${r.QValues}" }.join("; ")
    def polyg  = rows.collect { r -> r.PolyG_Percent as Double }.max()

    log.info "Illumina quality check: Phred scores look ${qtype} (${qvals}); reads with 3' poly-G tails: ${polyg}%"

    def suggested = binned ? "novaseq" : "miseq"
    if (params.qc_binned == null) {
        log.warn "No sequencing-platform preset was selected. Based on the quality scores, consider using `-profile ${suggested}`"
    } else if ((params.qc_binned as Boolean) != binned) {
        log.warn "The selected preset expects ${params.qc_binned ? 'binned' : 'continuous'} Phred scores, but the data look ${qtype}. Consider using `-profile ${suggested}`"
    }
    if (!binned && polyg > 1 && !params.qc_polyglen) {
        log.warn "More than 1% of reads end with a poly-G tail (${polyg}%), consider enabling poly-G trimming (`--qc_polyglen 10`)"
    }
}


workflow ILLUMINA_PE {

  take:
    ch_multiplexed   // tuple(R1, R2) with multiplexed reads (empty if the data are demultiplexed)
    ch_persample     // tuple(sampleID, [R1, R2]) per sample (empty if the data are multiplexed)
    ch_tags          // validated tags (single or symmetric dual tags)
    ch_tags_dual     // per-sample forward and reverse tags (dual tags) or a dummy file

  main:

    is_demultiplexed = params.demultiplexed

    if ( !is_demultiplexed ) {

      // Quality-score check on the raw data
      illumina_qcheck(ch_multiplexed)

      // Demultiplexing
      demux_illumina(ch_multiplexed, ch_tags, ch_tags_dual)

      // Per-sample read pairs, tuple(sampleID, [R1, R2])
      ch_pairs = demux_illumina.out.samples_demux
        .flatten()
        .map { f -> tuple(f.name.replaceAll(/_R[12]\.fq\.gz$/, ''), f) }
        .groupTuple(size: 2, sort: true)

      ch_demux_totals  = demux_illumina.out.totals
      ch_demux_summary = demux_illumina.out.summary

    } else {

      // Quality-score check on the first sample
      illumina_qcheck(ch_persample.take(1).map { _id, reads -> tuple(reads[0], reads[1]) })

      ch_pairs         = ch_persample
      ch_demux_totals  = channel.empty()
      ch_demux_summary = channel.empty()
    }

    illumina_qcheck.out.tsv
      .splitCsv(header: true, sep: '\t')
      .toList()
      .subscribe { rows -> report_quality_type(rows) }

    // Reorient read pairs (R1 = forward-primer strand)
    reorient_pe(ch_pairs)

    // Merge read pairs
    merge_pe(reorient_pe.out.reads)

    // Join non-merged read pairs (optional)
    if ( params.illumina_keep_notmerged ) {
      join_pe(merge_pe.out.notmerged)
      ch_joined = join_pe.out.joined
    } else {
      ch_joined = channel.empty()
    }

    // Per-sample stats
    ch_reorient_stats = reorient_pe.out.stats
      .collectFile(name: "Reorient_stats.tsv", keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/01_Demux")

    ch_merge_stats = merge_pe.out.stats
      .collectFile(name: "Merge_stats.tsv", keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/01_Demux")

  emit:
    merged         = merge_pe.out.merged      // per-sample merged reads, `{sampleID}.fq.gz`
    joined         = ch_joined                // per-sample joined reads, `{sampleID}_JoinedPE.fq.gz`
    reorient_stats = ch_reorient_stats        // per-sample input and reoriented pairs
    merge_stats    = ch_merge_stats           // per-sample merged and non-merged pairs
    demux_totals   = ch_demux_totals          // run totals of demultiplexing
    demux_summary  = ch_demux_summary         // per-sample demultiplexing stats
    qcheck         = illumina_qcheck.out.tsv  // quality-score check
}
