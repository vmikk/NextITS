/*
 * Step-1 results used in Step-2 (one or multiple sequencing runs)
 *
 * For each run, finds:
 *   - rRNA parts extracted by ITSx/ITSx2
 *       `03_ITSx_PooledParts/rRNA_parts.parquet`  (`parts_mode = parquet`)
 *       `03_ITSx_PooledParts/*.fasta.gz`          (`parts_mode = pooled`,    NextITS < 1.2.0)
 *       `03_ITSx/<sample>.<region>.fasta.gz`      (`parts_mode = persample`, older versions)
 *   - per-read sequence qualities (`09_DB/SeqQualities.parquet`)
 *   - Step-1 parameters (`pipeline_info/pipeline_params.tsv`) - which rRNA region was used for clustering
 */

// Column of the rRNA parts table with the region used in Step-2 (see `ITS_EXTRACTION` subworkflow)
def regionColumn(its_region, itsx_tool, itsx_partial) {
    def col = [
        "full": "ITS", "ITS1_5.8S_ITS2": "ITSnf",
        "ITS1": "ITS1", "ITS2": "ITS2", "SSU": "SSU", "LSU": "LSU"
        ][its_region]
    def use_partial = (itsx_tool == "ITSx") && itsx_partial.toString().toInteger() != 0
    if( col && use_partial && its_region in ["ITS1", "ITS2", "SSU", "LSU"] ){
        col = col + "_partial"
    }
    return col ?: "none"
}

// Find Step-1 runs, returns a list of maps (one per run)
def step1Runs(String base) {

    def single = file("${base}/07_SeqTable/Seqs.parquet")
    def dirs = single.exists()
        ? [ file(base) ]
        : files("${base}/*/07_SeqTable/Seqs.parquet").collect { it.parent.parent }.sort { it.name }

    def runs = []
    def ids  = [] as Set

    dirs.each { dir ->

        // Run ID (made unique, used for file names)
        def id = dir.name.replaceAll(/[^A-Za-z0-9_.-]/, '_')
        if( id in ids ){ id = "${id}_${runs.size() + 1}" }
        ids << id

        // Step-1 parameters (headerless `key<TAB>value`)
        def pp = dir.resolve("pipeline_info/pipeline_params.tsv")
        def s1 = pp.exists()
            ? pp.readLines().findAll { it.contains('\t') }.collectEntries { line -> def t = line.split('\t', 2); [ (t[0]): t[1] ] }
            : [:]
        def known = s1.containsKey("its_region")
        def its_region   = s1.its_region   ?: "full"
        def itsx_tool    = s1.itsx_tool    ?: "ITSx"
        def itsx_partial = s1.ITSx_partial?.isInteger() ? s1.ITSx_partial : "0"

        // rRNA parts
        def parts_parquet = dir.resolve("03_ITSx_PooledParts/rRNA_parts.parquet")
        def parts_pooled  = files("${dir}/03_ITSx_PooledParts/*.fasta.gz")
        def parts_sample  = files("${dir}/03_ITSx/*.{full,SSU,ITS1,5_8S,ITS2,LSU}.fasta.gz")
        def parts_mode = "none"
        def parts = null
        if( parts_parquet.exists() ){
            parts_mode = "parquet"
            parts = [ parts_parquet ]
        } else if( parts_pooled ){
            parts_mode = "pooled"
            parts = parts_pooled
        } else if( parts_sample ){
            parts_mode = "persample"
            parts = [ dir.resolve("03_ITSx") ]
        }

        // Sequence qualities
        def quals = dir.resolve("09_DB/SeqQualities.parquet")
        quals = quals.exists() ? quals : null

        runs << [
            id: id, dir: dir,
            its_region: its_region, itsx_tool: itsx_tool, ITSx_partial: itsx_partial, params_known: known,
            region_col: regionColumn(its_region, itsx_tool, itsx_partial),
            parts_mode: parts_mode, parts: parts, quals: quals
        ]
    }
    return runs
}

// Warnings about the Step-1 results (used for rRNA parts export)
// Only cases that may need user action: data missing for some runs (but present for others),
// unknown Step-1 parameters of runs with rRNA parts, different extraction settings
// (if data are missing for all runs, see `step1RunsInfo`)
def step1RunsWarnings(runs) {
    def msg = []
    def with_parts = runs.findAll { it.parts_mode != "none" }
    with_parts.findAll { !it.params_known }.each {
        msg << "Step-1 parameters not found for `${it.dir}` (pipeline_info/pipeline_params.tsv), assuming `--its_region full` for its rRNA parts"
    }
    if( with_parts.collect { [it.its_region, it.itsx_tool, it.ITSx_partial] }.unique().size() > 1 ){
        msg << "Step-1 runs differ in the rRNA region or extractor settings: " +
            with_parts.collect { "${it.id} (${it.its_region}, ${it.itsx_tool}, partial=${it.ITSx_partial})" }.join(", ")
    }
    def no_parts = runs.findAll { it.parts_mode == "none" }
    if( no_parts && with_parts ){
        msg << "No rRNA parts found for ${no_parts.size()} of ${runs.size()} Step-1 runs, their sequences will not contribute rRNA parts: " +
            no_parts.collect { it.its_region == "none" ? "${it.id} (`--its_region none`)" : it.id }.join(", ")
    }
    def no_quals = runs.findAll { !it.quals }
    if( no_quals && no_quals.size() < runs.size() ){
        msg << "No sequence qualities (09_DB/SeqQualities.parquet) found for ${no_quals.size()} of ${runs.size()} Step-1 runs, full-length sequences will not be available for them: " +
            no_quals.collect { it.id }.join(", ")
    }
    return msg
}

// Informational messages about the Step-1 results (data missing for all runs, e.g. older Step-1 versions)
def step1RunsInfo(runs) {
    def has_parts = runs.any { it.parts_mode != "none" }
    def has_quals = runs.any { it.quals }
    def all_none  = runs && runs.every { it.its_region == "none" }
    def why_parts = all_none ? " (Step-1 was run with `--its_region none`)" : ""
    if( !has_parts && !has_quals ){
        return [ "Export of rRNA parts skipped: the Step-1 results contain no rRNA parts or sequence qualities" + why_parts ]
    }
    if( !has_parts ){
        return [ "The Step-1 results contain no rRNA parts" + why_parts + ", only full-length sequences of OTUs will be exported" ]
    }
    if( !has_quals ){
        return [ "The Step-1 results contain no sequence qualities, full-length sequences of OTUs will not be exported" ]
    }
    return []
}

// Per-run summary of the Step-1 results (for the debug log)
def step1RunsDebug(runs) {
    return runs.collect { r ->
        "Step-1 run ${r.id}: dir=${r.dir}, its_region=${r.its_region}, itsx_tool=${r.itsx_tool}, ITSx_partial=${r.ITSx_partial}, " +
        "params_known=${r.params_known}, parts=${r.parts_mode}, quals=${r.quals ? 'yes' : 'no'}"
    }
}
