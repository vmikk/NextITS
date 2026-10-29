#!/bin/bash

## Prepare rRNA parts and sequence qualities of a single Step-1 run for the Step-2 rRNA parts export

## Usage:
# rrna_parts_prep.sh \
#   -r RunID \
#   -p parquet -i parts/ \
#   -c ITS \
#   -q quals/SeqQualities.parquet \
#   -t 4 -m 8GB

## Input data:
# - rRNA parts (`-i`), depending on the mode (`-p`):
#     parquet   = directory with `rRNA_parts.parquet` (NextITS >= 1.2.0)
#     pooled    = directory with pooled FASTA files from `03_ITSx_PooledParts` (e.g., `SSU.fasta.gz`, `ITS_Full.fasta.gz`)
#     persample = directory with `03_ITSx` subdirectory (per-sample FASTA files, e.g. `<sample>.SSU.fasta.gz`)
#     none      = no rRNA parts
# - Column with the rRNA region used for clustering (`-c`; e.g., ITS, ITSnf, ITS2, ITS2_partial).
#   If the parts table contains this information in its metadata, it takes precedence.
# - Sequence qualities from Step-1 (`-q`, `09_DB/SeqQualities.parquet`, optional)

## Outputs:
# - `<RunID>.parts.parquet` = rRNA parts table + `RunID` and `RegionKey` (sequence of the region used for clustering)
# - `<RunID>.quals.parquet` = one row per sequence hash:
#     SeqID, RunID, NReads (number of reads), and the read with the highest quality
#     (SampleID, ReadID, AvgPhredScore, MaxEE, MEEP, FullSeq = full-length primer-trimmed sequence)

usage() {
    echo "Usage: $0 -r RUN_ID -p MODE -i PARTS_DIR -c REGION_COLUMN [-q QUALITIES] [-t THREADS] [-m MEMORY]"
    exit 1
}

RUNID=""
MODE="none"
PARTS=""
REGCOL="none"
QUALS=""
THREADS=""
MEMORY=""

while getopts "r:p:i:c:q:t:m:h" opt; do
    case $opt in
        r) RUNID="$OPTARG" ;;
        p) MODE="$OPTARG" ;;
        i) PARTS="$OPTARG" ;;
        c) REGCOL="$OPTARG" ;;
        q) QUALS="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        m) MEMORY="$OPTARG" ;;
        *) usage ;;
    esac
done

if [[ -z "$RUNID" ]]; then
    echo -e "Error: Run ID is required!\n"
    usage
fi

SETTINGS="SET preserve_insertion_order = false; SET temp_directory = 'duckdb_tmp';"
if [[ -n "$THREADS" ]]; then SETTINGS+=" SET threads = ${THREADS};"; fi
if [[ -n "$MEMORY"  ]]; then SETTINGS+=" SET memory_limit = '${MEMORY}';"; fi

PARTSARGS=( -t "${THREADS:-1}" )
if [[ -n "$MEMORY" ]]; then PARTSARGS+=( -m "$MEMORY" ); fi

echo -e "Preparing rRNA parts of run ${RUNID} (mode: ${MODE})"

## Convert rRNA parts to Parquet (FASTA files from older NextITS versions)
SRC=""
case "$MODE" in
    parquet)
        SRC="${PARTS}/rRNA_parts.parquet"
        ;;
    pooled)
        SRC="tmp_parts.parquet"
        ITSx_parts_to_parquet.sh -o "$SRC" "${PARTSARGS[@]}" \
          -r ITS="${PARTS}/ITS_Full.fasta.gz" \
          -r SSU="${PARTS}/SSU.fasta.gz" \
          -r ITS1="${PARTS}/ITS1.fasta.gz" \
          -r S58="${PARTS}/5_8S.fasta.gz" \
          -r ITS2="${PARTS}/ITS2.fasta.gz" \
          -r LSU="${PARTS}/LSU.fasta.gz" \
          -r SSU_partial="${PARTS}/SSU_full_and_partial.fasta.gz" \
          -r ITS1_partial="${PARTS}/ITS1_full_and_partial.fasta.gz" \
          -r S58_partial="${PARTS}/5_8S_full_and_partial.fasta.gz" \
          -r ITS2_partial="${PARTS}/ITS2_full_and_partial.fasta.gz" \
          -r LSU_partial="${PARTS}/LSU_full_and_partial.fasta.gz"
        ;;
    persample)
        SRC="tmp_parts.parquet"
        D="${PARTS}/03_ITSx"
        ITSx_parts_to_parquet.sh -o "$SRC" "${PARTSARGS[@]}" \
          -r ITS="${D}/*.full.fasta.gz" \
          -r SSU="${D}/*.SSU.fasta.gz" \
          -r ITS1="${D}/*.ITS1.fasta.gz" \
          -r S58="${D}/*.5_8S.fasta.gz" \
          -r ITS2="${D}/*.ITS2.fasta.gz" \
          -r LSU="${D}/*.LSU.fasta.gz" \
          -r SSU_partial="${D}/*.SSU.full_and_partial.fasta.gz" \
          -r ITS1_partial="${D}/*.ITS1.full_and_partial.fasta.gz" \
          -r S58_partial="${D}/*.5_8S.full_and_partial.fasta.gz" \
          -r ITS2_partial="${D}/*.ITS2.full_and_partial.fasta.gz" \
          -r LSU_partial="${D}/*.LSU.full_and_partial.fasta.gz" \
          -r ITSnf="${D}/*_ITS1_58S_ITS2.fasta.gz"
        ;;
    none)
        ;;
    *)
        echo "Error: unknown mode '${MODE}'"
        exit 1
        ;;
esac

if [[ -n "$SRC" ]]; then

    ## The region column stored in the table metadata takes precedence
    META=$(duckdb -noheader -list -c "SELECT value::VARCHAR FROM parquet_kv_metadata('${SRC}') WHERE key::VARCHAR = 'region_column'" 2>/dev/null | head -n 1)
    if [[ -n "$META" ]]; then
        REGCOL="$META"
    fi
    echo "..Region used for clustering: ${REGCOL}"

    ## Sequence of the region used for clustering
    ## (ITS1-5.8S-ITS2 is not stored in the older versions, use the concatenated regions instead)
    case "$REGCOL" in
        SSU|ITS1|S58|ITS2|LSU|ITS|SSU_partial|ITS1_partial|S58_partial|ITS2_partial|LSU_partial)
            REGKEY="${REGCOL}" ;;
        ITSnf)
            REGKEY="coalesce(ITSnf, ITS1 || S58 || ITS2)" ;;
        *)
            REGKEY="NULL::VARCHAR" ;;
    esac

    duckdb -c "
    ${SETTINGS}
    COPY (
      SELECT '${RUNID}' AS RunID, *, ${REGKEY} AS RegionKey
      FROM read_parquet('${SRC}')
    ) TO '${RUNID}.parts.parquet' (FORMAT PARQUET, COMPRESSION 'ZSTD', COMPRESSION_LEVEL 3);
    SELECT count(*) AS Sequences, count(RegionKey) AS With_region_sequence FROM '${RUNID}.parts.parquet';
    "

    if [[ -f tmp_parts.parquet ]]; then rm tmp_parts.parquet; fi
fi


## Summarize sequence qualities by hash (keep the read with the highest quality)
if [[ -n "$QUALS" && -s "$QUALS" ]]; then

    echo -e "\nSummarizing sequence qualities"

    ## Columns could be missing in the older versions
    COLS=$(duckdb -noheader -list -c "SELECT string_agg(column_name, ',') FROM (DESCRIBE SELECT * FROM read_parquet('${QUALS}'))")
    col_or_null() {
        if [[ ",${COLS}," == *",$1,"* ]]; then echo "$1"; else echo "NULL::$2"; fi
    }
    SAMPLE=$(col_or_null SampleID VARCHAR)
    READID=$(col_or_null PacBioID VARCHAR)
    PHRED=$(col_or_null AvgPhredScore DOUBLE)
    MAXEE=$(col_or_null MaxEE DOUBLE)
    MEEP=$(col_or_null MEEP DOUBLE)

    duckdb -c "
    ${SETTINGS}
    COPY (
      WITH b AS (
        SELECT
          Hash AS SeqID,
          count(*) AS NReads,
          arg_max(
            {'SampleID': ${SAMPLE}, 'ReadID': ${READID}, 'AvgPhredScore': ${PHRED}, 'MaxEE': ${MAXEE}, 'MEEP': ${MEEP}, 'FullSeq': upper(Sequence)},
            (coalesce(${PHRED}, 0), -coalesce(${MAXEE}, 0), ${READID})
          ) AS Best
        FROM read_parquet('${QUALS}')
        GROUP BY Hash
      )
      SELECT SeqID, '${RUNID}' AS RunID, NReads, unnest(Best)
      FROM b
    ) TO '${RUNID}.quals.parquet' (FORMAT PARQUET, COMPRESSION 'ZSTD', COMPRESSION_LEVEL 3);
    SELECT count(*) AS Hashes, sum(NReads) AS Reads FROM '${RUNID}.quals.parquet';
    "
fi

if [[ -d duckdb_tmp ]]; then rm -rf duckdb_tmp; fi

echo "..Done"
