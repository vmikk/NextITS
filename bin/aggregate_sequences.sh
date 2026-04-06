#!/bin/bash

## Aggregate sequence parquet files, optionally filter de novo chimeras,
## recover likely false-positive chimeras, and export FASTA + Parquet outputs

# NB. currently, script creates a temporary DB in /tmp/aggregate_sequences.qWBByQ/aggregate.duckdb
## to do - move it from temp to current dir?

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: aggregate_sequences.sh --seqtabs DIR [options]

Required:
  --seqtabs DIR            Directory containing sequence tables in Parquet format

Options:
  --maxchim NUM|NA         Maximum de novo chimera score to keep (default: 0.6)
  --recoverdenovo BOOL     Recover de novo chimeras seen as non-chimeras elsewhere (default: TRUE)
  --output PREFIX          Output prefix (default: Seqs)
  --threads INT            Number of CPU threads to use (default: 4)
  --memory SIZE            DuckDB memory limit, e.g. 50GB (optional)
  --compression INT        Parquet ZSTD compression level 0-22 (default: 10)
  -h, --help               Show this help message
EOF
    exit 1
}

to_lower() {
    printf "%s" "$1" | tr '[:upper:]' '[:lower:]'
}

normalize_na() {
    local value
    value="$(to_lower "${1:-}")"
    case "$value" in
        ""|"na"|"null")
            printf ""
            ;;
        *)
            printf "%s" "$1"
            ;;
    esac
}

normalize_bool() {
    local value
    value="$(to_lower "${1:-}")"
    case "$value" in
        "true"|"t"|"1"|"yes")
            printf "true"
            ;;
        "false"|"f"|"0"|"no")
            printf "false"
            ;;
        *)
            echo "Error: Invalid logical value: ${1:-}" >&2
            usage
            ;;
    esac
}

escape_sql() {
    printf "%s" "$1" | sed "s/'/''/g"
}

calc_percent() {
    awk -v numerator="$1" -v denominator="$2" 'BEGIN {
        if (denominator == 0) {
            printf "0.0"
        } else {
            printf "%.1f", (numerator / denominator) * 100
        }
    }'
}

require_command() {
    if ! command -v "$1" >/dev/null 2>&1; then
        echo "Error: Required command not found: $1" >&2
        exit 1
    fi
}

## Defaults
SEQTABS=""
MAXCHIM="0.6"
RECOVERDENOVO="TRUE"
OUTPUT="Seqs"
THREADS="4"
MEMORY=""
COMPRESSION="10"

## Parse command-line options
while [[ $# -gt 0 ]]; do
    case "$1" in
        --seqtabs)
            [[ $# -ge 2 ]] || usage
            SEQTABS="$2"
            shift 2
            ;;
        --maxchim)
            [[ $# -ge 2 ]] || usage
            MAXCHIM="$2"
            shift 2
            ;;
        --recoverdenovo)
            [[ $# -ge 2 ]] || usage
            RECOVERDENOVO="$2"
            shift 2
            ;;
        --output)
            [[ $# -ge 2 ]] || usage
            OUTPUT="$2"
            shift 2
            ;;
        --threads)
            [[ $# -ge 2 ]] || usage
            THREADS="$2"
            shift 2
            ;;
        --memory)
            [[ $# -ge 2 ]] || usage
            MEMORY="$2"
            shift 2
            ;;
        --compression)
            [[ $# -ge 2 ]] || usage
            COMPRESSION="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Error: Unknown argument: $1" >&2
            usage
            ;;
    esac
done

## Normalize values
SEQTABS="$(normalize_na "$SEQTABS")"
MAXCHIM="$(normalize_na "$MAXCHIM")"
MEMORY="$(normalize_na "$MEMORY")"
RECOVERDENOVO="$(normalize_bool "$RECOVERDENOVO")"

## Validation
if [[ -z "$SEQTABS" ]]; then
    echo -e "Error: Input directory with quality-filtered sequences is not specified!\n" >&2
    usage
fi

if [[ ! -d "$SEQTABS" ]]; then
    echo "Error: Input directory not found: $SEQTABS" >&2
    exit 1
fi

if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -le 0 ]]; then
    echo "Error: Threads must be a positive integer!" >&2
    usage
fi

if ! [[ "$COMPRESSION" =~ ^[0-9]+$ ]] || [[ "$COMPRESSION" -lt 0 ]] || [[ "$COMPRESSION" -gt 22 ]]; then
    echo "Error: Compression level must be an integer between 0 and 22!" >&2
    usage
fi

if [[ -n "$MAXCHIM" ]] && ! [[ "$MAXCHIM" =~ ^-?[0-9]+([.][0-9]+)?$ ]]; then
    echo "Error: --maxchim must be numeric or NA/null!" >&2
    usage
fi

require_command duckdb
require_command pigz

## View user-supplied parameters
echo -e "\nInput parameters:"
echo "Path to sequence tables: $SEQTABS"
if [[ -n "$MAXCHIM" ]]; then
    echo "Max de novo chimera score: $MAXCHIM"
else
    echo "Max de novo chimera score: NA"
fi
echo "De novo chimera recovery: $RECOVERDENOVO"
echo "Output prefix: $OUTPUT"
echo "CPU threads: $THREADS"
if [[ -n "$MEMORY" ]]; then
    echo "DuckDB memory limit: $MEMORY"
fi
echo "Parquet compression level (ZSTD): $COMPRESSION"

echo
echo "Setting number of CPU threads to: $THREADS"

## Collect parquet files recursively
echo -e "\n..Looking for sequence tables"
shopt -s nullglob globstar
PARQUET_FILES=( "$SEQTABS"/**/*.parquet )
shopt -u nullglob globstar

if [[ "${#PARQUET_FILES[@]}" -eq 0 ]]; then
    echo "Error: No parquet files found in: $SEQTABS" >&2
    exit 1
fi

echo "... Tables found: ${#PARQUET_FILES[@]}"

PARQUET_SQL_LIST=""
for file_path in "${PARQUET_FILES[@]}"; do
    escaped_path="$(escape_sql "$file_path")"
    if [[ -n "$PARQUET_SQL_LIST" ]]; then
        PARQUET_SQL_LIST+=", "
    fi
    PARQUET_SQL_LIST+="'${escaped_path}'"
done
PARQUET_SQL_LIST="[${PARQUET_SQL_LIST}]"

TMPDIR_ROOT="$(mktemp -d -t aggregate_sequences.XXXXXX)"
DB_PATH="${TMPDIR_ROOT}/aggregate.duckdb"
trap 'rm -rf "$TMPDIR_ROOT"' EXIT

run_sql() {
    duckdb "$DB_PATH" -init /dev/null -batch -c "$1"
}

query_scalar() {
    duckdb "$DB_PATH" -init /dev/null -batch -noheader -list -c "$1"
}

OUTPUT_PARQUET="${OUTPUT}.parquet"
OUTPUT_PARQUET_SQL="$(escape_sql "$OUTPUT_PARQUET")"

SQL_SETUP=""
SQL_SETUP+="SET threads TO ${THREADS};"
if [[ -n "$MEMORY" ]]; then
    SQL_SETUP+="SET memory_limit = '$(escape_sql "$MEMORY")';"
fi

echo -e "\n..Loading sequence tables"
run_sql "
${SQL_SETUP}
CREATE OR REPLACE TABLE seq_table AS
SELECT
    SeqID___SampleID,
    SampleID,
    SeqID,
    Abundance,
    SeqLen,
    PhredScore,
    MaxEE,
    MEEP,
    DeNovo_Chimera,
    DeNovo_Chimera_Score,
    Sequence
FROM read_parquet(${PARQUET_SQL_LIST});
"

TOTAL_RECORDS="$(query_scalar "SELECT COUNT(*) FROM seq_table;")"
TOTAL_SEQUENCES="$(query_scalar "SELECT COUNT(DISTINCT Sequence) FROM seq_table;")"
TOTAL_SAMPLES="$(query_scalar "SELECT COUNT(DISTINCT SampleID) FROM seq_table;")"
TOTAL_READS="$(query_scalar "SELECT COALESCE(SUM(Abundance), 0) FROM seq_table;")"
MAX_SCORE_OBSERVED="$(query_scalar "SELECT COALESCE(CAST(MAX(DeNovo_Chimera_Score) AS VARCHAR), 'NA') FROM seq_table;")"

echo "... Total number of records: $TOTAL_RECORDS"
echo "... Total number unique sequences: $TOTAL_SEQUENCES"
echo "... Total number unique samples (fastq files): $TOTAL_SAMPLES"

if [[ -n "$MAXCHIM" ]]; then
    echo -e "\n..Filtering data by max de novo chimera score"
    echo "... Max de novo chimera score observed: $MAX_SCORE_OBSERVED"

    PUTATIVE_CHIMERA_RECORDS="$(query_scalar "SELECT COUNT(*) FROM seq_table WHERE DeNovo_Chimera_Score >= ${MAXCHIM};")"
    PUTATIVE_CHIMERA_SEQUENCES="$(query_scalar "SELECT COUNT(DISTINCT Sequence) FROM seq_table WHERE DeNovo_Chimera_Score >= ${MAXCHIM};")"
    PUTATIVE_CHIMERA_READS="$(query_scalar "SELECT COALESCE(SUM(Abundance), 0) FROM seq_table WHERE DeNovo_Chimera_Score >= ${MAXCHIM};")"

    echo "... Putative chimera records: $PUTATIVE_CHIMERA_RECORDS"
    echo "... Putative chimera sequences: $PUTATIVE_CHIMERA_SEQUENCES"
    echo "... Putative chimera reads: $PUTATIVE_CHIMERA_READS"

    if [[ "$RECOVERDENOVO" == "false" ]]; then
        run_sql "
        ${SQL_SETUP}
        CREATE OR REPLACE TABLE filtered_table AS
        SELECT
            SeqID___SampleID,
            SampleID,
            SeqID,
            Abundance,
            SeqLen,
            PhredScore,
            MaxEE,
            MEEP,
            DeNovo_Chimera,
            DeNovo_Chimera_Score,
            Sequence
        FROM seq_table
        WHERE DeNovo_Chimera_Score < ${MAXCHIM}
           OR DeNovo_Chimera_Score IS NULL;
        "
    else
        run_sql "
        ${SQL_SETUP}
        CREATE OR REPLACE TABLE chimera_rows AS
        SELECT
            SeqID___SampleID,
            DeNovo_Chimera_Score,
            Sequence,
            Abundance
        FROM seq_table
        WHERE DeNovo_Chimera_Score >= ${MAXCHIM};

        CREATE OR REPLACE TABLE nonchimera_rows AS
        SELECT s.*
        FROM seq_table AS s
        LEFT JOIN chimera_rows AS c
          ON s.SeqID___SampleID = c.SeqID___SampleID
        WHERE c.SeqID___SampleID IS NULL;

        CREATE OR REPLACE TABLE recovered_sequences AS
        SELECT DISTINCT c.Sequence
        FROM (SELECT DISTINCT Sequence FROM chimera_rows) AS c
        INNER JOIN (SELECT DISTINCT Sequence FROM nonchimera_rows) AS n
          ON c.Sequence = n.Sequence;

        CREATE OR REPLACE TABLE true_chimera_sequences AS
        SELECT DISTINCT c.Sequence
        FROM chimera_rows AS c
        LEFT JOIN recovered_sequences AS r
          ON c.Sequence = r.Sequence
        WHERE r.Sequence IS NULL;

        CREATE OR REPLACE TABLE filtered_table AS
        SELECT s.*
        FROM seq_table AS s
        LEFT JOIN true_chimera_sequences AS t
          ON s.Sequence = t.Sequence
        WHERE t.Sequence IS NULL;
        "

        RECOVERED_SEQUENCES="$(query_scalar "SELECT COUNT(*) FROM recovered_sequences;")"
        RECOVERED_RECORDS="$(query_scalar "SELECT COUNT(*) FROM chimera_rows WHERE Sequence IN (SELECT Sequence FROM recovered_sequences);")"
        RECOVERED_READS="$(query_scalar "SELECT COALESCE(SUM(Abundance), 0) FROM chimera_rows WHERE Sequence IN (SELECT Sequence FROM recovered_sequences);")"
        TRUE_CHIMERA_SEQUENCES="$(query_scalar "SELECT COUNT(*) FROM true_chimera_sequences;")"

        if [[ "$RECOVERED_SEQUENCES" -gt 0 ]]; then
            echo ".... Probably there are a few false-positive chimeras"
            echo ".... Recovering $RECOVERED_SEQUENCES sequences"
            echo ".... Recovered chimera records: $RECOVERED_RECORDS"
            echo ".... Recovered chimera reads: $RECOVERED_READS"
        fi
        echo "... True chimera sequences removed: $TRUE_CHIMERA_SEQUENCES"
    fi
else
    echo -e "\n..Skipping de novo chimera filtering (max score is NA)"
    run_sql "
    ${SQL_SETUP}
    CREATE OR REPLACE TABLE filtered_table AS
    SELECT
        SeqID___SampleID,
        SampleID,
        SeqID,
        Abundance,
        SeqLen,
        PhredScore,
        MaxEE,
        MEEP,
        DeNovo_Chimera,
        DeNovo_Chimera_Score,
        Sequence
    FROM seq_table;
    "
fi

FILTERED_RECORDS="$(query_scalar "SELECT COUNT(*) FROM filtered_table;")"
FILTERED_SEQUENCES="$(query_scalar "SELECT COUNT(DISTINCT Sequence) FROM filtered_table;")"
FILTERED_SAMPLES="$(query_scalar "SELECT COUNT(DISTINCT SampleID) FROM filtered_table;")"
FILTERED_READS="$(query_scalar "SELECT COALESCE(SUM(Abundance), 0) FROM filtered_table;")"

RECORDS_REMOVED=$(( TOTAL_RECORDS - FILTERED_RECORDS ))
READS_REMOVED=$(( TOTAL_READS - FILTERED_READS ))

echo "... Records removed: ${RECORDS_REMOVED} ($(calc_percent "$RECORDS_REMOVED" "$TOTAL_RECORDS")%)"
echo "... Reads removed: ${READS_REMOVED} ($(calc_percent "$READS_REMOVED" "$TOTAL_READS")%)"
echo "... Filtered records retained: $FILTERED_RECORDS"
echo "... Filtered unique sequences: $FILTERED_SEQUENCES"
echo "... Filtered unique samples: $FILTERED_SAMPLES"

echo -e "\n..Sorting table by abundance, quality score"
run_sql "
${SQL_SETUP}
CREATE OR REPLACE TABLE filtered_sorted AS
SELECT
    SeqID___SampleID,
    SampleID,
    SeqID,
    Abundance,
    SeqLen,
    PhredScore,
    MaxEE,
    MEEP,
    DeNovo_Chimera,
    DeNovo_Chimera_Score,
    Sequence
FROM filtered_table
ORDER BY Abundance DESC, PhredScore DESC, SeqID;
"

echo "..Preparing FASTA file"
echo "..Exporting FASTA file with filtered sequences"
duckdb "$DB_PATH" -init /dev/null -batch -noheader -list -c "
${SQL_SETUP}
SELECT
    '>' || SeqID || ';size=' || CAST(Abundance AS VARCHAR) || chr(10) || Sequence
FROM filtered_sorted
ORDER BY Abundance DESC, PhredScore DESC, SeqID;
" | pigz -p"$THREADS" > "${OUTPUT}.fa.gz"

echo "..Exporting filtered table"
run_sql "
${SQL_SETUP}
COPY (
    SELECT
        SeqID___SampleID,
        SampleID,
        SeqID,
        Abundance,
        SeqLen,
        PhredScore,
        MaxEE,
        MEEP,
        DeNovo_Chimera,
        DeNovo_Chimera_Score,
        Sequence
    FROM filtered_sorted
) TO '${OUTPUT_PARQUET_SQL}'
(FORMAT PARQUET, COMPRESSION zstd, COMPRESSION_LEVEL ${COMPRESSION});
"

echo "All done."
