#!/bin/bash

## Pool rRNA regions extracted by ITSx (or ITSx2) into a single wide Parquet table

## Usage:
# ITSx_parts_to_parquet.sh \
#   -o rRNA_parts.parquet \
#   -r ITS=full -r SSU=ssu -r ITS1=its1 -r S58=58s -r ITS2=its2 -r LSU=lsu \
#   -r SSU_partial=ssu_partial \
#   -k its_region=full -k itsx_tool=ITSx -k ITSx_partial=0 -k region_column=ITS \
#   -t 4 -m 8GB

## Input data:
# - FASTA files with extracted regions, headers = "hash;size=N" (hash = SHA1 of the full primer-trimmed sequence)
# - `-r COL=PATH`, where PATH is
#     a directory (all `*.fasta.gz` files in it are used; placeholder files are ignored),
#     a file name pattern (e.g., `03_ITSx/*.SSU.fasta.gz`; the wildcard must be in the file name), or
#     a single file

## Output:
# - Parquet table with one row per sequence hash:
#   SeqID, Abundance, SSU, ITS1, S58, ITS2, LSU, ITS, ITSnf,
#   SSU_partial, ITS1_partial, S58_partial, ITS2_partial, LSU_partial
#   (NULL if the region was not detected)
# - `ITS` = ITSx `full` region, `ITSnf` = ITS1-5.8S-ITS2 cut by `extract_itsx_regions.R` (`--its_region ITS1_5.8S_ITS2`)
# - `Abundance` = number of reads (sum of `size=` across samples)
# - Key-value metadata (`-k key=value`) is stored in the Parquet file footer
#   (could be read with `SELECT * FROM parquet_kv_metadata('rRNA_parts.parquet')`)

## Notes
# - The same hash may be present in multiple samples. Region sequences are usually identical,
#   but the ITSx detection may differ between chunks (E-values depend on the chunk content).
#   In this case, the most frequent sequence (then, the most abundant, the longest) is kept.

## Function to display usage information
usage() {
    echo "Usage: $0 -o OUTPUT -r COL=PATH [-r COL=PATH ...] [-k KEY=VALUE ...] [-t THREADS] [-m MEMORY] [-x TEMP_DIR] [-z COMPRESSION]"
    echo "  -o OUTPUT      : Output Parquet file path"
    echo "  -r COL=PATH    : rRNA region (column name) and directory, file name pattern, or file with sequences"
    echo "                   COL = SSU, ITS1, S58, ITS2, LSU, ITS, ITSnf, SSU_partial, ITS1_partial, S58_partial, ITS2_partial, LSU_partial"
    echo "  -k KEY=VALUE   : Metadata to store in the Parquet file (optional)"
    echo "  -t THREADS     : Number of CPU threads to use (optional)"
    echo "  -m MEMORY      : Memory limit (e.g., '100GB') (optional)"
    echo "  -x TEMP_DIR    : Temporary directory for DuckDB (optional, default: ./duckdb_tmp)"
    echo "  -z COMPRESSION : ZSTD compression level (0-22) (optional, default: 12)"
    exit 1
}

VALID_COLS=("SSU" "ITS1" "S58" "ITS2" "LSU" "ITS" "ITSnf" "SSU_partial" "ITS1_partial" "S58_partial" "ITS2_partial" "LSU_partial")

## Initialize variables
OUTPUT=""
REGIONS=()
KVMETA=()
THREADS=""
MEMORY=""
TEMPDIR="duckdb_tmp"
COMPRESSION="12"

## Parse command-line options
while getopts "o:r:k:t:m:x:z:h" opt; do
    case $opt in
        o) OUTPUT="$OPTARG" ;;
        r) REGIONS+=("$OPTARG") ;;
        k) KVMETA+=("$OPTARG") ;;
        t) THREADS="$OPTARG" ;;
        m) MEMORY="$OPTARG" ;;
        x) TEMPDIR="$OPTARG" ;;
        z) COMPRESSION="$OPTARG" ;;
        *) usage ;;
    esac
done

## Validate input parameters
if [[ -z "$OUTPUT" || ${#REGIONS[@]} -eq 0 ]]; then
    echo -e "Error: Missing required parameters!\n"
    usage
fi

for reg in "${REGIONS[@]}"; do
    col="${reg%%=*}"
    if [[ ! " ${VALID_COLS[*]} " =~ " ${col} " ]] || [[ "$reg" != *=* ]]; then
        echo -e "Error: Invalid region specification '${reg}'\n"
        usage
    fi
done

if ! [[ "$COMPRESSION" =~ ^[0-9]+$ ]] || [ "$COMPRESSION" -lt 0 ] || [ "$COMPRESSION" -gt 22 ]; then
    echo -e "Error: Compression level must be an integer between 0 and 22!\n"
    usage
fi


## List FASTA files for a region (directory, file name pattern, or file)
list_files() {
    local path="$1"
    if [[ -d "$path" ]]; then
        find -L "$path" -maxdepth 1 -name "*.fasta.gz" -size +0 -print0
    elif [[ "$(basename "$path")" == *"*"* ]]; then
        find -L "$(dirname "$path")" -maxdepth 1 -name "$(basename "$path")" -size +0 -print0
    elif [[ -s "$path" ]]; then
        printf '%s\0' "$path"
    fi
}

## Convert sequences of all regions to a single tab-delimited stream (Col - ID - Sequence)
stream_regions() {
    for reg in "${REGIONS[@]}"; do
        col="${reg%%=*}"
        path="${reg#*=}"
        nf=$(list_files "$path" | tr -cd '\0' | wc -c)
        echo "..${col}: ${nf} file(s) [${path}]" >&2
        list_files "$path" \
          | xargs -0 -r seqkit fx2tab -i \
          | cut -f1,2 \
          | awk -v c="$col" 'BEGIN{OFS="\t"} {print c, $0}'
    done
}

## DuckDB settings
TEMPDIR_EXISTED=false
if [[ -e "$TEMPDIR" ]]; then
    TEMPDIR_EXISTED=true
fi
SETTINGS="SET preserve_insertion_order = false; SET temp_directory = '${TEMPDIR}';"
if [[ -n "$THREADS" ]]; then
    SETTINGS+=" SET threads = ${THREADS};"
fi
if [[ -n "$MEMORY" ]]; then
    SETTINGS+=" SET memory_limit = '${MEMORY}';"
fi

## Parquet metadata
KV=""
for kv in "${KVMETA[@]}"; do
    KV+="${KV:+, }${kv%%=*}: '${kv#*=}'"
done
if [[ -n "$KV" ]]; then
    KV=", KV_METADATA {${KV}}"
fi

## Wide table columns
WIDE=""
for col in "${VALID_COLS[@]}"; do
    WIDE+=", max(Seq) FILTER (WHERE Col = '${col}') AS ${col}"
done

echo -e "\nPooling rRNA regions into ${OUTPUT}"

stream_regions \
  | duckdb -c "
${SETTINGS}

CREATE TABLE lng AS
SELECT
  Col,
  split_part(Id, ';', 1) AS SeqID,
  coalesce(TRY_CAST(regexp_extract(Id, ';size=([0-9]+)', 1) AS BIGINT), 1) AS Size,
  upper(Seq) AS Seq
FROM read_csv('/dev/stdin',
  header = false, delim = '\t', quote = '', escape = '', auto_detect = false,
  columns = {'Col': 'VARCHAR', 'Id': 'VARCHAR', 'Seq': 'VARCHAR'});

CREATE TABLE cnt AS
SELECT Col, SeqID, Seq, count(*) AS N, sum(Size) AS S
FROM lng
GROUP BY ALL;

SELECT Col AS Region, count(DISTINCT SeqID) AS Sequences, sum(S) AS Reads
FROM cnt GROUP BY Col ORDER BY Col;

SELECT Col AS Region, count(*) AS Hashes_with_conflicting_sequences
FROM (SELECT Col, SeqID FROM cnt GROUP BY ALL HAVING count(*) > 1)
GROUP BY Col ORDER BY Col;

COPY (
  WITH pick AS (
    SELECT Col, SeqID, Seq
    FROM cnt
    QUALIFY row_number() OVER (PARTITION BY Col, SeqID ORDER BY N DESC, S DESC, length(Seq) DESC, Seq) = 1
  ),
  ab AS (
    SELECT SeqID, CAST(max(S) AS BIGINT) AS Abundance
    FROM (SELECT Col, SeqID, sum(S) AS S FROM cnt GROUP BY ALL)
    GROUP BY SeqID
  )
  SELECT SeqID, any_value(Abundance) AS Abundance ${WIDE}
  FROM pick JOIN ab USING (SeqID)
  GROUP BY SeqID
  ORDER BY SeqID
) TO '${OUTPUT}' (FORMAT PARQUET, COMPRESSION 'ZSTD', COMPRESSION_LEVEL ${COMPRESSION}, ROW_GROUP_SIZE 100000 ${KV});
"

## Remove the temporary directory, if it was created by DuckDB
if [[ "$TEMPDIR_EXISTED" == false && -d "$TEMPDIR" ]]; then
    rm -rf "$TEMPDIR"
fi

if [[ ! -f "$OUTPUT" ]]; then
    echo "Error: Output file was not created"
    exit 1
fi

echo "..Done"
