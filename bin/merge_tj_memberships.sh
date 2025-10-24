#!/bin/bash

## Usage:
# merge_tj_memberships.sh \
#   -d 'Dereplicated.parquet' \
#   -c 'Clustered.parquet' \
#   -o 'TJPreclust.uc.parquet' \
#   -t 4

## Input data:
# - Parsed UC file from dereplication (`Dereplicated.parquet`)
# - Parsed UC file from clustering (`Clustered.parquet`)

## Function to display usage information
usage() {
    echo "Usage: $0 -d DEREP -c CLUST -o OUTPUT [-t THREADS]"
    echo "  -d DEREP    : Parquet file from dereplication"
    echo "  -c CLUST    : Parquet file from clustering"
    echo "  -o OUTPUT   : Output Parquet file path"
    echo "  -t THREADS  : Number of CPU threads to use (optional)"
    exit 1
}

## Initialize variables
DEREP=""
CLUST=""
OUTPUT="TJPreclust.uc.parquet"   # default output file name
THREADS=""

## Parse command-line options
while getopts "d:c:o:t:" opt; do
    case $opt in
        d) DEREP="$OPTARG" ;;
        c) CLUST="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        *) usage ;;
    esac
done


## Validate input parameters
if [[ -z "$DEREP" || -z "$CLUST" ]]; then
    echo -e "Error: Missing required parameters!\n"
    usage
fi

## Threads should be a positive integer
if [[ -n "$THREADS" && "$THREADS" -le 0 ]]; then
    echo -e "Error: Threads must be a positive integer!\n"
    usage
fi

## View user-supplied parameters
echo -e "\nInput parameters:"
echo "Parquet file from dereplication: $DEREP"
echo "Parquet file from clustering: $CLUST"
echo "Output file: $OUTPUT"
if [[ -n "$THREADS" ]]; then
    echo "Threads: $THREADS"
fi

SQL_COMMAND=""

## Add configuration settings (if provided)
if [[ -n "$THREADS" ]]; then
    SQL_COMMAND+="
SET threads TO ${THREADS};
"
fi


SQL_COMMAND+="
COPY (
  SELECT
    d.query   AS SeqID,
    c.target  AS OTU
  FROM read_parquet('${DEREP}') AS d
  LEFT JOIN read_parquet('${CLUST}') AS c
    ON d.target = c.query
) TO '${OUTPUT}' (FORMAT PARQUET, COMPRESSION 'ZSTD', COMPRESSION_LEVEL 8);
"

## Execute the SQL command
echo -e "\nExecuting DuckDB command"

duckdb -c "${SQL_COMMAND}"
