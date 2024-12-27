#!/bin/bash

## Usage:
# merge_hash_tables.sh \
#   -i '/path/to/input/directory' \
#   -o '/path/to/output.parquet' \
#   -t 4

## Input data:
# - Tab-delimited tables with columns:
#   SampleID - Hash - PacBioID - AvgPhredScore - MaxEE - MEEP - Sequence - Quality - Length

## Function to display usage information
usage() {
    echo "Usage: $0 -i INPUTDIR -o OUTPUT [-t THREADS] [-m MEMORY] [-x TEMP_DIR] [-z COMPRESSION]"
    echo "  -i INPUTDIR    : Input directory with text files"
    echo "  -o OUTPUT      : Output Parquet file path"
    echo "  -t THREADS     : Number of CPU threads to use (optional)"
    echo "  -z COMPRESSION: ZSTD compression level (0-22) (optional, default: 14)"
    exit 1
}

## Initialize variables
INPUT=""
OUTPUT=""
THREADS=""
COMPRESSION="14"

## Parse command-line options
while getopts "i:o:t:z:" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        z) COMPRESSION="$OPTARG" ;;
        *) usage ;;
    esac
done


## Validate input parameters
if [[ -z "$INPUT" || -z "$OUTPUT" ]]; then
    echo -e "Error: Missing required parameters!\n"
    usage
fi

## Threads should be a positive integer
if [[ -n "$THREADS" && "$THREADS" -le 0 ]]; then
    echo -e "Error: Threads must be a positive integer!\n"
    usage
fi

## Validate compression level
if ! [[ "$COMPRESSION" =~ ^[0-9]+$ ]] || [ "$COMPRESSION" -lt 0 ] || [ "$COMPRESSION" -gt 22 ]; then
    echo -e "Error: Compression level must be an integer between 0 and 22!\n"
    usage
fi

## View user-supplied parameters
echo -e "\nInput parameters:"
echo "Input directory: $INPUT"
echo "Output file: $OUTPUT"
if [[ -n "$THREADS" ]]; then
    echo "Threads: $THREADS"
fi
echo "Parquet compression level (ZSTD): $COMPRESSION"


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
        column0 as SampleID,
        column1 as Hash,
        column2 as PacBioID,
        column3 as AvgPhredScore,
        column4 as MaxEE,
        column5 as MEEP,
        column6 as Sequence,
        column7 as Quality,
        column8 as Length
    FROM read_csv('${INPUT}/*.txt.gz', 
        header  = false, 
        delim   = '\t',
        columns = {
            'column0': 'VARCHAR',
            'column1': 'VARCHAR',
            'column2': 'VARCHAR',
            'column3': 'DOUBLE',
            'column4': 'DOUBLE',
            'column5': 'DOUBLE',
            'column6': 'VARCHAR',
            'column7': 'VARCHAR',
            'column8': 'INTEGER'
        }
    )
) TO '${OUTPUT}' (FORMAT PARQUET, COMPRESSION 'ZSTD', COMPRESSION_LEVEL ${COMPRESSION});
"

## Execute the SQL command
echo -e "\nExecuting DuckDB command"

duckdb -c "${SQL_COMMAND}"

