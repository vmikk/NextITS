#!/bin/bash

## Import sequences into DuckDB-compatible tables

## Input = FASTA formatted sequences (header = "hash;size=...")
## Output = Table in DuckDB-native format or Parquet

# Define usage function
usage() {
    echo "Usage: $0 [-i input_file] [-o output_file] [-f format]"
    echo "  -i : Input FASTA file (required)"
    echo "  -o : Output file (optional, defaults to input filename with .db/.parquet extension)"
    echo "  -f : Output format (optional): 'duckdb' or 'parquet' (default, 'parquet')"
    exit 1
}


## Extract rRNA region name from filename
if [[ $input_file =~ ([^.]+)\.fasta\.gz$ ]]; then
    rRNA_part="${BASH_REMATCH[1]}"
else
    echo "Error in extracting rRNA region name from filename"
    rRNA_part="X"
fi

## Check if rRNA region name is valid
VALID_PARTS=("full" "SSU" "ITS1" "5_8S" "ITS2" "LSU")
if [[ ! " ${VALID_PARTS[@]} " =~ " ${rRNA_part} " ]]; then
    echo "..Error: Invalid rRNA region name. Supported names are: ${VALID_PARTS[*]}"
    rRNA_part="X"
fi

## 'full' is a reserved keyword in DuckDB, rename to ITS
if [ "$rRNA_part" == "full" ]; then
    rRNA_part="ITS"
fi

## DuckDB table name cannot start with a number
if [[ "$rRNA_part" == "5_8S" ]]; then
    rRNA_part="S58"
fi

## Extract sample name from filename
sample_name="${input_file/.fasta.gz/}"
db_file="${sample_name}.db"

## Check if input file exists
if [ ! -f "$input_file" ]; then
    echo "..Error: File $input_file not found"
    exit 1
fi

echo "..Importing $input_file into $db_file"

seqkit fx2tab "$input_file" \
  | sed 's/;size=/\t/'  \
  | duckdb "$db_file" \
"
DROP TABLE IF EXISTS ${rRNA_part};
CREATE TABLE ${rRNA_part} (
  SeqID VARCHAR PRIMARY KEY,
  Abundance INTEGER,
  Sequence VARCHAR
);

INSERT INTO ${rRNA_part}
SELECT * FROM read_csv(
  '/dev/stdin',
  header = false,
  columns = {
    'SeqID': 'VARCHAR',
    'Abundance': 'INTEGER',
    'Sequence': 'VARCHAR'
  }
);"

echo "..Data imported"


#### Check the data
# duckdb "$db_file"
# 
# -- Show all tables
# SHOW TABLES;
# SELECT * FROM information_schema.tables;
# 
# -- Show all column names and their types
# DESCRIBE ITS1;
# 
# -- Show first 10 rows
# SELECT * FROM ITS1 LIMIT 10;
# 
# -- Get count of rows
# SELECT COUNT(*) FROM ITS1;
#