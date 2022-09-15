#!/bin/bash

# $1 = input file
# $2 = text to add to the resulting file

zcat "$1" \
  | awk \
    -F '\t' -v OFS='\t' \
    -v fnm="$2" \
    'NR>1 { print fnm , $2 , $3 , $4 , $5 , $6 }' \
  | sed 's/_hash_table.txt//' \
  | sed '1i SampleID\tSeqID\tSeqLen\tPhredScore\tMaxEE\tMEEP'
