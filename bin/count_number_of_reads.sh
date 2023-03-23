#!/bin/bash

## Count number of reads in the dereplicated file
## Size annotations should be in USEARCH-style (e.g., size=100)

# $1 = input file
# $2 = text to add to the resulting file

seqkit seq --name "$1" \
  | grep -Po ';size=[0-9]+' \
  | sed 's/;size=//g' \
  | awk -F '\t' -v OFS='\t' -v fnm="$2" '{sum+=$1} END {print fnm , sum}'
