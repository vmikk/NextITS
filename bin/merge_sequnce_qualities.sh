#!/bin/bash

# $1 = input file
# $2 = text to add to the resulting file

zcat "$1" \
  | awk \
    -F '\t' -v OFS='\t' \
    -v fnm="$2" \
    '{ print $2 , $3 , $4 , fnm }' \
  | sed 's/_hash_table.txt//'
