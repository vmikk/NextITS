#!/bin/bash

# $1 = input file
# $2 = text to add to the resulting file

zcat "$1" \
  | awk \
    -F '\t' -v OFS='\t' \
    -v fnm="$2" \
    '$1 ~ /H/ { print fnm , $9 , $10 }' \
  | sed 's/_uch.uc//'

