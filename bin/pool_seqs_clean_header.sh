#!/bin/bash

# $1 = input file (full path)
# $2 = sample ID (e.g., basename of the input file)

zcat "${1}" \
  | sed -r '/^>/ s/;sample=[^;]*/;/g ; s/;;/;/g' \
  | sed "s/>.*/&;sample=${2}; / ; s/_NoChimera.fa//g ; s/_RescuedChimera.fa//g ; s/_JoinedPE//g ; s/Rescued_Chimeric_sequences.part_//g" \
  | sed -r '/^>/ s/;;/;/g'


