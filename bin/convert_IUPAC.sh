#!/bin/bash

## Function to convert IUPAC codes in primers
# Based on PipeCraft2 scripts
# https://github.com/SuvalineVana/pipecraft/blob/main/src/pipecraft-core/service_scripts/submodules/framework.functions.sh
# Git commit 5650545 (Jun 9, 2022)
# Author - Sten Anslan

echo "$1" | \
if grep -q -E "R|Y|S|W|K|M|B|D|H|V|N|I" ; then
    
    ## Define IUPAC codes
    R=$"[AG]"
    Y=$"[CT]"
    S=$"[GC]"
    W=$"[AT]"
    K=$"[GT]"
    M=$"[AC]"
    B=$"[CGT]"
    D=$"[AGT]"
    H=$"[ACT]"
    V=$"[ACG]"
    N=$"[ATGC]"
    I=$"[ATGC]"
    
    ## Replace IUPAC codes
    primer=$(echo "$1" | \
    sed -e "s/R/$R/g; s/Y/$Y/g; \
    s/S/$S/g; s/W/$W/g; s/K/$K/g; \
    s/M/$M/g; s/B/$B/g; s/D/$D/g; \
    s/H/$H/g; s/V/$V/g; s/N/$N/g; \
    s/I/$I/g")
    
    ## Return convered primer
    echo "$primer"
else
    ## Return original primer when no IUPAC codes were detected
    echo "$1"
fi

## Example:
# ./convert_IUPAC.sh "CGACCWGCGGARGGATCATTA"  # CGACC[AT]GCGGA[AG]GGATCATTA
