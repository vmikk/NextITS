#!/bin/bash

## Function to reverse-complement DNA sequences (with the support of IUPAC codes)

echo "$1" \
  | tr \
    "[ATGCUatgcuNnYyRrSsWwKkMmBbDdHhVv]" \
    "[TACGAtacgaNnRrYySsWwMmKkVvHhDdBb]" \
  | rev
