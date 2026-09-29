#!/bin/bash

## Merge per-chunk ITSx / ITSx2 summary reports into a single per-sample report
##
## The pipeline runs the ITS extractor on chunks of the dereplicated sequences, so every chunk produces its own `*.summary.txt`
## Simply concatenating them would leave several copies of each counter in the file,
## and any downstream parser (e.g. the Step-1 report) would pick up the counts of the first chunk only
## Instead, the numeric counters are summed across chunks, while the layout of the report is taken from the first chunk
##
## Both summary formats are supported (any line of the form "<key>:<whitespace><integer>" is treated as a counter):
##   ITSx  - "Number of sequences in input file:  51", "  Fungi:  51", ...
##   ITSx2 - "Sequences processed: 51", ...
##
## Usage:
##   merge_itsx_summaries.sh -o Sample.summary.txt chunk0.summary.txt chunk1.summary.txt ...

set -euo pipefail

usage() {
    echo "Usage: $0 -o <output file> <summary file> [<summary file> ...]"
    echo "  -o : Output file (required)"
    exit 1
}

output_file=""

while getopts "o:h" opt; do
    case $opt in
        o) output_file="$OPTARG" ;;
        h) usage ;;
        ?) usage ;;
    esac
done
shift $((OPTIND - 1))

if [ -z "${output_file}" ]; then
    echo "Error: Output file is required"
    usage
fi

if [ "$#" -lt 1 ]; then
    echo "Error: At least one input summary file is required"
    usage
fi

for f in "$@"; do
    if [ ! -f "$f" ]; then
        echo "Error: File ${f} not found"
        exit 1
    fi
done

## A single chunk needs no merging
if [ "$#" -eq 1 ]; then
    cp -- "$1" "${output_file}"
    exit 0
fi

echo "..Merging $# summary reports into ${output_file}"

awk -v first="$1" '
  ## Split a counter line into its key and its value.
  ## The value is the trailing integer, the key is everything before the last colon.
  function counter_key(line,   i, ci) {
    if (line !~ /^.*:[ \t]*[0-9]+[ \t]*$/) return ""
    ci = 0
    for (i = length(line); i > 0; i--) {
      if (substr(line, i, 1) == ":") { ci = i; break }
    }
    if (ci == 0) return ""
    return substr(line, 1, ci - 1)
  }
  function counter_value(line,   i, ci, v) {
    ci = 0
    for (i = length(line); i > 0; i--) {
      if (substr(line, i, 1) == ":") { ci = i; break }
    }
    v = substr(line, ci + 1)
    gsub(/[ \t]/, "", v)
    return v + 0
  }

  {
    key = counter_key($0)
    if (key != "") sum[key] += counter_value($0)

    ## Keep the layout of the first report as the template
    if (FILENAME == first) template[++n] = $0

    ## The run end time should come from the chunk that finished last
    if ($0 ~ /^ITSx run finished at /) finished = $0
  }

  END {
    for (i = 1; i <= n; i++) {
      line = template[i]
      key  = counter_key(line)
      if (key != "") {
        printf "%s:\t%d\n", key, sum[key]
      } else if (line ~ /^ITSx run finished at / && finished != "") {
        print finished
      } else {
        print line
      }
    }
  }
' "$@" > "${output_file}"

echo "..Done"
