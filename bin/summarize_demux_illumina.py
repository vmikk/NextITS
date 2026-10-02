#!/usr/bin/env python3
"""
Summarize cutadapt-based demultiplexing of Illumina paired-end reads

Inputs:
- `seqkit stats --tabular` table for R1 files in `Strict/`, `StrictB/`, `Rescued/` and `Demux/`
  (`Rescued/<sample>~1_R1.fq.gz` and `<sample>~2` are pooled per sample)
- list of sample names
- cutadapt JSON reports of the demultiplexing passes

Outputs:
- per-sample table: `SampleID, Strict_Pairs, Rescued_Pairs, Demultiplexed_Pairs`
- run totals: input pairs, strict, rescued, tag-jumped (or unknown tag combinations), too short, unassigned

Exits with an error if the combined number of pairs differs from strict + rescued
"""

from __future__ import annotations
import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path


def load_json(path: str):
    if not path:
        return None
    return json.loads(Path(path).read_text())


def too_short(rep) -> int:
    if rep is None:
        return 0
    filt = rep["read_counts"].get("filtered") or {}
    return int(filt.get("too_short") or 0)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--pass1", nargs="+", required=True, help="JSON reports of strict demultiplexing (pass 1 and 1b)")
    ap.add_argument("--pass2", default="", help="JSON report of tag-jump removal")
    ap.add_argument("--pass3", default="", help="JSON report of the rescue pass")
    ap.add_argument("--mode", required=True)
    ap.add_argument("--revtag-orient", default="NA")
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--out-totals", required=True)
    a = ap.parse_args(argv)

    samples = [s.strip() for s in Path(a.samples).read_text().split("\n") if s.strip()]

    ## Per-sample counts from the files
    counts = defaultdict(lambda: defaultdict(int))
    with open(a.counts) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            path = Path(row["file"])
            name = re.sub(r"_R1\.fq\.gz$", "", path.name)
            name = re.sub(r"~[12]$", "", name)
            counts[path.parent.name][name] += int(row["num_seqs"].replace(",", ""))

    rows = []
    bad = []
    for s in samples:
        strict = counts["Strict"][s] + counts["StrictB"][s]
        resc = counts["Rescued"][s]
        total = counts["Demux"][s]
        if total != strict + resc:
            bad.append((s, strict, resc, total))
        rows.append((s, strict, resc, total))

    with open(a.out_summary, "w") as fh:
        fh.write("SampleID\tStrict_Pairs\tRescued_Pairs\tDemultiplexed_Pairs\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    if bad:
        for s, st, rs, tt in bad:
            print(f"MISMATCH: {s} strict={st} rescued={rs} combined={tt}", file=sys.stderr)
        sys.exit("ERROR: combined number of read pairs differs from strict + rescued")

    ## Run totals
    p1 = [load_json(p) for p in a.pass1 if p]
    p2 = load_json(a.pass2)
    p3 = load_json(a.pass3)

    n_input = int(p1[0]["read_counts"]["input"])
    n_strict = sum(r[1] for r in rows)
    n_resc = sum(r[2] for r in rows)
    n_short = sum(too_short(r) for r in p1) + too_short(p3)
    if p2 is not None:
        n_jump = int(p2["read_counts"]["input"]) - int(p2["read_counts"]["output"])
    else:
        n_jump = 0
    n_unas = n_input - n_strict - n_resc - n_jump - n_short

    def pct(x):
        return f"{100 * x / n_input:.2f}" if n_input else "0.00"

    with open(a.out_totals, "w") as fh:
        fh.write("Metric\tValue\tPercent\n")
        fh.write(f"Tag_layout\t{a.mode}\t\n")
        fh.write(f"Reverse_tag_orientation\t{a.revtag_orient}\t\n")
        fh.write(f"Input_Pairs\t{n_input}\t100.00\n")
        fh.write(f"Strict_Pairs\t{n_strict}\t{pct(n_strict)}\n")
        fh.write(f"Rescued_Pairs\t{n_resc}\t{pct(n_resc)}\n")
        fh.write(f"Demultiplexed_Pairs\t{n_strict + n_resc}\t{pct(n_strict + n_resc)}\n")
        fh.write(f"TagJump_Pairs\t{n_jump}\t{pct(n_jump)}\n")
        fh.write(f"TooShort_Pairs\t{n_short}\t{pct(n_short)}\n")
        fh.write(f"Unassigned_Pairs\t{n_unas}\t{pct(n_unas)}\n")
        fh.write(f"Samples_With_Reads\t{sum(1 for r in rows if r[3] > 0)}\t\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
