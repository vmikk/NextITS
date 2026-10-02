#!/usr/bin/env python3
"""
Prepare tag (barcode) FASTA files for cutadapt-based demultiplexing of Illumina paired-end reads

Tags are searched at the 5' end of the mates, within a restricted window:
each tag is prefixed with `XN{window}` (non-internal 5' adapter, preceded by up to `window` bases)
If `--window 0`, tags are anchored at the read start (`^TAG`)

Modes:
- `single`          one tag per sample, on one of the mates
- `dual_symmetric`  one tag per sample, the same tag on both mates
- `dual_asymmetric` two tags per sample (`--fwd` and `--rev` files, the same sample order);
                    `--rev-orient revcomp` reverse-complements the reverse tags

Outputs (in `--outdir`):
- `T1.fasta`    tags searched at the 5' end of the first mate (named by sample)
- `T2.fasta`    tags searched at the 5' end of the second mate (the same order and names as T1; dual modes)
- `Tall.fasta`  all distinct tags (tag-jump detection; dual modes)
- `Tresc.fasta` tags used by a single sample only, named `<sample>~1` / `<sample>~2`
                (rescue of pairs with one readable tag; dual modes; may be empty)
"""

from __future__ import annotations
import argparse
import sys
from pathlib import Path

COMP = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn", "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")


def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def read_fasta(path: str):
    records = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                records.append([line[1:].split()[0], ""])
            else:
                if not records:
                    sys.exit(f"ERROR: malformed FASTA file {path}")
                records[-1][1] += line.upper()
    return [(name, seq) for name, seq in records]


def write_fasta(path: Path, records, prefix: str):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{prefix}{seq}\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", required=True, choices=["single", "dual_symmetric", "dual_asymmetric"])
    ap.add_argument("--tags", help="FASTA with one tag per sample (single, dual_symmetric)")
    ap.add_argument("--fwd", help="FASTA with forward tags (dual_asymmetric)")
    ap.add_argument("--rev", help="FASTA with reverse tags (dual_asymmetric)")
    ap.add_argument("--rev-orient", default="forward", choices=["forward", "revcomp"])
    ap.add_argument("--window", type=int, default=30)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args(argv)

    prefix = f"XN{{{a.window}}}" if a.window > 0 else "^"
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    if a.mode in ("single", "dual_symmetric"):
        if not a.tags:
            sys.exit("ERROR: --tags is required for single and dual_symmetric modes")
        tags = read_fasta(a.tags)
        fwd, rev = tags, tags
    else:
        if not (a.fwd and a.rev):
            sys.exit("ERROR: --fwd and --rev are required for dual_asymmetric mode")
        fwd = read_fasta(a.fwd)
        rev = read_fasta(a.rev)
        if [n for n, _ in fwd] != [n for n, _ in rev]:
            sys.exit("ERROR: forward and reverse tag files must list the same samples in the same order")
        if a.rev_orient == "revcomp":
            rev = [(n, revcomp(s)) for n, s in rev]

    if not fwd:
        sys.exit("ERROR: no tags found")

    write_fasta(out / "T1.fasta", fwd, prefix)

    if a.mode == "single":
        return 0

    write_fasta(out / "T2.fasta", rev, prefix)

    ## All distinct tags
    seen = {}
    for _, s in fwd + rev:
        seen.setdefault(s, f"tag{len(seen) + 1:04d}")
    write_fasta(out / "Tall.fasta", [(n, s) for s, n in seen.items()], prefix)

    ## Tags that identify a sample on their own
    users = {}
    for (n, s1), (_, s2) in zip(fwd, rev):
        users.setdefault(s1, set()).add(n)
        users.setdefault(s2, set()).add(n)

    resc = []
    for (n, s1), (_, s2) in zip(fwd, rev):
        if len(users[s1]) == 1:
            resc.append((f"{n}~1", s1))
        if s2 != s1 and len(users[s2]) == 1:
            resc.append((f"{n}~2", s2))
    write_fasta(out / "Tresc.fasta", resc, prefix)

    print(f"Samples: {len(fwd)}; distinct tags: {len(seen)}; tags usable for rescue: {len(resc)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
