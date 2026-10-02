#!/usr/bin/env python3
"""
Add sequences excluded from DADA2 denoising (low-abundance or with ambiguities)
to the DADA2 results, based on their matches to ASVs (VSEARCH `--usearch_global`)

Inputs:
- ASVs (`DADA2_denoised.fa.gz`), headers `SeqID;size=ABUNDANCE`
- Pseudo-UC file produced by DADA2 or PAPA2 (`DADA2_denoised.uc.gz`)
- Sequences excluded from denoising, headers `SeqID;size=ABUNDANCE`
- UC file with matches of excluded sequences to ASVs (`--usearch_global`)

Outputs:
- FASTA with ASVs (abundances include the mapped sequences), sorted by abundance;
  unmapped sequences are appended as separate sequences if `--unmapped keep`
- Pseudo-UC file (DADA2 records + `H` records for mapped sequences + `S` records for unmapped sequences if `--unmapped keep`)
- Summary lines appended to the DADA2 summary file
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path

_SIZE_RE = re.compile(r";size=[0-9]+;?")


def _open_text(path: str, mode: str = "rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="ascii", newline="\n")
    return open(path, mode, encoding="ascii", newline="\n")


def _is_empty(path: str) -> bool:
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return True
    with _open_text(path) as fh:
        for line in fh:
            if line.strip():
                return False
    return True


def _strip_size(label: str) -> str:
    return _SIZE_RE.sub("", label).rstrip(";")


def _parse_size(label: str) -> int:
    m = re.search(r";size=([0-9]+)", label)
    if m is None:
        raise ValueError(f"Expected ';size=' in sequence header, got: {label[:120]!r}")
    return int(m.group(1))


def read_fasta(path: str):
    """Yield (SeqID, abundance, sequence) tuples."""
    if _is_empty(path):
        return
    with _open_text(path) as fh:
        header = None
        chunks = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield _strip_size(header), _parse_size(header), "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield _strip_size(header), _parse_size(header), "".join(chunks)


def read_mapback_uc(path: str) -> dict:
    """Parse VSEARCH UC (usearch_global), return {query_id: (target_id, strand)}."""
    hits = {}
    if _is_empty(path):
        return hits
    with _open_text(path) as fh:
        for line in fh:
            if not line.startswith("H\t"):
                continue
            cols = line.rstrip("\n").split("\t")
            query = _strip_size(cols[8])
            target = _strip_size(cols[9])
            if query not in hits:
                hits[query] = (target, cols[4] or "+")
    return hits


def _uc_line(record_type: str, strand: str, query: str, target: str) -> str:
    # RecordType, ClustNum, SeqLen, Ident, Strand, V6, V7, ALN, DerepSeqID, ASV
    return "\t".join([record_type, "", "", "", strand, "", "", ".", query, target])


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--asvs", required=True, help="FASTA with ASVs")
    p.add_argument("--uc", required=True, help="Pseudo-UC file from DADA2/PAPA2")
    p.add_argument("--excluded", required=True, help="FASTA with sequences excluded from denoising")
    p.add_argument("--mapback", required=True, help="UC file with matches of excluded sequences to ASVs")
    p.add_argument("--unmapped", choices=["keep", "discard"], default="keep",
                   help="What to do with excluded sequences that do not match any ASV")
    p.add_argument("--summary", required=True, help="Summary file (lines are appended)")
    p.add_argument("--out_fasta", required=True, help="Output FASTA (uncompressed)")
    p.add_argument("--out_uc", required=True, help="Output pseudo-UC file (uncompressed)")
    args = p.parse_args(argv)

    print("\nAdding sequences excluded from denoising to the DADA2 results")

    ## ASVs
    asv_seqs = {}
    asv_abund = {}
    for sid, abund, seq in read_fasta(args.asvs):
        asv_seqs[sid] = seq
        asv_abund[sid] = abund
    print(f"..Number of ASVs: {len(asv_seqs)}")

    ## Matches
    hits = read_mapback_uc(args.mapback)
    unknown = {t for t, _ in hits.values() if t not in asv_abund}
    if unknown:
        print(f"ERROR: {len(unknown)} match targets are not in the ASV file "
              f"(e.g., {sorted(unknown)[0]})", file=sys.stderr)
        return 1

    ## Existing pseudo-UC records
    uc_lines = []
    if not _is_empty(args.uc):
        with _open_text(args.uc) as fh:
            uc_lines = [line.rstrip("\n") for line in fh if line.strip()]

    n_excl = n_excl_reads = 0
    n_mapped = n_mapped_reads = 0
    n_unmapped = n_unmapped_reads = 0
    kept = {}
    for sid, abund, seq in read_fasta(args.excluded):
        n_excl += 1
        n_excl_reads += abund
        hit = hits.get(sid)
        if hit is not None:
            target, strand = hit
            asv_abund[target] += abund
            uc_lines.append(_uc_line("H", strand, sid, target))
            n_mapped += 1
            n_mapped_reads += abund
        else:
            n_unmapped += 1
            n_unmapped_reads += abund
            if args.unmapped == "keep":
                if sid in asv_seqs or sid in kept:
                    print(f"ERROR: duplicated sequence ID {sid}", file=sys.stderr)
                    return 1
                kept[sid] = (abund, seq)
                uc_lines.append(_uc_line("S", "+", sid, sid))

    for sid, (abund, seq) in kept.items():
        asv_seqs[sid] = seq
        asv_abund[sid] = abund

    ## Export ASVs (sorted by abundance, then by SeqID)
    order = sorted(asv_seqs, key=lambda s: (-asv_abund[s], s))
    with open(args.out_fasta, "w", encoding="ascii", newline="\n") as out:
        for sid in order:
            out.write(f">{sid};size={asv_abund[sid]}\n{asv_seqs[sid]}\n")

    with open(args.out_uc, "w", encoding="ascii", newline="\n") as out:
        for line in uc_lines:
            out.write(line + "\n")

    def perc(x, total):
        return round(x / total * 100, 2) if total else 0.0

    rows = [
        ("Number of sequences excluded from denoising", n_excl),
        ("Number of reads of sequences excluded from denoising", n_excl_reads),
        ("Number of excluded sequences mapped to ASVs", n_mapped),
        ("Number of reads of excluded sequences mapped to ASVs", n_mapped_reads),
        (f"Number of excluded sequences not mapped to ASVs ({args.unmapped})", n_unmapped),
        (f"Number of reads of excluded sequences not mapped to ASVs ({args.unmapped})", n_unmapped_reads),
        ("Percentage of excluded reads mapped to ASVs", perc(n_mapped_reads, n_excl_reads)),
        ("Number of sequences in the output", len(asv_seqs)),
    ]
    with open(args.summary, "a", encoding="utf-8") as fh:
        for name, val in rows:
            fh.write(f"{name}\t{val}\n")

    for name, val in rows:
        print(f"..{name}: {val}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
