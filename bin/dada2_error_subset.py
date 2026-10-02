#!/usr/bin/env python3
"""
Prepare a subset of dereplicated sequences for DADA2 error-rate learning

Buckets are processed in random order (fixed seed) until the number of read bases
(abundance x length, as `nbases` in DADA2 `learnErrors`) reaches the target.
Within each bucket:
  1. If the bucket has more than `nbases / minbuckets` read bases, all its reads are
     subsampled with the same probability, so that no bucket
     (e.g., a bucket built around a super-abundant sequence)
     contributes more than 1/minbuckets of the subset
  2. Sequences are clustered with VSEARCH (`--cluster_size`),
     so that each cluster contains a parent sequence together with its error variants
  3. Clusters are packed into groups of up to `--groupsize` sequences;
     the group label is added to sequence headers (`SeqID;size=N;grp=LABEL`),
     and each group is used as a separate sample during error learning
     (error variants are within the same cluster as their parent,
     so comparisons between unrelated sequences are avoided)

Subsampling is done at the level of reads (new abundance ~ Binomial(abundance, p)),
with the same probability for all members of a bucket. This keeps the ratio
between parent sequences and their error variants, which the error model is learned from.

Input:
- Bucket FASTA files, headers `SeqID;size=ABUNDANCE`

Outputs:
- FASTA with the subset (gzip-compressed), headers `SeqID;size=ABUNDANCE;grp=LABEL`
- TSV with per-bucket statistics
"""

from __future__ import annotations

import argparse
import gzip
import os
import random
import re
import subprocess
import sys
import tempfile
from collections import defaultdict

import numpy as np

_SIZE_RE = re.compile(r";size=([0-9]+)")


def read_fasta(path: str):
    """Return lists of SeqIDs, abundances, and sequences."""
    ids, sizes, seqs = [], [], []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        header, chunks = None, []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs.append("".join(chunks))
                header = line[1:]
                m = _SIZE_RE.search(header)
                if m is None:
                    raise ValueError(f"Expected ';size=' in header: {header[:120]!r}")
                ids.append(header.split(";size=")[0])
                sizes.append(int(m.group(1)))
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            seqs.append("".join(chunks))
    return ids, np.array(sizes, dtype=np.int64), seqs


def cluster_sequences(ids, sizes, seqs, cluster_id, threads, tmpdir):
    """Cluster sequences with VSEARCH, return a cluster index per sequence."""
    fasta = os.path.join(tmpdir, "seqs.fa")
    uc = os.path.join(tmpdir, "clusters.uc")
    with open(fasta, "w") as out:
        for i in np.argsort(-sizes, kind="stable"):
            out.write(f">{ids[i]};size={sizes[i]}\n{seqs[i]}\n")
    subprocess.run(
        [
            "vsearch", "--cluster_size", fasta,
            "--id", str(cluster_id),
            "--sizein",
            "--strand", "both",
            "--threads", str(threads),
            "--uc", uc,
            "--quiet",
        ],
        check=True,
    )
    index = {sid: i for i, sid in enumerate(ids)}
    clust = np.full(len(ids), -1, dtype=np.int64)
    with open(uc) as fh:
        for line in fh:
            if line[0] not in "SH":
                continue
            cols = line.split("\t")
            clust[index[cols[8].split(";size=")[0]]] = int(cols[1])
    os.remove(fasta)
    os.remove(uc)
    if np.any(clust < 0):
        raise RuntimeError("Some sequences were not assigned to clusters")
    return clust


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--buckets", nargs="+", required=True, help="Bucket FASTA files")
    p.add_argument("--nbases", type=float, default=5e8,
                   help="Target number of read bases (default: 5e8)")
    p.add_argument("--minbuckets", type=int, default=5,
                   help="Minimum number of buckets to spread the target over; each bucket contributes "
                        "at most nbases/minbuckets read bases (default: 5)")
    p.add_argument("--clusterid", type=float, default=0.97,
                   help="Sequence identity for clustering (default: 0.97)")
    p.add_argument("--groupsize", type=int, default=2000,
                   help="Maximum number of sequences per group of clusters (default: 2000)")
    p.add_argument("--seed", type=int, default=111)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--output", default="DADA2_error_subset.fa.gz")
    p.add_argument("--stats", default="DADA2_error_subset_buckets.tsv")
    args = p.parse_args(argv)
    sys.stdout.reconfigure(line_buffering=True)

    rng = np.random.default_rng(args.seed)
    buckets = sorted(args.buckets)
    random.Random(args.seed).shuffle(buckets)

    print(f"Target number of read bases: {args.nbases:.0f}")
    bucket_budget = args.nbases / max(1, args.minbuckets)
    print(f"Maximum number of read bases per bucket: {bucket_budget:.0f}")
    print(f"Clustering identity: {args.clusterid}")
    print(f"Number of buckets available: {len(buckets)}\n")

    total_bases = 0
    rows = []
    with gzip.open(args.output, "wt", compresslevel=6) as out, \
            tempfile.TemporaryDirectory(dir=".") as tmpdir:
        for bucket_num, path in enumerate(buckets, start=1):
            if total_bases >= args.nbases:
                break
            ids, sizes, seqs = read_fasta(path)
            if not ids:
                continue
            lens = np.array([len(s) for s in seqs], dtype=np.int64)
            in_seqs = len(ids)
            in_reads = int(sizes.sum())
            in_bases = int((sizes * lens).sum())

            ## 1. Limit the contribution of the bucket
            p0 = min(1.0, bucket_budget / in_bases)
            if p0 < 1.0:
                sizes = rng.binomial(sizes, p0)
                keep = sizes > 0
                ids = [x for x, k in zip(ids, keep) if k]
                seqs = [x for x, k in zip(seqs, keep) if k]
                sizes, lens = sizes[keep], lens[keep]

            ## 2. Cluster sequences
            if len(ids) > 1:
                clust = cluster_sequences(ids, sizes, seqs, args.clusterid, args.threads, tmpdir)
            else:
                clust = np.zeros(len(ids), dtype=np.int64)
            n_clusters = len(set(clust.tolist()))

            ## 3. Pack clusters into groups
            cseqs = defaultdict(int)
            for c, sz in zip(clust, sizes):
                if sz > 0:
                    cseqs[c] += 1
            group_of = {}
            g, gn = 1, 0
            for c in sorted(cseqs):
                if gn > 0 and gn + cseqs[c] > args.groupsize:
                    g, gn = g + 1, 0
                group_of[c] = g
                gn += cseqs[c]

            out_seqs = out_reads = out_bases = 0
            for sid, size, seq, ln, c in zip(ids, sizes, seqs, lens, clust):
                if size > 0:
                    out.write(f">{sid};size={size};grp={bucket_num}_{group_of[c]}\n{seq}\n")
                    out_seqs += 1
                    out_reads += int(size)
                    out_bases += int(size) * int(ln)
            total_bases += out_bases

            n_groups = len(set(group_of.values()))
            rows.append((os.path.basename(path), in_seqs, in_reads, in_bases,
                         round(p0, 6), n_clusters, n_groups,
                         out_seqs, out_reads, out_bases))
            print(f"..{os.path.basename(path)}: {in_reads} reads -> {out_reads} reads "
                  f"({out_seqs} sequences; {n_clusters} clusters; {n_groups} groups); "
                  f"total read bases: {total_bases}")

    with open(args.stats, "w") as fh:
        fh.write("\t".join(["File", "NumSeqs", "NumReads", "NumReadBases", "BucketSamplingProb",
                            "NumClusters", "NumGroups",
                            "OutSeqs", "OutReads", "OutReadBases"]) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    print(f"\nBuckets used: {len(rows)}")
    print(f"Groups of clusters: {sum(r[6] for r in rows)}")
    print(f"Sequences: {sum(r[7] for r in rows)}")
    print(f"Reads: {sum(r[8] for r in rows)}")
    print(f"Read bases: {total_bases}")
    if total_bases < args.nbases:
        print("NB! all buckets were used, but the target number of read bases was not reached")
    return 0


if __name__ == "__main__":
    sys.exit(main())
