#!/usr/bin/env python3
"""
Learn error rates (no-quality model) from NextITS-style dereplicated FASTA/FASTQ.

The whole input is used for learning (as R `learnErrors` does with a single sample);
the input size is controlled upstream (e.g., by selecting whole buckets).

Outputs:
- DADA2_ErrorRates_noqualErrfun.npz (NumPy array, compressed)
"""

from __future__ import annotations

import argparse
import os
import random
import sys
import time
from pathlib import Path

## Load co-located helper module (should be in `bin/` on PATH, as other scripts)
_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

import numpy as np
import papa2_io


def _configure_stdio() -> None:
    """Force line-buffered console output for batch/HPC log files."""
    for stream_name in ("stdout", "stderr"):
        stream = getattr(sys, stream_name, None)
        reconfigure = getattr(stream, "reconfigure", None)
        if callable(reconfigure):
            reconfigure(line_buffering=True, write_through=True)


def _parse_bool(s: str) -> bool:
    x = str(s).strip().upper()
    if x in ("TRUE", "T", "1", "YES"):
        return True
    if x in ("FALSE", "F", "0", "NO"):
        return False
    raise argparse.ArgumentTypeError(f"Invalid boolean: {s!r}")


def _parse_hpgap(s: str | None):
    if s is None or str(s).strip() == "":
        return None
    t = str(s).strip().upper()
    if t in ("NA", "NULL", "NONE"):
        return None
    return float(s)


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Learn DADA2 error rates (noqual_errfun) via papa2."
    )
    p.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input dereplicated FASTA/FASTQ (gzip-compressed data supported)",
    )
    p.add_argument("-b", "--bandsize", type=float, default=16.0)
    p.add_argument(
        "-s",
        "--detectsingletons",
        type=_parse_bool,
        default=True,
        help="TRUE/FALSE (default: TRUE)",
    )
    p.add_argument("-A", "--omegaA", type=float, default=1e-20)
    p.add_argument("-C", "--omegaC", type=float, default=1e-40)
    p.add_argument("-P", "--omegaP", type=float, default=1e-4)
    p.add_argument("-x", "--maxconsist", type=int, default=10)
    p.add_argument("--match", type=float, default=4.0)
    p.add_argument("--mismatch", type=float, default=-5.0)
    p.add_argument("--gappenalty", type=float, default=-8.0)
    p.add_argument("--hpgap", type=_parse_hpgap, default=None)
    p.add_argument(
        "-t",
        "--threads",
        type=int,
        default=4,
        help="OMP threads for C core (set before papa2 import at runtime)",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    _configure_stdio()
    args = _build_parser().parse_args(argv)
    start = time.time()

    print("\nParsing input options and arguments...\n")

    if not Path(args.input).exists():
        print(f"Input file not found: {args.input}", file=sys.stderr)
        return 1

    print("Parameters specified:")
    print(f"Input file: {args.input}")
    print(f"Band size for the Needleman-Wunsch alignment: {args.bandsize}")
    print(f"Singleton detection: {args.detectsingletons}")
    print(f"OMEGA_A: {args.omegaA}")
    print(f"OMEGA_C: {args.omegaC}")
    print(f"OMEGA_P: {args.omegaP}")
    print(f"Number of iterations of the self-consistency loop: {args.maxconsist}")
    print(f"Alignment for matches: {args.match}")
    print(f"Alignment for mismatches: {args.mismatch}")
    print(f"Gap penalty: {args.gappenalty}")
    print(f"Homopolymer gap penalty: {args.hpgap}")
    print(f"Number of CPU threads to use: {args.threads}")
    print()

    ## papa2 sizes its OpenMP pool from DADA2_CORES / DADA2_OMP_THREADS
    ## (falls back to os.cpu_count(), ignoring OMP_NUM_THREADS)
    nthreads = str(max(1, args.threads))
    os.environ["OMP_NUM_THREADS"] = nthreads
    os.environ["DADA2_CORES"] = nthreads

    random.seed(111)
    np.random.seed(111)

    print("\nLoading input data")
    ## Sequences may be split into groups (`;grp=` in headers, groups of whole sequence clusters);
    ## each group is processed as a separate sample (transition counts are pooled across groups),
    ## which avoids comparing unrelated sequences and reduces the computation time
    groups = papa2_io.load_nextits_derep_groups(args.input)
    dereps = [d for _, d, _ in groups]
    metas = [m for _, _, m in groups]

    if len(dereps) == 1:
        # Single sample: use all threads within the sample, avoid extra process pools
        os.environ["DADA2_OMP_THREADS"] = nthreads
        os.environ["DADA2_WORKERS"] = "1"
    # Otherwise, papa2 splits the threads between parallel samples and within-sample threads

    ## Import after setting the number of threads
    import papa2
    from papa2.dada import dada

    print("Loading papa2", papa2.__version__)

    num_seqs = sum(m["num_seqs"] for m in metas)
    num_singl = sum(m["num_singl"] for m in metas)
    num_reads = sum(m["num_reads"] for m in metas)
    perc_ns = round((num_seqs - num_singl) / num_seqs * 100, 2) if num_seqs else 0.0

    print("\n")
    print(f"Number of unique sequences detected: {num_seqs}")
    print(f"Number of singleton sequences: {num_singl}")
    print(f"Total abundance of sequences: {num_reads}")
    print(f"Percentage of non-singleton sequences: {perc_ns}")

    if perc_ns < 10:
        print(
            "WARNING: <10% of reads are duplicates of other reads,\n"
            "         meaning that DADA2 might not be the right algorithmic choice"
        )

    bases_used = sum(
        int(d["abundances"][i]) * len(d["seqs"][i])
        for d in dereps
        for i in range(len(d["seqs"]))
    )
    print(f"Total read bases used for error learning: {bases_used}")
    print(f"Number of sequence groups (processed as separate samples): {len(dereps)}")
    if len(dereps) > 1:
        gsizes = sorted(len(d["seqs"]) for d in dereps)
        print(f"Unique sequences per group: min {gsizes[0]}, "
              f"median {gsizes[len(gsizes) // 2]}, max {gsizes[-1]}")

    # NW / gap scores are integers in the C API (ctypes c_int); argparse gives float.
    hpgap = args.hpgap
    dada_kw = dict(
        BAND_SIZE=int(args.bandsize),
        DETECT_SINGLETONS=bool(args.detectsingletons),
        OMEGA_A=args.omegaA,
        # R learnErrors forces OMEGA_C=0 during learning; match papa2 learn_errors
        OMEGA_C=0.0,
        OMEGA_P=args.omegaP,
        MAX_CONSIST=int(args.maxconsist),
        MATCH=int(args.match),
        MISMATCH=int(args.mismatch),
        GAP_PENALTY=int(args.gappenalty),
        HOMOPOLYMER_GAP_PENALTY=None if hpgap is None else int(hpgap),
        USE_QUALS=False,
    )

    ## papa2 calls the error function once per self-consistency round, so wrap it to report progress
    round_start = time.time()
    round_num = 0

    def errfun_with_progress(trans):
        nonlocal round_start, round_num
        err = papa2_io.noqual_errfun_pc(trans)
        subst = np.delete(err[:, 0], [0, 5, 10, 15])
        print(
            f"..Round {round_num}: {(time.time() - round_start) / 60.0:.2f} min, "
            f"substitutions observed: {papa2_io.count_substitutions(trans)}, "
            f"mean substitution rate: {subst.mean():.4e}"
        )
        round_num += 1
        round_start = time.time()
        return err

    print("\nEstimating error rates (self-consistency, noqual_errfun)")
    print("(round 0 is the initialization pass with a single cluster)")
    results = dada(
        dereps if len(dereps) > 1 else dereps[0],
        err=None,
        error_estimation_function=errfun_with_progress,
        self_consist=True,
        verbose=False,
        **dada_kw,
    )

    if isinstance(results, dict):
        results = [results]
    err = np.asarray(results[0]["err_out"], dtype=np.float64)
    ncol = max(np.asarray(r["trans"]).shape[1] for r in results)
    trans = np.zeros((16, ncol), dtype=np.float64)
    for r in results:
        t = np.asarray(r["trans"], dtype=np.float64)
        trans[:, : t.shape[1]] += t
    nsubs = papa2_io.count_substitutions(trans)
    print(f"\nSelf-consistency rounds: {len(results[0]['err_in'])}")
    print(f"Number of ASVs in the last round: {sum(len(r['cluster_seqs']) for r in results)}")
    print(f"Observed substitutions (read-weighted): {nsubs}")
    print(f"Observed transitions (read-weighted): {int(round(trans.sum()))}")
    print(f"Learned error matrix shape: {err.shape} (16 x nQual)")
    print("Substitution rates:")
    print("\n".join(papa2_io.format_error_rates(err)))

    ## With no observed substitutions the model only contains the pseudocount, and would mark nearly every variant as a new ASV
    if nsubs == 0:
        print(
            "\nERROR: no substitutions were observed during error learning, "
            "the error model is degenerate.\n"
            "The learning input probably contains no error variants "
            "(e.g., only the most abundant sequences).",
            file=sys.stderr,
        )
        return 1

    ## Output path is relative to process cwd
    out_npz = Path.cwd() / "DADA2_ErrorRates_noqualErrfun.npz"
    print(f"\nExporting error rates to {out_npz}")
    np.savez_compressed(
        str(out_npz),
        err=err,
        trans=trans,
        input_path=np.array([args.input], dtype=object),
    )

    elapsed = (time.time() - start) / 60.0
    print(f"\nElapsed time: {elapsed:.4f} minutes")
    print("\nAll done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
