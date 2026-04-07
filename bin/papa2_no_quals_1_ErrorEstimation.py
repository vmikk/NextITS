#!/usr/bin/env python3
"""
Learn error rates (no-quality model) from NextITS-style dereplicated FASTA/FASTQ.

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


def _parse_nonnegative_int(s: str) -> int:
    value = int(s)
    if value < 0:
        raise argparse.ArgumentTypeError(
            f"Expected a non-negative integer, got {s!r}"
        )
    return value


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
    p.add_argument(
        "-n",
        "--nbases",
        type=float,
        default=1e8,
        help="Target bases for error learning (default: 1e8)",
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
        "--maxreadsperseq",
        type=_parse_nonnegative_int,
        default=1000,
        help=(
            "Maximum effective reads per unique sequence during error learning "
            "(0 disables capping; default: 1000)"
        ),
    )
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
    print(f"Number of bases to use for error rate learning: {args.nbases}")
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
    print(
        "Maximum effective reads per unique for error learning: "
        f"{args.maxreadsperseq}"
    )
    print(f"Number of CPU threads to use: {args.threads}")
    print()

    os.environ["OMP_NUM_THREADS"] = str(max(1, args.threads))
    # Single-sample learning: avoid extra process pools
    os.environ.setdefault("DADA2_WORKERS", "1")

    random.seed(111)
    np.random.seed(111)

    ## Import after OMP_NUM_THREADS
    import papa2
    from papa2.dada import dada

    print("Loading papa2", papa2.__version__)

    print("\nLoading input data")
    derep_full, meta = papa2_io.load_nextits_derep(args.input)

    num_seqs = meta["num_seqs"]
    num_singl = meta["num_singl"]
    num_reads = meta["num_reads"]
    perc_ns = meta["perc_nonsingleton"]

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

    print("\nPreparing derep-class object (papa2 dict)")
    derep_learn = papa2_io.subsample_derep_by_nbases(
        derep_full,
        float(args.nbases),
        max_reads_per_seq=int(args.maxreadsperseq),
    )
    learn_reads = int(np.asarray(derep_learn["abundances"], dtype=np.int64).sum())
    bases_used = sum(
        int(derep_learn["abundances"][i]) * len(derep_learn["seqs"][i])
        for i in range(len(derep_learn["seqs"]))
    )
    print(
        f"Learning subset: {len(derep_learn['seqs'])} uniques, "
        f"{learn_reads} effective reads, "
        f"~{bases_used} effective bases (target nbases={args.nbases})"
    )

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

    print("\nEstimating error rates (self-consistency, noqual_errfun)")
    result = dada(
        derep_learn,
        err=None,
        error_estimation_function=papa2.noqual_errfun,
        self_consist=True,
        verbose=False,
        **dada_kw,
    )

    err = np.asarray(result["err_out"], dtype=np.float64)
    print(f"\nLearned error matrix shape: {err.shape} (16 x nQual)")
    print(f"Error rate min/max: {err.min():.4e} / {err.max():.4e}")

    ## Output path is relative to process cwd
    out_npz = Path.cwd() / "DADA2_ErrorRates_noqualErrfun.npz"
    print(f"\nExporting error rates to {out_npz}")
    np.savez_compressed(
        str(out_npz),
        err=err,
        nbases=np.array([args.nbases]),
        maxreadsperseq=np.array([args.maxreadsperseq], dtype=np.int64),
        input_path=np.array([args.input], dtype=object),
    )

    elapsed = (time.time() - start) / 60.0
    print(f"\nElapsed time: {elapsed:.4f} minutes")
    print("\nAll done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
