#!/usr/bin/env python3
"""
Denoise NextITS-style dereplicated FASTA/FASTQ with a pre-learned error model

Inputs:
- NextITS-style dereplicated FASTA/FASTQ: headers `SeqID;size=ABUNDANCE`
- Pre-learned error model: `DADA2_ErrorRates_noqualErrfun.npz`

Outputs:
- `DADA2_InferedSeqs_noqualErrfun.pkl.gz`
- `DADA2_denoised.fa.gz`
- `DADA2_denoised.uc.gz`
"""

from __future__ import annotations

import argparse
import gzip
import os
import pickle
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
        description="DADA2 sample inference (noqual_errfun) via papa2."
    )
    p.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input dereplicated FASTA/FASTQ (gzip-compressed data supported)",
    )
    p.add_argument(
        "-e",
        "--errors",
        required=True,
        help="Error model (.npz from papa2_no_quals_1_ErrorEstimation.py)",
    )
    p.add_argument("-b", "--bandsize", type=float, default=16.0)
    p.add_argument(
        "-s",
        "--detectsingletons",
        type=_parse_bool,
        default=True,
    )
    p.add_argument("-A", "--omegaA", type=float, default=1e-20)
    p.add_argument("-C", "--omegaC", type=float, default=1e-40)
    p.add_argument("-P", "--omegaP", type=float, default=1e-4)
    p.add_argument("-x", "--maxconsist", type=int, default=10)
    p.add_argument("--match", type=float, default=4.0)
    p.add_argument("--mismatch", type=float, default=-5.0)
    p.add_argument("--gappenalty", type=float, default=-8.0)
    p.add_argument("--hpgap", type=_parse_hpgap, default=None)
    p.add_argument("-t", "--threads", type=int, default=4)
    return p


def _load_errors(path: str) -> np.ndarray:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(path)
    if p.suffix.lower() == ".npz":
        z = np.load(str(p), allow_pickle=True)
        if "err" not in z.files:
            raise KeyError(f"{path} must contain array 'err'")
        return np.asarray(z["err"], dtype=np.float64)
    raise ValueError(f"Unsupported error model format: {path} (expected .npz)")


def _extend_err_to_maxq(err: np.ndarray, max_q: int) -> np.ndarray:
    if max_q + 1 <= err.shape[1]:
        return err
    extra = max_q + 1 - err.shape[1]
    tail = np.tile(err[:, -1:], (1, extra))
    return np.hstack([err, tail])


def main(argv: list[str] | None = None) -> int:
    _configure_stdio()
    args = _build_parser().parse_args(argv)
    start = time.time()

    if not Path(args.input).exists():
        print(f"Input file not found: {args.input}", file=sys.stderr)
        return 1
    if not Path(args.errors).exists():
        print(f"Error model file not found: {args.errors}", file=sys.stderr)
        return 1

    print("\nParsing input options and arguments...\n")
    print("Parameters specified:")
    print(f"Input file: {args.input}")
    print(f"Pre-computed error model: {args.errors}")
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

    os.environ["OMP_NUM_THREADS"] = str(max(1, args.threads))
    os.environ.setdefault("DADA2_WORKERS", "1")

    random.seed(111)
    np.random.seed(111)

    import papa2
    from papa2.dada import dada

    print("Loading papa2", papa2.__version__)

    print("\nLoading input data")
    derep, meta = papa2_io.load_nextits_derep(args.input)
    num_seqs = meta["num_seqs"]
    num_singl = meta["num_singl"]
    num_reads = meta["num_reads"]
    perc_ns = meta["perc_nonsingleton"]
    seq_ids_file = meta["seq_ids_file_order"]
    abunds_file = meta["abundances_file_order"]
    uidx_file = meta["unique_index_file_order"]

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

    print("\nLoading pre-computed error rates")
    err = _load_errors(args.errors)
    max_q = 0
    if derep["quals"].size and not np.all(np.isnan(derep["quals"])):
        max_q = int(np.nanmax(derep["quals"]))
    err = _extend_err_to_maxq(err, max_q)
    print(f"Error matrix shape (after Q extension): {err.shape}")

    # NW / gap scores are integers in the C API (ctypes c_int); argparse gives float.
    hpgap = args.hpgap
    dada_kw = dict(
        BAND_SIZE=int(args.bandsize),
        DETECT_SINGLETONS=bool(args.detectsingletons),
        OMEGA_A=args.omegaA,
        OMEGA_C=args.omegaC,
        OMEGA_P=args.omegaP,
        MAX_CONSIST=int(args.maxconsist),
        MATCH=int(args.match),
        MISMATCH=int(args.mismatch),
        GAP_PENALTY=int(args.gappenalty),
        HOMOPOLYMER_GAP_PENALTY=None if hpgap is None else int(hpgap),
        USE_QUALS=False,
    )

    print("\nRunning sample inference")
    dadares = dada(
        derep,
        err=err,
        error_estimation_function=papa2.noqual_errfun,
        self_consist=False,
        verbose=False,
        **dada_kw,
    )

    out_pkl = Path("DADA2_InferedSeqs_noqualErrfun.pkl.gz")
    print(f"\nExporting DADA2 object to {out_pkl}")
    with gzip.open(out_pkl, "wb") as fh:
        pickle.dump(dadares, fh, protocol=pickle.HIGHEST_PROTOCOL)

    cluster_seqs = dadares["cluster_seqs"]
    cluster_abunds = np.asarray(dadares["cluster_abunds"])
    n_asv = len(cluster_seqs)

    ## First SeqID seen in file order for each sequence
    seq_to_first_sqid: dict = {}
    for sid, seq in zip(seq_ids_file, meta["sequences_file_order"]):
        if seq not in seq_to_first_sqid:
            seq_to_first_sqid[seq] = sid

    ## One row per ASV in dadares cluster order; SeqNumID == .I == 1..n (R)
    res_rows = []
    for i in range(n_asv):
        seq = cluster_seqs[i]
        res_rows.append(
            {
                "Sequence": seq,
                "Abundance": int(cluster_abunds[i]),
                "SeqNumID": i + 1,
                "SeqID": seq_to_first_sqid.get(seq, ""),
            }
        )
    seqnum_to_rep_seqid = {r["SeqNumID"]: r["SeqID"] for r in res_rows}

    ## Sort: -Abundance, SeqID, NA last (missing SeqID last)
    res_sorted = sorted(
        res_rows,
        key=lambda r: (
            -r["Abundance"],
            r["SeqID"] is None or r["SeqID"] == "",
            r["SeqID"] or "",
        ),
    )

    print("Preparing resulting table")
    print(
        f"ASVs inferred: {n_asv}, total reads in ASV table: "
        f"{int(cluster_abunds.sum())}"
    )

    ## Pseudo-UC: one row per input record (file order)
    map_u = np.asarray(dadares["map"], dtype=np.int64)
    uc_pre = []
    for k in range(len(seq_ids_file)):
        uidx = int(uidx_file[k])
        mid = int(map_u[uidx])
        if mid < 0:
            m1 = None
            rep_sid = None
        else:
            ## papa2 map: 0-based cluster index -> R SeqNumID is 1-based
            m1 = mid + 1
            rep_sid = seqnum_to_rep_seqid.get(m1, "")
        uc_pre.append(
            {
                "DerepSeqID": seq_ids_file[k],
                "SeqNumID": m1,
                "Abundance": int(abunds_file[k]),
                "ASV": rep_sid,
            }
        )

    # uc_pkl = Path("DADA2_UC.pkl.gz")
    # print(f"Exporting pre-UC table to {uc_pkl}")
    # with gzip.open(uc_pkl, "wb") as fh:
    #     pickle.dump(uc_pre, fh, protocol=pickle.HIGHEST_PROTOCOL)

    ## Summary stats
    num_asvs = n_asv
    num_asvreads = int(cluster_abunds.sum())
    num_merged = sum(
        1
        for r in uc_pre
        if r["SeqNumID"] is not None
        and r["ASV"] is not None
        and r["DerepSeqID"] != r["ASV"]
    )
    num_dsc = sum(1 for r in uc_pre if r["SeqNumID"] is None)
    num_dscreads = sum(r["Abundance"] for r in uc_pre if r["SeqNumID"] is None)
    perc_dsc = round(num_dsc / num_seqs * 100, 2) if num_seqs else 0.0
    perc_dscreads = round(num_dscreads / num_reads * 100, 2) if num_reads else 0.0

    print("\nRun summary:")
    print(f"Number of ASVs infered: {num_asvs}")
    print(f"Number of reads in ASV table: {num_asvreads}")
    print(f"Number of sequences merged into ASVs: {num_merged}")
    print(
        f"Number of discarded sequences (%): {num_dsc} ( {perc_dsc}% )"
    )
    print(
        f"Number of reads of discarded sequences (%): {num_dscreads} "
        f"( {perc_dscreads}% )"
    )

    ## FASTA (gzip, no line wrapping)
    print("\nExporting denoised sequences")
    fa_path = Path("DADA2_denoised.fa.gz")
    with gzip.open(fa_path, "wt", encoding="ascii", newline="\n") as out:
        for r in res_sorted:
            out.write(f">{r['SeqID']};size={r['Abundance']}\n")
            out.write(r["Sequence"] + "\n")

    ## Pseudo-UC file (tsv, gzip)
    print("\nExporting pseudo-UC file")
    uc_lines = []
    for r in uc_pre:
        if r["SeqNumID"] is None:
            continue
        rt = "S" if r["DerepSeqID"] == r["ASV"] else "H"
        # RecordType, ClustNum, SeqLen, Ident, Strand, V6, V7, ALN, DerepSeqID, ASV
        line = "\t".join(
            [
                rt,
                "",
                "",
                "",
                "+",
                "",
                "",
                ".",
                r["DerepSeqID"],
                r["ASV"] or "",
            ]
        )
        uc_lines.append(line)

    uc_gz = Path("DADA2_denoised.uc.gz")
    with gzip.open(uc_gz, "wt", encoding="utf-8", newline="\n") as out:
        out.write("\n".join(uc_lines))
        if uc_lines:
            out.write("\n")

    print("Exporting run statistics")
    smr_path = Path("DADA2_denoising_summary.txt")
    smr_rows = [
        ("Number of unique sequences (prior denoising)", num_seqs),
        ("Number of singleton sequences (prior denoising)", num_singl),
        ("Total abundance of sequences (prior denoising)", num_reads),
        ("Percentage of non-singleton sequences (prior denoising)", perc_ns),
        ("Number of ASVs infered", num_asvs),
        ("Number of reads in ASV table", num_asvreads),
        (
            "Number of sequences merged into ASVs (excluding representative seqs)",
            num_merged,
        ),
        ("Number of discarded sequences", num_dsc),
        ("Percentage of discarded sequences", perc_dsc),
        ("Number of reads of discarded sequences", num_dscreads),
        ("Percentage of reads of discarded sequences", perc_dscreads),
    ]
    with smr_path.open("w", encoding="utf-8") as fh:
        for name, val in smr_rows:
            fh.write(f"{name}\t{val}\n")

    elapsed = (time.time() - start) / 60.0
    print(f"\nElapsed time: {elapsed:.4f} minutes")
    print("\nAll done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
