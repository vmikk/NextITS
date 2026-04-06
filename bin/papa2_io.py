"""
Build a papa2-compatible derep dict (`seqs`, `abundances`, `quals`, `map`) 
with unique sequences sorted by descending abundance
(same convention as `papa2.derep_fastq`)

Input:
- NextITS-style dereplicated FASTQ: headers `SeqID;size=ABUNDANCE`
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from typing import BinaryIO, Iterator, List, Tuple

import numpy as np


@dataclass
class NextITSRecord:
    """One FASTQ record after parsing the NextITS header."""

    seq_id: str
    abundance: int
    sequence: str
    qual_ascii: bytes  # raw quality string (same length as sequence)


def _open_fastq(path: str) -> BinaryIO:
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def _parse_header_size(header_line: bytes) -> Tuple[str, int]:
    """Parse `@SeqID;size=N` or `SeqID;size=N` (strip @)."""
    h = header_line.strip()
    if h.startswith(b"@"):
        h = h[1:]
    text = h.decode("ascii", errors="replace")
    if ";size=" not in text:
        raise ValueError(
            f"Expected ';size=' in FASTQ header, got: {text[:120]!r}"
        )
    seq_id, rest = text.split(";size=", 1)
    rest = rest.split()[0] if rest.split() else rest
    abundance = int(float(rest))
    return seq_id, abundance


def _iter_fastq_records(path: str) -> Iterator[NextITSRecord]:
    with _open_fastq(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq_line = fh.readline()
            plus = fh.readline()
            qual_line = fh.readline()
            if not qual_line:
                break
            seq_id, abundance = _parse_header_size(header)
            seq = seq_line.strip().decode("ascii").upper()
            q = qual_line.rstrip(b"\n\r")
            if len(q) != len(seq):
                raise ValueError(
                    f"Qual length {len(q)} != seq length {len(seq)} for {seq_id!r}"
                )
            yield NextITSRecord(
                seq_id=seq_id,
                abundance=abundance,
                sequence=seq,
                qual_ascii=q,
            )


def _qual_bytes_to_floats(q: bytes) -> np.ndarray:
    return np.frombuffer(q, dtype=np.uint8).astype(np.float64) - 33.0


def load_nextits_derep_fastq(path: str) -> Tuple[dict, dict]:
    """Load NextITS dereplicated FASTQ into a papa2 derep dict plus metadata.

    Returns:
        derep: dict with keys ``seqs``, ``abundances``, ``quals``, ``map``.
        meta: dict with:
            - ``num_seqs``, ``num_singl``, ``num_reads``, ``perc_nonsingleton``
            - ``seq_ids_file_order``: list of SeqID per FASTQ record (file order)
            - ``abundances_file_order``: np.ndarray abundances per record
            - ``sequences_file_order``: list of sequences per record
            - ``unique_index_file_order``: for each file record, index into
              sorted uniques (after merge + abundance sort)
    """
    records: List[NextITSRecord] = list(_iter_fastq_records(path))
    if not records:
        empty = {
            "seqs": [],
            "abundances": np.array([], dtype=np.int32),
            "quals": np.zeros((0, 0), dtype=np.float64),
            "map": np.array([], dtype=np.int32),
        }
        meta = {
            "num_seqs": 0,
            "num_singl": 0,
            "num_reads": 0,
            "perc_nonsingleton": 0.0,
            "seq_ids_file_order": [],
            "abundances_file_order": np.array([], dtype=np.int64),
            "sequences_file_order": [],
            "unique_index_file_order": np.array([], dtype=np.int32),
        }
        return empty, meta

    seq_ids_file = [r.seq_id for r in records]
    abunds_file = np.array([r.abundance for r in records], dtype=np.int64)
    seqs_file = [r.sequence for r in records]

    # Stats in R sense: one row per FASTQ record (unique sequence line)
    num_seqs = len(records)
    num_singl = int(np.sum(abunds_file < 2))
    num_reads = int(abunds_file.sum())
    perc_nonsingleton = (
        round((num_seqs - num_singl) / num_seqs * 100, 2) if num_seqs else 0.0
    )

    # Merge identical sequences (sum abundances, weighted mean qual)
    seq_to_idx: dict = {}
    merged_abund: List[int] = []
    merged_qual_sum: List[np.ndarray] = []
    merged_qual_weight: List[int] = []

    for r in records:
        s = r.sequence
        qf = _qual_bytes_to_floats(r.qual_ascii)
        if s not in seq_to_idx:
            idx = len(merged_abund)
            seq_to_idx[s] = idx
            merged_abund.append(0)
            slen = len(s)
            merged_qual_sum.append(np.zeros(slen, dtype=np.float64))
            merged_qual_weight.append(0)
        idx = seq_to_idx[s]
        merged_abund[idx] += r.abundance
        merged_qual_sum[idx][: len(qf)] += qf * r.abundance
        merged_qual_weight[idx] += r.abundance

    n_u = len(merged_abund)
    uniq_seqs = [None] * n_u  # type: ignore
    for s, idx in seq_to_idx.items():
        uniq_seqs[idx] = s

    maxlen = max(len(s) for s in uniq_seqs) if uniq_seqs else 0
    quals = np.full((n_u, maxlen), np.nan, dtype=np.float64)
    abundances = np.zeros(n_u, dtype=np.int32)
    for i in range(n_u):
        abundances[i] = int(merged_abund[i])
        slen = len(uniq_seqs[i])
        if merged_qual_weight[i] > 0:
            quals[i, :slen] = merged_qual_sum[i][:slen] / merged_qual_weight[i]
        else:
            quals[i, :slen] = np.nan

    # Sort by abundance descending (papa2 / derep_fastq convention)
    order = np.argsort(-abundances, kind="mergesort")
    seqs = [uniq_seqs[int(j)] for j in order]
    abundances = abundances[order]
    quals = quals[order, :]

    sorted_idx_of_seq = {seqs[j]: j for j in range(len(seqs))}
    unique_index_file_order = np.array(
        [sorted_idx_of_seq[s] for s in seqs_file], dtype=np.int32
    )

    derep = {
        "seqs": seqs,
        "abundances": abundances,
        "quals": quals,
        # Unused by C core for inference; present for API parity with derep_fastq
        "map": np.arange(len(seqs), dtype=np.int32),
    }

    meta = {
        "num_seqs": num_seqs,
        "num_singl": num_singl,
        "num_reads": num_reads,
        "perc_nonsingleton": perc_nonsingleton,
        "seq_ids_file_order": seq_ids_file,
        "abundances_file_order": abunds_file,
        "sequences_file_order": seqs_file,
        "unique_index_file_order": unique_index_file_order,
    }
    return derep, meta


def subsample_derep_by_nbases(
    derep: dict, target_bases: float, max_reads_per_seq: int = 0
) -> dict:
    """Take a subset of uniques (highest abundance first) until >= target_bases.

    Approximates DADA2 ``learnErrors(..., nbases=...)`` read/budget behavior
    for a single pre-dereplicated sample.
    Optionally caps the effective abundance of each unique during learning 
    to avoid the budget being consumed by a handful of extremely abundant sequences.
    """
    seqs = derep["seqs"]
    abunds = derep["abundances"]
    quals = derep["quals"]
    if not seqs:
        return dict(derep)
    if max_reads_per_seq < 0:
        raise ValueError("max_reads_per_seq must be >= 0")

    learn_abunds = np.asarray(abunds, dtype=np.int32)
    if max_reads_per_seq > 0:
        learn_abunds = np.minimum(learn_abunds, int(max_reads_per_seq)).astype(
            np.int32, copy=False
        )

    total_bases = float(
        sum(int(learn_abunds[i]) * len(seqs[i]) for i in range(len(seqs)))
    )
    if total_bases <= target_bases:
        max_col = max(len(s) for s in seqs)
        new_quals = np.asarray(quals, dtype=np.float64)
        if new_quals.shape[1] > max_col:
            new_quals = new_quals[:, :max_col]
        elif new_quals.shape[1] < max_col:
            pad = max_col - new_quals.shape[1]
            new_quals = np.hstack(
                [new_quals, np.full((new_quals.shape[0], pad), np.nan)]
            )
        return {
            "seqs": list(seqs),
            "abundances": np.asarray(learn_abunds, dtype=np.int32),
            "quals": new_quals,
            "map": np.arange(len(seqs), dtype=np.int32),
        }

    cum = 0
    keep_idx: List[int] = []
    for i, seq in enumerate(seqs):
        cum += int(learn_abunds[i]) * len(seq)
        keep_idx.append(i)
        if cum >= target_bases:
            break
    new_seqs = [seqs[i] for i in keep_idx]
    new_abunds = np.array([int(learn_abunds[i]) for i in keep_idx], dtype=np.int32)
    new_quals = np.asarray(quals[keep_idx, :], dtype=np.float64)
    if new_seqs:
        max_col = max(len(s) for s in new_seqs)
        if new_quals.shape[1] > max_col:
            new_quals = new_quals[:, :max_col]
        elif new_quals.shape[1] < max_col:
            pad = max_col - new_quals.shape[1]
            new_quals = np.hstack(
                [new_quals, np.full((new_quals.shape[0], pad), np.nan)]
            )
    return {
        "seqs": new_seqs,
        "abundances": new_abunds,
        "quals": new_quals,
        "map": np.arange(len(new_seqs), dtype=np.int32),
    }
