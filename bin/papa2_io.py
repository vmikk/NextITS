"""
Build a papa2-compatible derep dict (`seqs`, `abundances`, `quals`, `map`) 
with unique sequences sorted by descending abundance
(same convention as `papa2.derep_fastq`)

Input:
- NextITS-style dereplicated FASTA/FASTQ: headers `SeqID;size=ABUNDANCE`
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from typing import BinaryIO, Iterator, List, Tuple

import numpy as np


@dataclass
class NextITSRecord:
    """One dereplicated input record after parsing the NextITS header."""

    seq_id: str
    abundance: int
    sequence: str
    qual_ascii: bytes | None = None


def _open_input(path: str) -> BinaryIO:
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def _parse_header_size(header_line: bytes) -> Tuple[str, int]:
    """Parse `@SeqID;size=N`, `>SeqID;size=N`, or `SeqID;size=N`."""
    h = header_line.strip()
    if h.startswith((b"@", b">")):
        h = h[1:]
    text = h.decode("ascii", errors="replace")
    if ";size=" not in text:
        raise ValueError(
            f"Expected ';size=' in dereplicated header, got: {text[:120]!r}"
        )
    seq_id, rest = text.split(";size=", 1)
    rest = rest.split()[0] if rest.split() else rest
    abundance = int(float(rest))
    return seq_id, abundance


def _detect_seq_format(path: str) -> str:
    with _open_input(path) as fh:
        while True:
            line = fh.readline()
            if not line:
                raise ValueError(f"Input file is empty: {path}")
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(b">"):
                return "fasta"
            if stripped.startswith(b"@"):
                return "fastq"
            raise ValueError(
                "Could not auto-detect dereplicated input format from the "
                f"first non-empty line of {path!r}: {stripped[:120]!r}"
            )


def _iter_fastq_records(path: str) -> Iterator[NextITSRecord]:
    with _open_input(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            if not header.strip():
                continue
            seq_line = fh.readline()
            plus = fh.readline()
            qual_line = fh.readline()
            if not qual_line:
                raise ValueError("Incomplete FASTQ record at end of file")
            if not plus.startswith(b"+"):
                raise ValueError(
                    f"Expected '+' line in FASTQ record, got: {plus[:120]!r}"
                )
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


def _iter_fasta_records(path: str) -> Iterator[NextITSRecord]:
    with _open_input(path) as fh:
        seq_id: str | None = None
        abundance: int | None = None
        seq_chunks: List[bytes] = []

        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(b">"):
                if seq_id is not None:
                    sequence = b"".join(seq_chunks).decode("ascii").upper()
                    yield NextITSRecord(
                        seq_id=seq_id,
                        abundance=int(abundance),
                        sequence=sequence,
                    )
                seq_id, abundance = _parse_header_size(stripped)
                seq_chunks = []
                continue
            if seq_id is None:
                raise ValueError(
                    f"Expected FASTA header line starting with '>', got: {stripped[:120]!r}"
                )
            seq_chunks.append(stripped)

        if seq_id is not None:
            sequence = b"".join(seq_chunks).decode("ascii").upper()
            yield NextITSRecord(
                seq_id=seq_id,
                abundance=int(abundance),
                sequence=sequence,
            )


def _qual_bytes_to_floats(q: bytes) -> np.ndarray:
    return np.frombuffer(q, dtype=np.uint8).astype(np.float64) - 33.0


def _empty_nextits_derep() -> Tuple[dict, dict]:
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


def _build_constant_quals(seqs: List[str], q_value: float = 40.0) -> np.ndarray:
    maxlen = max((len(s) for s in seqs), default=0)
    quals = np.full((len(seqs), maxlen), np.nan, dtype=np.float64)
    for i, seq in enumerate(seqs):
        quals[i, : len(seq)] = q_value
    return quals


def load_nextits_derep(path: str) -> Tuple[dict, dict]:
    """Load NextITS dereplicated FASTA/FASTQ into a papa2 derep dict plus metadata.

    Returns:
        derep: dict with keys ``seqs``, ``abundances``, ``quals``, ``map``.
        meta: dict with:
            - ``num_seqs``, ``num_singl``, ``num_reads``, ``perc_nonsingleton``
            - ``seq_ids_file_order``: list of SeqID per input record (file order)
            - ``abundances_file_order``: np.ndarray abundances per input record
            - ``sequences_file_order``: list of sequences per input record
            - ``unique_index_file_order``: for each input record, index into
              sorted uniques (after merge + abundance sort)
    """
    input_format = _detect_seq_format(path)
    record_iter = _iter_fasta_records if input_format == "fasta" else _iter_fastq_records
    records: List[NextITSRecord] = list(record_iter(path))
    if not records:
        return _empty_nextits_derep()

    seq_ids_file = [r.seq_id for r in records]
    abunds_file = np.array([r.abundance for r in records], dtype=np.int64)
    seqs_file = [r.sequence for r in records]

    # Stats in R sense: one row per input record (unique sequence line)
    num_seqs = len(records)
    num_singl = int(np.sum(abunds_file < 2))
    num_reads = int(abunds_file.sum())
    perc_nonsingleton = (
        round((num_seqs - num_singl) / num_seqs * 100, 2) if num_seqs else 0.0
    )

    # Merge identical sequences (sum abundances, weighted mean qual for FASTQ)
    seq_to_idx: dict = {}
    merged_abund: List[int] = []
    merged_qual_sum: List[np.ndarray] = []
    merged_qual_weight: List[int] = []

    for r in records:
        s = r.sequence
        if s not in seq_to_idx:
            idx = len(merged_abund)
            seq_to_idx[s] = idx
            merged_abund.append(0)
            if input_format == "fastq":
                slen = len(s)
                merged_qual_sum.append(np.zeros(slen, dtype=np.float64))
                merged_qual_weight.append(0)
        idx = seq_to_idx[s]
        merged_abund[idx] += r.abundance
        if input_format == "fastq":
            qf = _qual_bytes_to_floats(r.qual_ascii or b"")
            merged_qual_sum[idx][: len(qf)] += qf * r.abundance
            merged_qual_weight[idx] += r.abundance

    n_u = len(merged_abund)
    uniq_seqs = [None] * n_u  # type: ignore
    for s, idx in seq_to_idx.items():
        uniq_seqs[idx] = s

    abundances = np.zeros(n_u, dtype=np.int32)
    for i in range(n_u):
        abundances[i] = int(merged_abund[i])
    if input_format == "fastq":
        maxlen = max(len(s) for s in uniq_seqs) if uniq_seqs else 0
        quals = np.full((n_u, maxlen), np.nan, dtype=np.float64)
        for i in range(n_u):
            slen = len(uniq_seqs[i])
            if merged_qual_weight[i] > 0:
                quals[i, :slen] = (
                    merged_qual_sum[i][:slen] / merged_qual_weight[i]
                )
            else:
                quals[i, :slen] = np.nan
    else:
        quals = _build_constant_quals(uniq_seqs)

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
