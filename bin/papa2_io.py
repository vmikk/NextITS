"""
Build a papa2-compatible derep dict (`seqs`, `abundances`, `quals`, `map`) 
with unique sequences sorted by descending abundance
(same convention as `papa2.derep_fastq`)

Input:
- NextITS-style dereplicated FASTA/FASTQ: headers `SeqID;size=ABUNDANCE`

Also provides a quality-free error function with a pseudocount
(port of R dada2 `noqualErrfun`) and error-model diagnostics
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
    group: str = ""


def _open_input(path: str) -> BinaryIO:
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def _parse_header_size(header_line: bytes) -> Tuple[str, int, str]:
    """Parse `@SeqID;size=N`, `>SeqID;size=N`, or `SeqID;size=N`.

    An optional group label (`;grp=LABEL`, added to the error-learning subset)
    is returned as the third element (empty string if absent).
    """
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
    group = ""
    if ";grp=" in rest:
        group = rest.split(";grp=", 1)[1].split(";")[0]
    abundance = int(float(rest.split(";")[0]))
    return seq_id, abundance, group


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
            seq_id, abundance, group = _parse_header_size(header)
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
                group=group,
            )


def _iter_fasta_records(path: str) -> Iterator[NextITSRecord]:
    with _open_input(path) as fh:
        seq_id: str | None = None
        abundance: int | None = None
        group = ""
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
                        group=group,
                    )
                seq_id, abundance, group = _parse_header_size(stripped)
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
                group=group,
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


def load_nextits_derep_groups(path: str) -> List[Tuple[str, dict, dict]]:
    """Load NextITS dereplicated FASTA/FASTQ, split by group label (`;grp=`).

    Returns a list of (group, derep, meta) tuples, one per group,
    or a single tuple with an empty group label if the input has no labels.
    Used for error learning, where groups are treated as separate samples.
    """
    input_format = _detect_seq_format(path)
    record_iter = _iter_fasta_records if input_format == "fasta" else _iter_fastq_records
    by_group: dict = {}
    for r in record_iter(path):
        by_group.setdefault(r.group, []).append(r)
    return [
        (g, *load_nextits_derep(path, records=recs, input_format=input_format))
        for g, recs in sorted(by_group.items())
    ]


def load_nextits_derep(
    path: str, records: List[NextITSRecord] | None = None, input_format: str | None = None
) -> Tuple[dict, dict]:
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
    if records is None:
        input_format = _detect_seq_format(path)
        record_iter = _iter_fasta_records if input_format == "fasta" else _iter_fastq_records
        records = list(record_iter(path))
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


_NT = "ACGT"


def noqual_errfun_pc(trans, pseudocount: float = 1.0) -> np.ndarray:
    """Estimate error rates ignoring quality scores (constant across Q columns).

    Port of R dada2 ``noqualErrfun`` (``R/errorModels.R``), including its pseudocount.
    ``papa2.noqual_errfun`` has no pseudocount, so a transition that was never observed gets a rate of exactly 0
    then lambda = 0 and every sequence with that substitution is split off as a new ASV.
    """
    trans = np.asarray(trans, dtype=np.float64)
    ncol = trans.shape[1]
    obs = trans.sum(axis=1) + pseudocount
    err = np.zeros((16, ncol), dtype=np.float64)
    for nti in range(4):
        rows = range(nti * 4, nti * 4 + 4)
        tot = sum(obs[r] for r in rows)
        for r in rows:
            if r != nti * 5:
                err[r, :] = obs[r] / tot
        err[nti * 5, :] = 1.0 - sum(err[r, 0] for r in rows if r != nti * 5)
    return err


def count_substitutions(trans) -> int:
    """Total number of observed substitutions (off-diagonal transitions)."""
    trans = np.asarray(trans, dtype=np.float64)
    diag = [0, 5, 10, 15]
    return int(round(trans.sum() - trans[diag, :].sum()))


def format_error_rates(err) -> List[str]:
    """One line per substitution type with its rate (first Q column)."""
    err = np.asarray(err, dtype=np.float64)
    lines = []
    for nti in range(4):
        for ntj in range(4):
            if nti != ntj:
                lines.append(f"  {_NT[nti]}2{_NT[ntj]}: {err[nti * 4 + ntj, 0]:.4e}")
    return lines
