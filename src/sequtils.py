"""Sequence + metric utilities. Pure numpy/scipy — no DL framework, so these run and are tested
(see tests/test_core.py). Everything the model-scoring code needs that isn't the model itself.

Conventions:
  - one-hot channel order is A, C, G, T (ChromBPNet convention); N and other chars -> all-zero.
  - ChromBPNet input length is 2114 bp; profile head length is 1000 bp.
"""
from __future__ import annotations
import numpy as np

BASES = "ACGT"
_LOOKUP = {b: i for i, b in enumerate(BASES)}


def one_hot(seq: str) -> np.ndarray:
    """(L, 4) float32 one-hot. Case-insensitive; non-ACGT -> zero row."""
    x = np.zeros((len(seq), 4), dtype=np.float32)
    for i, b in enumerate(seq.upper()):
        j = _LOOKUP.get(b)
        if j is not None:
            x[i, j] = 1.0
    return x


def one_hot_batch(seqs) -> np.ndarray:
    """(N, L, 4). All seqs must share length L."""
    seqs = list(seqs)
    L = len(seqs[0])
    if any(len(s) != L for s in seqs):
        raise ValueError("all sequences must have equal length")
    out = np.zeros((len(seqs), L, 4), dtype=np.float32)
    for k, s in enumerate(seqs):
        out[k] = one_hot(s)
    return out


def revcomp(seq: str) -> str:
    return seq.upper().translate(str.maketrans("ACGT", "TGCA"))[::-1]


def fetch_window(fasta, chrom: str, center0: int, length: int) -> str:
    """Fetch a `length`-bp window centered on 0-based coordinate `center0`.

    `fasta` is a pyfaidx.Fasta (or anything indexable as fasta[chrom][start:end] -> seq-like).
    Out-of-bounds windows are edge-padded with N so the returned string is always `length` bp.
    """
    half = length // 2
    start = center0 - half
    end = start + length
    chrom_len = len(fasta[chrom])
    lpad = max(0, -start)
    rpad = max(0, end - chrom_len)
    seq = str(fasta[chrom][max(0, start):min(chrom_len, end)])
    return "N" * lpad + seq + "N" * rpad


def apply_snv(seq: str, offset: int, ref: str, alt: str) -> str:
    """Return `seq` with the single base at `offset` swapped ref->alt.

    Validates that seq[offset] matches `ref` (case-insensitive); raises on mismatch so a
    wrong-strand or wrong-coordinate variant fails loudly instead of silently scoring garbage.
    Only single-nucleotide variants are handled here (indels need length-aware windows).
    """
    if len(ref) != 1 or len(alt) != 1:
        raise ValueError(f"apply_snv handles SNVs only, got ref={ref!r} alt={alt!r}")
    if seq[offset].upper() != ref.upper():
        raise ValueError(f"ref mismatch at offset {offset}: seq has {seq[offset]!r}, expected {ref!r}")
    return seq[:offset] + alt.upper() + seq[offset + 1:]


# --- metrics on model outputs ------------------------------------------------

def softmax(logits: np.ndarray, axis: int = -1) -> np.ndarray:
    z = logits - np.max(logits, axis=axis, keepdims=True)
    e = np.exp(z)
    return e / np.sum(e, axis=axis, keepdims=True)


def profile_jsd(ref_logits: np.ndarray, alt_logits: np.ndarray) -> float:
    """Jensen-Shannon divergence (base 2, in [0,1]) between two profile-shape distributions.

    Inputs are raw profile logits (length = profile head, e.g. 1000). Softmax'd to probabilities
    first. scipy.jensenshannon returns the JS *distance* (sqrt of divergence); we square it to
    report divergence, matching the JSD values quoted in the ChromBPNet paper.
    """
    from scipy.spatial.distance import jensenshannon
    p = softmax(np.asarray(ref_logits, dtype=np.float64))
    q = softmax(np.asarray(alt_logits, dtype=np.float64))
    dist = jensenshannon(p, q, base=2)
    return float(dist ** 2)


def logfc_counts(ref_logcount: float, alt_logcount: float) -> float:
    """log2 fold-change in predicted counts.

    ChromBPNet's counts head emits natural-log counts, so the log2 fold-change is the
    difference in nat-log space divided by ln(2). Sign: positive => alt more accessible.
    """
    return float((alt_logcount - ref_logcount) / np.log(2.0))
