"""Tests for the framework-free core: sequence handling, metrics, and PIP enrichment.

These run without TensorFlow, a GPU, or any downloaded data — they exercise every piece of logic
that isn't the model forward pass. Run:  python -m pytest tests/ -q  (or python tests/test_core.py)
"""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import sequtils            # noqa: E402


class FakeContig:
    def __init__(self, seq): self.seq = seq
    def __len__(self): return len(self.seq)
    def __getitem__(self, sl): return self.seq[sl]


class FakeFasta:
    def __init__(self, d): self.d = {k: FakeContig(v) for k, v in d.items()}
    def __getitem__(self, k): return self.d[k]


def test_one_hot_order_and_n():
    x = sequtils.one_hot("ACGTN")
    assert x.shape == (5, 4)
    assert np.array_equal(x[0], [1, 0, 0, 0])   # A
    assert np.array_equal(x[3], [0, 0, 0, 1])   # T
    assert np.array_equal(x[4], [0, 0, 0, 0])   # N -> zero
    assert x.sum() == 4


def test_revcomp():
    assert sequtils.revcomp("AACG") == "CGTT"


def test_fetch_window_centered_and_padded():
    fa = FakeFasta({"chr1": "AAAACCCCGGGGTTTT"})   # length 16
    w = sequtils.fetch_window(fa, "chr1", 8, 4)     # centered on index 8, half=2 -> [6:10]
    assert w == "CCGG" and len(w) == 4
    left = sequtils.fetch_window(fa, "chr1", 0, 6)  # runs off the left edge -> N-padded
    assert left.startswith("NNN") and len(left) == 6


def test_apply_snv_and_mismatch_guard():
    seq = "ACGTACGT"
    assert sequtils.apply_snv(seq, 2, "G", "A") == "ACATACGT"
    try:
        sequtils.apply_snv(seq, 2, "T", "A")        # wrong ref -> must raise
    except ValueError:
        pass
    else:
        raise AssertionError("apply_snv did not guard against ref mismatch")


def test_logfc_counts_sign_and_scale():
    # doubling counts => +1 in log2 space
    lc = sequtils.logfc_counts(np.log(10.0), np.log(20.0))
    assert abs(lc - 1.0) < 1e-9


def test_profile_jsd_bounds():
    a = np.array([10.0, 0.0, 0.0, 0.0])
    assert sequtils.profile_jsd(a, a) < 1e-12          # identical -> 0
    b = np.array([0.0, 0.0, 0.0, 10.0])
    j = sequtils.profile_jsd(a, b)                      # disjoint -> ~1 (base 2)
    assert 0.99 < j <= 1.0


def test_pip_enrichment_monotone():
    import pandas as pd
    from src import pip_enrichment
    # construct data where high PIP clearly rides with high |logfc|
    rng = np.random.default_rng(0)
    n = 2000
    pip = rng.uniform(0, 1, n)
    logfc = rng.normal(0, 0.1, n) + pip * 2.0          # score grows with PIP
    df = pd.DataFrame({"PIP": pip, "logfc_counts": logfc})
    cfg = {"variant_scoring": {"high_score_quantile": 0.90},
           "pip": {"bins": [0.0, 0.1, 0.5, 0.9, 1.0]}}
    res = pip_enrichment.enrichment(cfg, df=df)
    top = res.iloc[-1]["odds_ratio_vs_rest"]
    assert top > 1.0                                    # highest PIP bin enriched for high scores


def _run_all():
    fns = [v for k, v in globals().items() if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
        print(f"ok  {fn.__name__}")
    print(f"\n{len(fns)} tests passed")


if __name__ == "__main__":
    _run_all()
