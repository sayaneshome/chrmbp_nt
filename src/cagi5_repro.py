"""Correctness check: reproduce the paper's HepG2 CAGI5 MPRA correlations at SORT1 / LDLR / F9.

Score every CAGI5 saturation-MPRA SNV at the control loci with the HepG2 DNase model, correlate
predicted effect (logfc_counts) against measured MPRA activity, per locus. Run FIRST — if the
control loci don't reproduce the Fig. 7d / supplement values, nothing downstream is trustworthy.
rs12740374 (the C/EBP-site-creating SORT1 variant) is the canonical anchor.

Expects the satMutMPRA export lifted to hg38 (or exported as hg38 directly). CONFIRM the portal's
column names and set COLS below; the portal's effect column is typically 'Value' (log2 expression).
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from .config import load_config, REPO_ROOT
from . import sequtils, variant_scoring

# satMutMPRA export columns. CONFIRM against the downloaded TSV; override as needed.
COLS = dict(chrom="Chromosome", pos="Position", ref="Ref", alt="Alt",
            locus="Element", effect="Value")


def _score_variants(cfg, sub: pd.DataFrame, model, fasta, L: int) -> np.ndarray:
    """Return predicted logfc_counts for each SNV row in `sub` (NaN where it can't be scored)."""
    half = L // 2
    ref_seqs, alt_seqs, idx = [], [], []
    for i, v in sub.iterrows():
        try:
            ref_win = sequtils.fetch_window(fasta, str(v[COLS["chrom"]]), int(v[COLS["pos"]]) - 1, L)
            alt_win = sequtils.apply_snv(ref_win, half, str(v[COLS["ref"]]), str(v[COLS["alt"]]))
        except (KeyError, ValueError):
            continue
        ref_seqs.append(ref_win); alt_seqs.append(alt_win); idx.append(i)
    pred = pd.Series(np.nan, index=sub.index)
    if not idx:
        return pred.to_numpy()
    _, ref_lc = model.predict(sequtils.one_hot_batch(ref_seqs))
    _, alt_lc = model.predict(sequtils.one_hot_batch(alt_seqs))
    pred.loc[idx] = [sequtils.logfc_counts(r, a) for r, a in zip(ref_lc, alt_lc)]
    return pred.to_numpy()


def run(cfg) -> pd.DataFrame:
    from pyfaidx import Fasta
    fasta = Fasta(str(REPO_ROOT / cfg["genome"]["fasta"]))
    model = variant_scoring.load_model(cfg, "dnase")   # CAGI5 controls are matched to DNase
    L = cfg["variant_scoring"]["input_len"]

    mpra = pd.read_csv(REPO_ROOT / cfg["cagi5"]["mpra"], sep="\t")
    mpra = mpra[(mpra[COLS["ref"]].str.len() == 1) & (mpra[COLS["alt"]].str.len() == 1)]

    rows = []
    for locus in cfg["cagi5"]["control_loci"]:
        sub = mpra[mpra[COLS["locus"]] == locus].copy()
        if sub.empty:
            rows.append({"locus": locus, "n": 0, "pearson": np.nan, "spearman": np.nan})
            continue
        sub["pred"] = _score_variants(cfg, sub, model, fasta, L)
        d = sub.dropna(subset=["pred", COLS["effect"]])
        rows.append({
            "locus": locus, "n": len(d),
            "pearson": pearsonr(d["pred"], d[COLS["effect"]])[0] if len(d) > 2 else np.nan,
            "spearman": spearmanr(d["pred"], d[COLS["effect"]])[0] if len(d) > 2 else np.nan,
        })
    return pd.DataFrame(rows)


def main():
    cfg = load_config()
    res = run(cfg)
    out = REPO_ROOT / cfg["paths"]["tables"] / "cagi5_reproduction.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out, index=False)
    print(res.to_string(index=False))
    print(f"\nCompare against paper Fig. 7d / supplement (Zenodo 14735906) -> {out}")


if __name__ == "__main__":
    main()
