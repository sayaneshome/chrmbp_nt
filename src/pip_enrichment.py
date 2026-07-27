"""Concordance: do high ChromBPNet variant scores enrich for high-PIP fine-mapped variants?

For each PIP bin, compute the fraction of variants that are "high-scoring" (|logfc_counts| above
the configured quantile) and the enrichment of that bin vs ALL OTHER variants (background) as a
2x2 odds ratio with a Fisher exact p-value. Bin-vs-rest is used rather than bin-vs-lowest so the
odds ratio stays defined even when a sparse bin has zero high-scorers. If the model tracks real
regulatory effects, the high-PIP bins are enriched for high scores.

Uses DNase scores by default (the CAGI5-validated assay). Pure pandas/scipy — runs and is tested.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from .config import load_config, REPO_ROOT


def enrichment(cfg, assay="dnase", df: pd.DataFrame = None) -> pd.DataFrame:
    if df is None:
        proc = REPO_ROOT / cfg["paths"]["processed"]
        df = pd.read_parquet(proc / f"variant_scores_{assay}.parquet")
    df = df.copy()

    thr = df["logfc_counts"].abs().quantile(cfg["variant_scoring"]["high_score_quantile"])
    df["high_score"] = df["logfc_counts"].abs() >= thr

    edges = cfg["pip"]["bins"]
    df["pip_bin"] = pd.cut(df["PIP"], bins=edges, include_lowest=True)
    bins = list(df["pip_bin"].cat.categories)

    rows = []
    for b in bins:
        in_bin = df["pip_bin"] == b
        g, rest = df[in_bin], df[~in_bin]
        hi, lo = int(g["high_score"].sum()), int((~g["high_score"]).sum())
        rhi, rlo = int(rest["high_score"].sum()), int((~rest["high_score"]).sum())
        # Haldane-Anscombe 0.5 correction keeps the OR finite when a cell is empty.
        odds = ((hi + 0.5) * (rlo + 0.5)) / ((lo + 0.5) * (rhi + 0.5))
        p = fisher_exact([[hi, lo], [rhi, rlo]], alternative="greater")[1] if len(g) else np.nan
        rows.append({"pip_bin": str(b), "n": len(g),
                     "high_score_rate": g["high_score"].mean() if len(g) else np.nan,
                     "odds_ratio_vs_rest": odds, "fisher_p": p})
    return pd.DataFrame(rows)


def main():
    cfg = load_config()
    res = enrichment(cfg)
    out = REPO_ROOT / cfg["paths"]["tables"] / "pip_enrichment.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out, index=False)
    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
