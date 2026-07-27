"""The new bit: do DNase and ATAC HepG2 models agree on COUNTERFACTUAL variant effects?

Reference-sequence agreement after bias correction is reported in the paper (profile JSD
0.81 -> 0.26). Agreement on *variant* predictions at fine-mapped GWAS variants is a different
question and, as far as I can tell, unreported. This module joins the two per-assay score tables
and flags disagreers — candidates for either model instability or a real profile-shape-vs-coverage
distinction.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from .config import load_config, REPO_ROOT

KEY = ["chrom", "pos", "ref", "alt", "variant_id"]


def merge_assays(cfg) -> pd.DataFrame:
    proc = REPO_ROOT / cfg["paths"]["processed"]
    d = pd.read_parquet(proc / "variant_scores_dnase.parquet")
    a = pd.read_parquet(proc / "variant_scores_atac.parquet")
    cols = [c for c in KEY if c in d.columns]
    m = d.merge(a, on=cols, suffixes=("_dnase", "_atac"))
    # standardized disagreement on the magnitude head
    for col in ("logfc_counts", "jsd_profile"):
        dz = (m[f"{col}_dnase"] - m[f"{col}_dnase"].mean()) / m[f"{col}_dnase"].std()
        az = (m[f"{col}_atac"] - m[f"{col}_atac"].mean()) / m[f"{col}_atac"].std()
        m[f"disagree_{col}"] = np.abs(dz - az)
    return m


def summarize(cfg, m: pd.DataFrame) -> dict:
    return {
        "n": len(m),
        "pearson_logfc": pearsonr(m["logfc_counts_dnase"], m["logfc_counts_atac"])[0],
        "spearman_logfc": spearmanr(m["logfc_counts_dnase"], m["logfc_counts_atac"])[0],
        "pearson_jsd": pearsonr(m["jsd_profile_dnase"], m["jsd_profile_atac"])[0],
    }


def main():
    cfg = load_config()
    m = merge_assays(cfg)
    proc = REPO_ROOT / cfg["paths"]["processed"]
    m.to_parquet(proc / "assay_concordance.parquet")
    disagreers = m.sort_values("disagree_logfc_counts", ascending=False).head(25)
    disagreers.to_csv(REPO_ROOT / cfg["paths"]["tables"] / "assay_disagreers.csv", index=False)
    print(summarize(cfg, m))
    print(f"Top disagreers -> tables/assay_disagreers.csv")


if __name__ == "__main__":
    main()
