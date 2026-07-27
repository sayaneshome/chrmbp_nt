"""Figures. Standards (rcParams, colorblind-safe palette, labeled axes, tight_layout) follow the
repo CLAUDE.md and are applied here once so every figure inherits them.

Core figures:
  - pip_enrichment.png            : high-score enrichment across PIP bins (the headline)
  - atac_vs_dnase_variant_scores  : the assay-concordance question
  - cagi5_scatter                 : SORT1 / LDLR correctness check
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 12,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "legend.frameon": False,
})

OKABE_ITO = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442"]


def pip_enrichment_figure(enr: pd.DataFrame, out_path):
    """Bar of high-score odds ratio (bin vs rest) across PIP bins; * marks Fisher p<0.05."""
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.bar(enr["pip_bin"].astype(str), enr["odds_ratio_vs_rest"], color=OKABE_ITO[0])
    ax.axhline(1.0, ls="--", color="0.5", lw=1)
    for i, r in enr.reset_index(drop=True).iterrows():
        if pd.notna(r.get("fisher_p")) and r["fisher_p"] < 0.05:
            ax.text(i, r["odds_ratio_vs_rest"], "*", ha="center", va="bottom", fontsize=14)
    ax.set_xlabel("Fine-mapping posterior inclusion probability (PIP) bin")
    ax.set_ylabel("Odds ratio of high ChromBPNet score\n(bin vs. rest)")
    ax.set_title("HepG2 variant scores vs lipid fine-mapping")
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)


def assay_scatter(m: pd.DataFrame, out_path):
    """DNase vs ATAC counterfactual logfc; disagreers highlighted."""
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(m["logfc_counts_dnase"], m["logfc_counts_atac"], s=10, alpha=0.5,
               color=OKABE_ITO[0], label="fine-mapped variants")
    top = m.sort_values("disagree_logfc_counts", ascending=False).head(25)
    ax.scatter(top["logfc_counts_dnase"], top["logfc_counts_atac"], s=28,
               color=OKABE_ITO[3], label="top disagreers")
    lims = [m[["logfc_counts_dnase", "logfc_counts_atac"]].min().min(),
            m[["logfc_counts_dnase", "logfc_counts_atac"]].max().max()]
    ax.plot(lims, lims, "--", color="0.5", lw=1, zorder=0)
    ax.set_xlabel("DNase model: log2 fold-change in counts (alt/ref)")
    ax.set_ylabel("ATAC model: log2 fold-change in counts (alt/ref)")
    ax.set_title("Do the two assays agree on variant effects?")
    ax.legend(fontsize=9, loc="upper left")
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)
