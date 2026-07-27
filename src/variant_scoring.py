"""Counterfactual variant scoring with a HepG2 ChromBPNet model (one assay).

For each SNV: build ref and alt 2114 bp windows centered on the variant, run the bias-corrected
model, and record
    logfc_counts - log2 fold-change of predicted counts (magnitude of accessibility change)
    jsd_profile  - Jensen-Shannon divergence between ref and alt predicted profiles (shape change)

Only the model call needs the TF env; window building + metrics are the tested utilities in
sequtils.py. Both assays run this independently; assay_concordance.py compares them.

Output schema (stable across assays, matches the Kundaje variant-scorer if you swap that in):
    [variant_id, chrom, pos, ref, alt, trait, PIP, logfc_counts, jsd_profile, assay]
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from tqdm import tqdm

from .config import load_config, REPO_ROOT
from .chrombpnet_model import ChromBPNetModel
from . import sequtils

# Column names in the UKBB fine-mapping tables. CONFIRM against the release manifest; override
# here if the release uses different headers (e.g. 'variant', 'susie_pip').
COLS = dict(chrom="chromosome", pos="position", ref="allele1", alt="allele2",
            pip="pip", vid="variant")


def load_model(cfg, assay) -> ChromBPNetModel:
    return ChromBPNetModel(REPO_ROOT / cfg["models"][assay]["dir"])


def load_variants(cfg) -> pd.DataFrame:
    """Union of fine-mapped variants across lipid traits, normalized to the output schema.

    Assumes the per-trait tables have already been lifted hg19 -> hg38 (see liftover.py) and
    written next to the originals as *.hg38.tsv, OR that COLS point at hg38 coordinates.
    """
    frames = []
    for t in cfg["traits"]:
        path = REPO_ROOT / t["finemap"]
        df = pd.read_csv(path, sep="\t")
        std = pd.DataFrame({
            "variant_id": df[COLS["vid"]],
            "chrom": df[COLS["chrom"]].astype(str).str.replace("^chr", "", regex=True).radd("chr"),
            "pos": df[COLS["pos"]].astype(int),
            "ref": df[COLS["ref"]].astype(str),
            "alt": df[COLS["alt"]].astype(str),
            "PIP": df[COLS["pip"]].astype(float),
            "trait": t["name"],
        })
        frames.append(std)
    out = pd.concat(frames, ignore_index=True)
    # keep SNVs only; indels need length-aware windows and are out of scope here
    snv = (out["ref"].str.len() == 1) & (out["alt"].str.len() == 1)
    n_indel = int((~snv).sum())
    if n_indel:
        print(f"dropping {n_indel} non-SNV variants (indels out of scope)")
    return out[snv].reset_index(drop=True)


def score(cfg, assay: str) -> pd.DataFrame:
    from pyfaidx import Fasta
    fasta = Fasta(str(REPO_ROOT / cfg["genome"]["fasta"]))
    model = load_model(cfg, assay)
    variants = load_variants(cfg)
    L = cfg["variant_scoring"]["input_len"]
    half = L // 2

    ref_seqs, alt_seqs, kept_idx, skipped = [], [], [], 0
    for i, v in variants.iterrows():
        try:
            ref_win = sequtils.fetch_window(fasta, v["chrom"], v["pos"] - 1, L)  # pos is 1-based
            alt_win = sequtils.apply_snv(ref_win, half, v["ref"], v["alt"])
        except (KeyError, ValueError):
            skipped += 1                      # ref mismatch / missing contig -> skip loudly-counted
            continue
        ref_seqs.append(ref_win)
        alt_seqs.append(alt_win)
        kept_idx.append(i)
    if skipped:
        print(f"skipped {skipped} variants (ref-allele mismatch or missing contig)")

    ref_prof, ref_lc = model.predict(sequtils.one_hot_batch(ref_seqs))
    alt_prof, alt_lc = model.predict(sequtils.one_hot_batch(alt_seqs))

    rows = variants.loc[kept_idx].reset_index(drop=True).copy()
    rows["logfc_counts"] = [sequtils.logfc_counts(r, a) for r, a in zip(ref_lc, alt_lc)]
    rows["jsd_profile"] = [sequtils.profile_jsd(ref_prof[k], alt_prof[k])
                           for k in range(len(kept_idx))]
    rows["assay"] = assay
    return rows


def main(argv):
    cfg = load_config()
    assay = argv[1] if len(argv) > 1 else "dnase"
    df = score(cfg, assay)
    out = REPO_ROOT / cfg["paths"]["processed"] / f"variant_scores_{assay}.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out)
    print(f"{assay}: scored {len(df)} variants -> {out}")


if __name__ == "__main__":
    import sys
    main(sys.argv)
