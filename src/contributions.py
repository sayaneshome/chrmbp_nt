"""Mechanistic coherence: DeepSHAP contribution scores at top hits + HepG2 motif check.

For the top-scoring fine-mapped variants, compute base-resolution contribution scores around each
with DeepSHAP on the HepG2 model's COUNTS head, then check the sequence the model relies on matches
the HepG2 TF lexicon (HNF4A, FOXA, CEBP, HNF1A - paper Fig. 5a). Right motifs => coherence, not just
correlation.

ChromBPNet models are Keras, so attribution uses shap.DeepExplainer with dinucleotide-shuffled
references (the standard ChromBPNet contribution recipe), not tangermeme (which is PyTorch). Needs
the TF env; the top-hit selection and I/O here are framework-free.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd

from .config import load_config, REPO_ROOT
from . import sequtils, variant_scoring


def top_hits(cfg, assay="dnase") -> pd.DataFrame:
    df = pd.read_parquet(REPO_ROOT / cfg["paths"]["processed"] / f"variant_scores_{assay}.parquet")
    order = df["logfc_counts"].abs().sort_values(ascending=False).index
    return df.loc[order].head(cfg["contributions"]["n_top_loci"]).reset_index(drop=True)


def _dinuc_shuffled_refs(onehot: np.ndarray, n: int, seed: int) -> np.ndarray:
    """n dinucleotide-preserving shuffles of a single (L,4) sequence -> (n, L, 4) reference set."""
    from tangermeme.ersatz import dinucleotide_shuffle  # ersatz shuffle util; framework-agnostic array in/out
    import numpy as _np
    seq = onehot[None]  # (1, L, 4)
    shuf = dinucleotide_shuffle(_np.asarray(seq), n_shuffles=n, random_state=seed)
    return _np.asarray(shuf).reshape(n, onehot.shape[0], onehot.shape[1])


def deepshap(cfg, hits: pd.DataFrame, assay="dnase"):
    """Per-hit contribution scores (counts head) around each variant; saved as .npz."""
    import shap
    from pyfaidx import Fasta
    fasta = Fasta(str(REPO_ROOT / cfg["genome"]["fasta"]))
    model = variant_scoring.load_model(cfg, assay)._load()
    L = cfg["variant_scoring"]["input_len"]
    seed = cfg["random_seed"]
    outdir = REPO_ROOT / cfg["paths"]["processed"] / "contribs"
    outdir.mkdir(parents=True, exist_ok=True)

    # counts head is output index 1 (see chrombpnet_model.predict); explain that scalar
    counts_model = model  # DeepExplainer on the counts output
    for _, v in hits.iterrows():
        seq = sequtils.fetch_window(fasta, v["chrom"], int(v["pos"]) - 1, L)
        oh = sequtils.one_hot(seq)
        refs = _dinuc_shuffled_refs(oh, n=20, seed=seed)
        explainer = shap.DeepExplainer((counts_model.input, counts_model.outputs[1]), refs)
        shap_vals = explainer.shap_values(oh[None])[0][0]      # (L, 4)
        contrib = shap_vals * oh                                # projected onto the actual base
        np.savez(outdir / f"{v['variant_id']}.npz", contrib=contrib, onehot=oh)
    print(f"Saved contribution tracks for {len(hits)} hits -> {outdir}")


def motif_coherence(cfg):
    """Run modisco-lite over the saved contribution windows, then TOMTOM-match discovered motifs
    to the reference DB; report how many map to the expected HepG2 lexicon."""
    import modiscolite
    proc = REPO_ROOT / cfg["paths"]["processed"]
    files = sorted((proc / "contribs").glob("*.npz"))
    contribs = np.stack([np.load(f)["contrib"] for f in files])
    onehots = np.stack([np.load(f)["onehot"] for f in files])

    pos, neg = modiscolite.tfmodisco.TFMoDISco(
        hypothetical_contribs=contribs, one_hot=onehots,
        max_seqlets_per_metacluster=2000)
    h5 = proc / "modisco_results.h5"
    modiscolite.io.save_hdf5(h5, pos, neg)
    modiscolite.report.report_motifs(
        str(h5), str(proc / "modisco_report"),
        meme_motif_db=str(REPO_ROOT / cfg["contributions"]["motif_db"]))
    print(f"modisco report -> {proc/'modisco_report'}; "
          f"check matches to {cfg['contributions']['expected_motifs']}")


def main():
    cfg = load_config()
    hits = top_hits(cfg)
    out = REPO_ROOT / cfg["paths"]["processed"] / "top_hits.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    hits.to_parquet(out)
    print(f"Top {len(hits)} hits -> {out}")
    deepshap(cfg, hits)
    motif_coherence(cfg)


if __name__ == "__main__":
    main()
