"""hg19 -> hg38 liftOver for the three hg19 inputs (CAGI5, UKBB fine-mapping, ABC regions).

Uses pyliftover (pure-python, reads the UCSC .chain.gz directly — no external binary). Variants
that fail to lift, or that map to a different chromosome, are dropped and COUNTED; the drop count
is returned so it can be logged (see the Limitations note in the README).
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd

from .config import load_config, REPO_ROOT

_LO_CACHE = {}


def _lifter(cfg):
    """Cached pyliftover.LiftOver built from the chain in config.yaml."""
    key = cfg["genome"]["liftover_chain"]
    if key not in _LO_CACHE:
        from pyliftover import LiftOver
        _LO_CACHE[key] = LiftOver(str(REPO_ROOT / key))
    return _LO_CACHE[key]


def lift_positions(cfg, df: pd.DataFrame, chrom_col="chrom", pos_col="pos",
                   pos_is_one_based=True) -> tuple[pd.DataFrame, int]:
    """Lift a variant/position table hg19 -> hg38.

    Adds hg38 `chrom`/`pos` in place of the originals (originals kept as *_hg19). Rows that fail
    to lift, or that switch chromosome, are dropped. Returns (lifted_df, n_dropped).

    pyliftover is 0-based; narrowPeak/BED positions are already 0-based, VCF/fine-mapping are
    typically 1-based — set pos_is_one_based accordingly.
    """
    lo = _lifter(cfg)
    new_chrom, new_pos, keep = [], [], []
    for chrom, pos in zip(df[chrom_col].astype(str), df[pos_col].astype(int)):
        c = chrom if chrom.startswith("chr") else f"chr{chrom}"
        p0 = pos - 1 if pos_is_one_based else pos
        res = lo.convert_coordinate(c, p0)
        if res and res[0][0] == c:              # lifted, same chromosome
            new_chrom.append(res[0][0])
            new_pos.append(res[0][1] + 1 if pos_is_one_based else res[0][1])
            keep.append(True)
        else:
            new_chrom.append(None)
            new_pos.append(None)
            keep.append(False)
    out = df.copy()
    out[f"{chrom_col}_hg19"] = out[chrom_col]
    out[f"{pos_col}_hg19"] = out[pos_col]
    out[chrom_col] = new_chrom
    out[pos_col] = new_pos
    n_dropped = int((~pd.Series(keep)).sum())
    return out.loc[keep].reset_index(drop=True), n_dropped


def lift_bed(cfg, df: pd.DataFrame, chrom_col="chrom", start_col="start", end_col="end"):
    """Lift interval starts and ends (0-based, BED convention); drop intervals whose ends
    don't both lift onto the same chromosome. Returns (lifted_df, n_dropped)."""
    s, ds = lift_positions(cfg, df, chrom_col, start_col, pos_is_one_based=False)
    e, _ = lift_positions(cfg, s, chrom_col, end_col, pos_is_one_based=False)
    good = e[end_col] > e[start_col]
    return e.loc[good].reset_index(drop=True), int(ds + (~good).sum())


def lift_finemap(cfg):
    """Lift each UKBB lipid fine-mapping table hg19 -> hg38 in place of config's expected path.

    Reads the raw release table, lifts, writes to the path config.yaml points variant_scoring at.
    Column names come from variant_scoring.COLS; adjust there if the release differs.
    """
    from . import variant_scoring as vs
    total_drop = 0
    for t in cfg["traits"]:
        dest = REPO_ROOT / t["finemap"]
        raw = dest.with_suffix(dest.suffix + ".hg19")   # expect the raw hg19 table saved here
        df = pd.read_csv(raw, sep="\t")
        lifted, dropped = lift_positions(cfg, df, vs.COLS["chrom"], vs.COLS["pos"],
                                         pos_is_one_based=True)
        dest.parent.mkdir(parents=True, exist_ok=True)
        lifted.to_csv(dest, sep="\t", index=False)
        total_drop += dropped
        print(f"  {t['name']}: lifted {len(lifted)}, dropped {dropped} -> {dest}")
    print(f"fine-mapping liftOver total dropped: {total_drop}")


if __name__ == "__main__":
    cfg = load_config()
    chain = REPO_ROOT / cfg["genome"]["liftover_chain"]
    if not chain.exists():
        raise SystemExit(f"chain MISSING (run: python -m src.download genome): {chain}")
    # CAGI5 and ABC can be lifted with lift_positions / lift_bed inside their own steps; the
    # fine-mapping tables are the ones the scoring loop reads directly, so lift them here.
    lift_finemap(cfg)
