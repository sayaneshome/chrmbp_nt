"""Partitioned-heritability enrichment of lipid traits in HepG2 regulatory regions (S-LDSC).

Mirrors the paper's K562/blood step, transplanted to HepG2/lipids. If lipid heritability is
enriched in HepG2 regions above the LDSC baseline, the cell-type choice is earned quantitatively.

LDSC (github.com/bulik/ldsc) is an external tool with its own (python-2/older) environment, so this
module ORCHESTRATES it via subprocess rather than importing it. Point cfg['ldsc']['ldsc_dir'] at a
checkout with make_annot.py / ldsc.py on PATH or in that dir. All LDSC inputs are hg19 — the HepG2
regions must be the hg19 (pre-liftover) BED, since the baseline/weights/plink files are hg19.
"""
from __future__ import annotations
from pathlib import Path
import subprocess
import pandas as pd

from .config import load_config, REPO_ROOT

CHROMS = list(range(1, 23))


def _sh(cmd):
    print(f"$ {cmd}")
    subprocess.run(cmd, shell=True, check=True)


def _ldsc(cfg, script):
    d = cfg["ldsc"].get("ldsc_dir", "ldsc")
    return f"python {Path(d)/script}"


def make_annotation(cfg):
    """HepG2 regions (hg19) -> per-chrom .annot.gz + LD scores against baseline SNPs."""
    regions_hg19 = REPO_ROOT / cfg["ldsc"]["hepg2_regions_hg19"]   # BED, HepG2 peaks/ABC in hg19
    annot_prefix = REPO_ROOT / cfg["ldsc"]["annot_prefix"]
    annot_prefix.parent.mkdir(parents=True, exist_ok=True)
    bim = cfg["ldsc"]["plink_bfile"]      # e.g. 1000G_EUR_Phase3_plink/1000G.EUR.QC.{chrom}
    for c in CHROMS:
        _sh(f"{_ldsc(cfg,'make_annot.py')} --bed-file {regions_hg19} "
            f"--bimfile {bim.format(chrom=c)}.bim --annot-file {annot_prefix}.{c}.annot.gz")
        _sh(f"{_ldsc(cfg,'ldsc.py')} --l2 --bfile {bim.format(chrom=c)} --ld-wind-cm 1 "
            f"--annot {annot_prefix}.{c}.annot.gz --thin-annot "
            f"--out {annot_prefix}.{c} --print-snps {cfg['ldsc']['hapmap_snps']}")


def partitioned_h2(cfg) -> pd.DataFrame:
    """ldsc.py --h2 per lipid trait against baseline + the HepG2 annotation; parse .results."""
    annot_prefix = REPO_ROOT / cfg["ldsc"]["annot_prefix"]
    baseline = cfg["ldsc"]["baseline"]       # baselineLD chr prefix
    weights = cfg["ldsc"]["weights"]
    freq = cfg["ldsc"]["frqfile"]
    rows = []
    for t in cfg["traits"]:
        sumstats = f"{cfg['ldsc']['sumstats_dir']}/{t['trait_id']}.sumstats.gz"
        out = REPO_ROOT / cfg["paths"]["tables"] / f"ldsc_{t['trait_id']}"
        _sh(f"{_ldsc(cfg,'ldsc.py')} --h2 {sumstats} "
            f"--ref-ld-chr {baseline},{annot_prefix}. --w-ld-chr {weights} "
            f"--frqfile-chr {freq} --overlap-annot --print-coefficients --out {out}")
        res = pd.read_csv(f"{out}.results", sep="\t")
        hepg2 = res.iloc[-1]      # the HepG2 annotation is the last category appended
        rows.append({"trait": t["name"],
                     "enrichment": hepg2.get("Enrichment"),
                     "enrichment_p": hepg2.get("Enrichment_p"),
                     "coefficient_z": hepg2.get("Coefficient_z-score")})
    return pd.DataFrame(rows)


def main():
    cfg = load_config()
    make_annotation(cfg)
    res = partitioned_h2(cfg)
    out = REPO_ROOT / cfg["paths"]["tables"] / "ldsc_enrichment.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out, index=False)
    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
