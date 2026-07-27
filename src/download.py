"""Fetch models and data: HepG2 ChromBPNet (DNase + ATAC), UKB lipid fine-mapping,
CAGI5 MPRA, ABC regions, GENCODE/genome. Fill TODO accessions.

Usage:
    python -m src.download all
    python -m src.download models
"""
import sys

from .config import load_config, REPO_ROOT


def _run(cmd):
    import subprocess
    print(f"$ {cmd}")
    subprocess.run(cmd, shell=True, check=True)


def genome(cfg):
    out = REPO_ROOT / cfg["genome"]["fasta"]
    out.parent.mkdir(parents=True, exist_ok=True)
    _run(f"wget -O {out}.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz")
    _run(f"gunzip -f {out}.gz && samtools faidx {out}")
    # liftOver chain for the hg19 inputs (CAGI5, UKBB fine-mapping, ABC regions)
    chain = REPO_ROOT / cfg["genome"]["liftover_chain"]
    _run(f"wget -O {chain} https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/"
         f"hg19ToHg38.over.chain.gz")


def models(cfg):
    """Released HepG2 DNase + ATAC ChromBPNet models (ENCODE portal, hg38, public).

    Both assays are required — the assay-concordance analysis is the novel piece.
      DNase: ENCSR006CUK -> ENCFF615AKY.tar   (trained on DNase-seq ENCSR149XIL)
      ATAC:  ENCSR380YGX -> ENCFF137WCM.tar   (trained on ATAC-seq  ENCSR291GJU)
    """
    for assay in ("dnase", "atac"):
        d = REPO_ROOT / cfg["models"][assay]["dir"]
        d.mkdir(parents=True, exist_ok=True)
        url = cfg["models"][assay]["url"]
        tar = d / f"{cfg['models'][assay]['model_file']}.tar"
        _run(f"wget -O {tar} {url}")
        _run(f"tar -xf {tar} -C {d}")


def finemap(cfg):
    """UKBB 94-traits fine-mapping release 1.1 (Finucane lab); lipids are 3 of the 94. hg19.

    One tarball. After extracting, split/point the per-trait tables at the paths in config.yaml
    (trait_id LDLC/HDLC/TG). Variants are hg19 -> lift to hg38 before scoring (see liftover()).
    """
    dest = REPO_ROOT / "data/raw/finemap"
    dest.mkdir(parents=True, exist_ok=True)
    tar = dest / "UKBB_94traits_release1.1.tar.gz"
    # Public Dropbox share (finucanelab.org/data). dl=1 forces a direct download.
    _run(f"wget -O {tar} 'https://www.dropbox.com/s/cdsdgwxkxkcq8cn/"
         f"UKBB_94traits_release1.1.tar.gz?dl=1'")
    _run(f"tar -xzf {tar} -C {dest}")
    print("Extracted. Confirm lipid trait strings (LDLC/HDLC/TG) against the release manifest.")


def cagi5(cfg):
    """CAGI5 saturation-mutagenesis MPRA (Kircher 2019). SORT1/LDLR/F9 assayed in HepG2.

    The satMutMPRA portal is an interactive element selector, so this is not a single wget:
    go to https://kircherlab.bihealth.org/satMutMPRA/ , select SORT1 + LDLR + F9, choose the
    build in cfg['cagi5']['mpra_build'], and 'Download Selected Elements' to the path below.
    GEO mirror for scripted access: GSE126550.
    """
    out = REPO_ROOT / cfg["cagi5"]["mpra"]
    out.parent.mkdir(parents=True, exist_ok=True)
    if not out.exists():
        raise NotImplementedError(
            f"Manual step: download SORT1/LDLR/F9 from the satMutMPRA portal to {out} "
            f"(or fetch GSE126550 from GEO).")


def abc(cfg):
    """Engreitz-lab ABC predictions (Nasser 2021). Filter CellType == HepG2-Roadmap. hg19."""
    out = REPO_ROOT / cfg["regions"]["abc"]
    out.parent.mkdir(parents=True, exist_ok=True)
    _run(f"wget -O {out} {cfg['regions']['abc_url']}")
    print(f"Filter CellType == '{cfg['regions']['abc_celltype']}' downstream; lift hg19 -> hg38.")


ASSETS = {
    "genome": genome, "models": models, "finemap": finemap,
    "cagi5": cagi5, "abc": abc,
}


def main(argv):
    cfg = load_config()
    which = argv[1] if len(argv) > 1 else "all"
    targets = ASSETS if which == "all" else {which: ASSETS[which]}
    for name, fn in targets.items():
        print(f"\n=== {name} ===")
        fn(cfg)


if __name__ == "__main__":
    main(sys.argv)
