#!/usr/bin/env bash
# End-to-end pipeline. Fill config.yaml accessions first. Steps are ordered so the
# correctness check (CAGI5) runs before anything downstream is trusted.
set -euo pipefail
cd "$(dirname "$0")/.."

echo "== 0. Download models + data (accessions are filled in config.yaml) =="
python -m src.download all
# NOTE: CAGI5 is a manual portal export (see src/download.py); everything else is scripted.

echo "== 0b. liftOver hg19 fine-mapping tables -> hg38 to match the models =="
python -m src.liftover              # lifts UKBB lipid tables; logs drop counts
# (CAGI5 is lifted inside cagi5_repro if you exported hg19; ABC regions inside ldsc use hg19.)

echo "== 1. Correctness: reproduce SORT1 / LDLR CAGI5 (HepG2 DNase) =="
python -m src.cagi5_repro          # -> results/tables/cagi5_reproduction.csv  (check vs Table 4)

echo "== 2. Cell-type justification: LDSC lipid heritability enrichment in HepG2 =="
python -m src.ldsc_enrichment      # -> results/tables/ldsc_enrichment.csv

echo "== 3. Score fine-mapped lipid variants, both assays =="
python -m src.variant_scoring dnase
python -m src.variant_scoring atac

echo "== 4. Concordance with fine-mapping (PIP enrichment) =="
python -m src.pip_enrichment       # -> results/tables/pip_enrichment.csv

echo "== 5. NEW: DNase vs ATAC counterfactual agreement =="
python -m src.assay_concordance    # -> results/tables/assay_disagreers.csv

echo "== 6. Mechanistic coherence: DeepSHAP + HepG2 motif check at top hits =="
python -m src.contributions

echo "Done. Build figures in notebooks/analysis.ipynb (src/plots.py)."
