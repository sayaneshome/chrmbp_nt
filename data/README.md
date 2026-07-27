# Data provenance

All large files are gitignored; fetch with `python -m src.download <asset>`. Record what you
actually pulled here so the analysis is reproducible.

| Asset | config key | Source / accession | Build | Notes |
|---|---|---|---|---|
| hg38 fasta | `genome.fasta` | UCSC goldenPath hg38 | hg38 | |
| hg19→hg38 chain | `genome.liftover_chain` | UCSC goldenPath liftOver | — | for the hg19 assets below |
| HepG2 DNase model | `models.dnase` | ENCODE ENCSR006CUK → ENCFF615AKY.tar | hg38 | tier-1; DNase ENCSR149XIL |
| HepG2 ATAC model | `models.atac` | ENCODE ENCSR380YGX → ENCFF137WCM.tar | hg38 | tier-1; ATAC ENCSR291GJU |
| LDL-C / HDL-C / TG fine-mapping | `traits[*]` | UKBB 94-traits release 1.1 (Finucane) | hg19 | one tar; PIP per variant; lift to hg38 |
| CAGI5 MPRA | `cagi5.mpra` | satMutMPRA portal / GEO GSE126550 | hg19 or hg38 | SORT1/LDLR/F9; manual portal export |
| ABC regions (HepG2) | `regions.abc` | Nasser 2021 AllPredictions...V3.txt.gz | hg19 | filter CellType==HepG2-Roadmap; lift |
| LDSC baseline/weights | `ldsc.*` | Alkes-group LDSC resources | hg19 | separate tool/env |
| Motif DB | `contributions.motif_db` | HOCOMOCO (or JASPAR) | — | for the HepG2-lexicon check |

**Build note:** models are hg38; CAGI5 / UKBB / ABC are hg19 → all lifted to hg38 before scoring.
Record any liftOver drop counts here when you run it.
