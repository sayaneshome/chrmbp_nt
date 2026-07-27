# Scoring fine-mapped blood-lipid variants with released HepG2 ChromBPNet models

*Exploratory. One causal cell type (HepG2, hepatocyte model), released ATAC + DNase ChromBPNet
models, blood lipid traits (LDL-C / HDL-C / TG). A transplant of the ChromBPNet paper's
blood-trait analysis (Fig. 7e) to the hepatocyte / lipid axis, plus one thing the paper didn't
report.*

---

## Question

Non-coding GWAS variants for blood lipids act largely through hepatocyte regulatory sequence.
**Do released HepG2 ChromBPNet models, applied as counterfactual variant scorers, agree with
statistical fine-mapping for lipid traits — and do the DNase and ATAC models agree with each
other on the *variant* predictions?**

The cell type is not a proxy. Hepatocytes are the settled causal cell type for LDL cholesterol;
HepG2 is the standard hepatocyte model and one of the five ENCODE tier-1 lines with both ATAC and
DNase ChromBPNet models released. That removes the usual "is this the right cell type" argument
before it starts.

Two sub-questions, in order of increasing novelty:
1. **Correctness** — can we reproduce the paper's published HepG2 CAGI5 MPRA numbers at SORT1 /
   LDLR (Table 4) before doing anything new?
2. **Concordance** — do high ChromBPNet variant scores enrich for high-PIP fine-mapped lipid
   variants, above what regulatory annotation alone gives?
3. **New** — at fine-mapped variants, do the **DNase and ATAC** HepG2 models agree on the
   *counterfactual* effect? Reference-sequence agreement after bias correction is reported
   (JSD 0.81 → 0.26); agreement on variant effects is, as far as I can tell, not.

## Method

1. **Models** — released HepG2 ChromBPNet models, **both DNase and ATAC** (`models/hepg2_*`).
   No training. See `models/*/README.md` for provenance.
2. **Positive controls / correctness** — reproduce SORT1 (rs12740374) and LDLR CAGI5 MPRA
   correlations with the HepG2 DNase model; check against the paper's per-locus supplement values
   (main-text benchmark is **Fig. 7d**; per-locus r values are in the supplement, Zenodo
   `10.5281/zenodo.14735906`) before trusting anything downstream. `src/cagi5_repro.py`.
3. **Heritability enrichment (justifies the cell type quantitatively)** — LDSC partitioned
   heritability of lipid traits in HepG2 regulatory regions, mirroring the paper's K562 / blood
   step. If lipid heritability is enriched in HepG2 peaks, the cell-type choice is earned, not
   asserted. `src/ldsc_enrichment.py`.
4. **Variant scoring** — score fine-mapped lipid variants (UK Biobank fine-mapping, Finucane-lab
   resource — same resource the paper used for blood traits, different trait file) as ref-vs-alt
   counterfactuals: log fold-change of predicted counts + profile JSD, per assay.
   `src/variant_scoring.py`.
5. **PIP concordance** — enrichment of high-|score| variants against posterior-inclusion-probability
   thresholds, vs. matched background. `src/pip_enrichment.py`.
6. **Mechanistic coherence** — DeepSHAP contribution scores at top hits; check the disrupted motifs
   are the HepG2 lexicon (HNF4A, FOXA, CEBP, HNF1A — paper Fig. 5a). Enrichment that lands on the
   right TF motifs is coherence, not just correlation. `src/contributions.py`.
7. **Assay concordance (the new bit)** — DNase vs ATAC counterfactual scores at the same fine-mapped
   variants; flag disagreers (profile-shape vs coverage effects, or model instability).
   `src/assay_concordance.py`.

Reproduce with `scripts/run_all.sh` or step through `notebooks/analysis.ipynb`.

## One figure

`results/figures/pip_enrichment.png` — enrichment of high ChromBPNet variant scores across PIP bins
for lipid traits, with the SORT1 / LDLR positive controls marked; companion panel:
`results/figures/atac_vs_dnase_variant_scores.png` for the assay-concordance question.

*(Placeholders until the notebook is run.)*

## What surprised me

*(Fill in after running. Candidate: whether ATAC and DNase agree on reference accessibility but
diverge on specific variant effects, and whether the divergers are profile-shape variants.)*

## Limitations

- **HepG2 is a hepatoblastoma line, not primary hepatocytes.** Standard in the field, but it is not
  normal liver, and some regulatory wiring may be cancer-context specific.
- **Enrichment is concordance, not causality.** Agreement between model scores and statistical
  fine-mapping can arise because both track the same confounder (e.g. sequence conservation, GC).
  The LDSC step and the motif-coherence check are what keep the result from being circular.
- **One cell type, released models, no retraining** — an observation about where these specific
  models and lipid fine-mapping agree, not a benchmark of variant-effect prediction.
- **Genome-build mismatch is a real error source.** The HepG2 ChromBPNet models are hg38, but the
  CAGI5 MPRA, the UKBB fine-mapping release, and the ABC regions are hg19. Coordinates are lifted
  over (hg19→hg38) before scoring; any variant that fails to lift, or lifts ambiguously, is dropped
  and counted — check the drop log before reading the enrichment.

## Provenance quick-reference

| Asset | Source / ID | Build | Notes |
|---|---|---|---|
| HepG2 DNase ChromBPNet | ENCODE `ENCSR006CUK` → `ENCFF615AKY.tar` | hg38 | public; trained on DNase ENCSR149XIL |
| HepG2 ATAC ChromBPNet | ENCODE `ENCSR380YGX` → `ENCFF137WCM.tar` | hg38 | public; trained on ATAC ENCSR291GJU |
| SORT1 / LDLR / F9 controls | CAGI5 satMutMPRA portal; GEO `GSE126550` | hg19/hg38 | HepG2 loci; targets in Fig. 7d + suppl. |
| Fine-mapped lipid variants | UKBB 94-traits release 1.1 (Finucane lab) | hg19 | traits `LDLC` / `HDLC` / `TG` |
| Regulatory regions | ABC Nasser 2021 `AllPredictions...V3.txt.gz` | hg19 | filter `CellType == HepG2-Roadmap` |
| TF motif reference | paper Fig. 5a lexicon | — | HNF4A, FOXA, CEBP, HNF1A |
| Paper per-locus targets | Zenodo `10.5281/zenodo.14735906` | — | supplement to Fig. 7d benchmark |
| (alt project) SMC + CAD | ChromBPNet SMC model `syn59479965` | — | scATAC pseudobulk; Synapse login; weaker controls |

*All hg19 assets are lifted to hg38 to match the models (see Limitations).*

## Repo layout

```
chrombpnet-hepg2-lipids/
├── README.md
├── environment.yml
├── config.yaml               # traits, model paths, thresholds, region/annotation paths
├── data/{raw,processed}/     # models, fine-mapping, CAGI5, ABC regions (gitignored)
├── models/{hepg2_dnase,hepg2_atac}/
├── src/
│   ├── sequtils.py           # one-hot, window extraction, JSD, logFC — pure numpy (tested)
│   ├── liftover.py           # hg19 -> hg38 via pyliftover, with drop counts (tested-adjacent)
│   ├── chrombpnet_model.py   # Keras bias-corrected model wrapper (the only TF-only module)
│   ├── download.py           # models (both assays), UKB lipid fine-mapping, CAGI5, ABC regions
│   ├── variant_scoring.py    # ref-vs-alt counterfactual: logFC counts + profile JSD, per assay
│   ├── cagi5_repro.py        # SORT1 / LDLR CAGI5 correlation — correctness check
│   ├── ldsc_enrichment.py    # partitioned heritability of lipid traits in HepG2 regions
│   ├── pip_enrichment.py     # high-score enrichment vs PIP bins (Fisher OR, tested)
│   ├── contributions.py      # DeepSHAP at top hits + HepG2 motif coherence check
│   ├── assay_concordance.py  # DNase vs ATAC variant-effect agreement (novel)
│   └── plots.py              # figures (repo CLAUDE.md standards baked in)
├── tests/test_core.py        # framework-free core: seq utils, metrics, enrichment (7 tests, pass)
├── notebooks/analysis.ipynb
├── results/{figures,tables}/
└── scripts/run_all.sh
```

