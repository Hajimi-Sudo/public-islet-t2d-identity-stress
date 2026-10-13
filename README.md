# Locked β-cell identity and stress scores in public human islet transcriptomes

Reproducible analysis code for the manuscript:

**Discordant reconstruction of locked β-cell identity and stress programs in public human islet transcriptomes**

Target journal: *Frontiers in Endocrinology* (Diabetes: Molecular Mechanisms).

This repository reconstructs two **literature-locked** gene-set scores in public human-islet transcriptomes:

- `IDENTITY_CORE_v1` (21 genes): expected **down** with worse glycemia / recorded T2D
- `STRESS_CORE_v1` (18 genes): expected **up** in the same models

The locked dual-program protocol is kept even though stress-up is **not** supported in the primary tests. Gene membership is frozen before any T2D contrast. No gene is added or dropped from differential expression.

Code repository: https://github.com/Hajimi-Sudo/public-islet-t2d-identity-stress

## What is in this repository

| Path | Contents |
|---|---|
| `data/locked_gene_sets/` | Frozen identity/stress sets plus additional composition proxies (`composition_proxy_sets_v1.json`) |
| `data/phenotypes/` | Donor-level GEO phenotypes used in the paper |
| `data/overlays/` | Frozen Open Targets MONDO_0005148 top-500 snapshot |
| `src/` | Download, scoring, overlay, sensitivity, and figure scripts |
| `expected_results/` | Frozen numerical snapshots from the locked run, including `sensitivity_v1/` and `composition_sensitivity_v1/` |
| `manuscript_figures/` | Frozen Figure 1–4 PDFs used in the manuscript |

## What is **not** in git

Large public matrices are downloaded into `cache/` (or `$PUBLIC_ISLET_CACHE`) and are **not** shipped:

- GEO series matrices and the Fadista GSE50244 processed gene matrix
- NCBI GPL `*.annot.gz` files
- GTEx v8 gene median TPM GCT

Also excluded on purpose:

- `GSE50244_phenotype_with_transferred_t2d.csv` (GSE38642 labels must not be transferred onto GSE50244)
- Fallback probe-map TSVs
- FASTQ / RAW.tar / NHANES / MIMIC / GWAS Catalog dumps / two-sample MR

## Cohort roles (do not reverse)

| Series | Role | Endpoint |
|---|---|---|
| GSE50244 | Discovery | continuous HbA1c (n=77 complete); **no** GEO T2D field |
| GSE76894 | Independent validation | recorded T2D vs ND, adjusted for age, sex, BMI, **batch** |
| GSE38642 | Same-biobank companion | 53/89 donor overlap with GSE50244; **not** independent |
| GSE20966 | Orthogonal LCM check | 10 vs 10; underpowered |

GSE76894 identity coverage is **16/21** (missing `PDX1`, `MAFA`, `KCNJ11`, `UCN3`, `GJD2`). That missingness is documented; it is not a reason to rewrite the locked set.

## Primary four tests (locked family)

From `expected_results/locked_scores_v1/primary_four_tests.csv`:

| Series | Set | Model | n | β | p | Verdict |
|---|---|---|---:|---:|---:|---|
| GSE50244 | IDENTITY | score ~ HbA1c + age + sex + BMI | 77 | −0.200 | 0.035 | supported (down) |
| GSE50244 | STRESS | same | 77 | +0.031 | 0.71 | **not** up |
| GSE76894 | IDENTITY | score ~ T2D + age + sex + BMI + batch | 102 | −0.569 | 0.015 | supported (down) |
| GSE76894 | STRESS | same | 102 | −0.181 | 0.19 | **not** up |

Joint identity-down **and** stress-up language is blocked. Sensitivity analyses (drop-5 identity genes, drop INS, leave-one-gene-out, drop batch 6, 500 outcome permutations with seed 66) must not replace these four tests.

Additional bulk-composition checks score three short marker means that **do not overlap** `IDENTITY_CORE_v1`: a 4-gene β-cell abundance proxy (`G6PC2`, `INSM1`, `SCGN`, `SYT13`), an exocrine contamination proxy, and a non-β endocrine proxy. These are not deconvolution. On these matrices the β-cell abundance proxy is collinear with the identity score, so reduced identity scores may partly reflect bulk β-cell abundance.

## Requirements

Python 3.10+ and the packages in `requirements.txt`:

```bash
python3 -m pip install -r requirements.txt
```

No GPU is required. Runtime on a laptop CPU is typically a few minutes after the public files are cached.

## Quick start

```bash
# 1. Download public matrices into ./cache
python3 src/download_public_data.py

# 2. Locked scores + primary models
python3 src/run_locked_scores.py

# 3. Tissue-selectivity overlay (not a T2D test)
python3 src/run_gtex_overlay.py

# 4. Open Targets overlay (not MR)
python3 src/run_opentargets_overlay.py

# 5. Prespecified sensitivity (optional; do not promote over step 2)
python3 src/run_sensitivity.py
# smoke-test only: PUBLIC_ISLET_N_PERM=5 python3 src/run_sensitivity.py
# paper default is 500 permutations, seed 66

# 5b. Additional composition covariates (optional; not deconvolution)
python3 src/run_composition_sensitivity.py

# 6. Figures/tables into results/
python3 src/make_figures.py

# 7. Compare the four primary coefficients to expected_results/
python3 src/compare_to_expected.py
```

Or:

```bash
bash run_all.sh
```

If the GEO/GTEx files already live somewhere else (for example a previous analysis cache):

```bash
export PUBLIC_ISLET_CACHE=/path/to/cache   # must contain geo/, annot/, gtex/
bash run_all.sh
```

`download_public_data.py` skips a file when the local size matches the expected NCBI/GTEx byte count. If a fetch fails under a proxy, it retries once with proxy environment variables unset.

Outputs are written to `results/` (`locked_scores_v1/`, `gtex_overlay_v1/`, `gwas_overlay_v1/`, `sensitivity_v1/`, `composition_sensitivity_v1/`, `figures/`, `tables/`).

## Scoring rule

For each locked set, genes present in the matrix are row-z-scored across donors, then averaged. GSE50244 maps `INS` → `IGF2+INS-IGF2+INS` in the processed Fadista matrix. Array series collapse probes to genes by median of probes annotated to that symbol.

## Public sources

- GEO FTP series matrices: GSE50244, GSE38642, GSE76894, GSE20966
- GSE50244 processed matrix: `GSE50244_Genes_counts_TMM_NormLength_atLeastMAF5_expressed.txt.gz`
- NCBI GEO platform annot: GPL570, GPL6244, GPL1352
- GTEx Analysis v8 gene median TPM (bulk pancreas, **not** purified islets)
- Open Targets GraphQL snapshot for MONDO_0005148 (type 2 diabetes mellitus), top 500 targets

Exact URLs and expected byte sizes are in `src/download_public_data.py`. Phenotype provenance is in `data/phenotypes/PHENOTYPE_SOURCE.md`.

## License

MIT. See `LICENSE`.

## Citation

See `CITATION.cff`. Cite the manuscript when it is published; until then cite this analysis code and the original GEO/GTEx/Open Targets sources.

## Repository

Published at https://github.com/Hajimi-Sudo/public-islet-t2d-identity-stress.

Do **not** commit `results/` or `cache/` (already listed in `.gitignore`). Do not add GEO matrices, GPL annot files, GTEx GCT, NHANES, MIMIC, or transferred-T2D phenotypes.

