#!/usr/bin/env python3
"""Repository paths. Large GEO/GTEx matrices live in cache/, not in git."""
from __future__ import annotations

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CACHE_ROOT = Path(os.environ.get("PUBLIC_ISLET_CACHE", REPO_ROOT / "cache")).expanduser().resolve()

DATA = REPO_ROOT / "data"
SETS = DATA / "locked_gene_sets" / "locked_gene_sets_v1.json"
PHENO = DATA / "phenotypes"
OVERLAYS = DATA / "overlays"

GEO = CACHE_ROOT / "geo"
ANNOT = CACHE_ROOT / "annot"
GTEX_DIR = CACHE_ROOT / "gtex"

RESULTS = REPO_ROOT / "results"
LOCKED_OUT = RESULTS / "locked_scores_v1"
GTEX_OUT = RESULTS / "gtex_overlay_v1"
OT_OUT = RESULTS / "gwas_overlay_v1"
SENS_OUT = RESULTS / "sensitivity_v1"
FIG_OUT = RESULTS / "figures"
TAB_OUT = RESULTS / "tables"

EXPECTED = REPO_ROOT / "expected_results"

GSE50244_ID_MAP = {"INS": "IGF2+INS-IGF2+INS"}
PRIMARY_SETS = ["IDENTITY_CORE_v1", "STRESS_CORE_v1"]
SECONDARY_SETS = ["HALLMARK_PANCREAS_BETA_CELLS", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"]
ALL_SETS = PRIMARY_SETS + SECONDARY_SETS
IDENTITY_MISSING_ON_GSE76894 = ["PDX1", "MAFA", "KCNJ11", "UCN3", "GJD2"]
