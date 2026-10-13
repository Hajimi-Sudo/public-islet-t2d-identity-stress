#!/usr/bin/env python3
"""Additional PERK/IRE1/ATF6 split of the locked stress/UPR core.

Does not replace the primary four tests.
Does not add genes from T2D differential expression.
"""
from __future__ import annotations

import json

import pandas as pd
from scipy import stats

from locked_scoring import (
    load_locked_sets,
    load_matrices,
    log,
    ols_fit,
    prepare_discovery,
    prepare_validation,
    score_sets,
    welch,
)
from paths import DATA, RESULTS

ARMS = DATA / "locked_gene_sets" / "upr_arm_sets_v1.json"
OUT = RESULTS / "upr_arm_split_v1"


def focus(fit: dict, term: str) -> dict:
    rec = {
        "n": fit.get("n"),
        "r2": fit.get("r2"),
        "error": fit.get("error"),
        "focus_term": term,
        "beta": fit.get(f"beta_{term}"),
        "se": fit.get(f"se_{term}"),
        "p": fit.get(f"p_{term}"),
        "ci95_low": fit.get(f"ci95_low_{term}"),
        "ci95_high": fit.get(f"ci95_high_{term}"),
    }
    for key, val in fit.items():
        if key.startswith(("beta_", "p_", "se_", "ci95_")):
            rec[key] = val
    return rec


def spearman(df, a, b):
    work = df[[a, b]].dropna()
    if len(work) < 8:
        return {"n": len(work), "error": "too_few"}
    rho, p = stats.spearmanr(work[a], work[b])
    return {"n": int(len(work)), "rho": float(rho), "p": float(p)}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    locked = load_locked_sets()
    arms = json.loads(ARMS.read_text())
    stress = locked["STRESS_CORE_v1"]
    arm_sets = {k: v["genes"] for k, v in arms["sets"].items()}
    extra = sorted(set().union(*arm_sets.values()) - set(stress))
    if extra:
        raise SystemExit(f"arm genes not in STRESS_CORE_v1: {extra}")
    gene_sets = {"STRESS_CORE_v1": stress, **arm_sets}
    data = load_matrices(gene_sets)
    s502, cov502 = score_sets(data["e502"], gene_sets)
    s768, cov768 = score_sets(data["e768"], gene_sets)
    cov502.insert(0, "series", "GSE50244")
    cov768.insert(0, "series", "GSE76894")
    coverage = pd.concat([cov502, cov768], ignore_index=True)
    coverage.to_csv(OUT / "gene_coverage.csv", index=False)
    d502 = prepare_discovery(data["p502"], s502)
    d768, batch_cols = prepare_validation(data["p768"], s768)
    d502.to_csv(OUT / "GSE50244_upr_arm_scores.csv", index=False)
    d768.to_csv(OUT / "GSE76894_upr_arm_scores.csv", index=False)
    outcomes = ["STRESS_CORE_v1"] + list(arm_sets)
    rows = []
    base502 = ["hba1c_n", "age_n", "sex_m", "bmi_n"]
    base768 = ["t2d", "age_n", "sex_m", "bmi_n"] + batch_cols
    for outcome in outcomes:
        sp = spearman(d502, outcome, "hba1c_n")
        rows.append({"analysis": "spearman_hba1c", "series": "GSE50244", "set": outcome, **{f"stat_{k}": v for k, v in sp.items()}})
        rows.append({"analysis": "primary_like", "series": "GSE50244", "set": outcome, "model": "ols_hba1c_age_sex_bmi", **focus(ols_fit(d502, outcome, base502), "hba1c_n")})
        w = welch(d768, outcome, "t2d_label", "T2D", "ND")
        rows.append({"analysis": "welch_t2d", "series": "GSE76894", "set": outcome, **{f"stat_{k}": v for k, v in w.items()}})
        rows.append({"analysis": "primary_like", "series": "GSE76894", "set": outcome, "model": "ols_t2d_age_sex_bmi_batch", **focus(ols_fit(d768, outcome, base768), "t2d")})
    rdf = pd.DataFrame(rows)
    rdf.to_csv(OUT / "upr_arm_models.csv", index=False)
    compact = rdf.loc[rdf["analysis"] == "primary_like", ["series", "set", "n", "beta", "ci95_low", "ci95_high", "p", "r2"]].copy()
    compact.to_csv(OUT / "upr_arm_compact.csv", index=False)
    notes = {
        "role": "additional_upr_arm_split",
        "replaces_primary_four_tests": False,
        "parent_set": "STRESS_CORE_v1",
        "arms": arm_sets,
        "shared_gene": "HSPA5",
        "note": "Arm scores are subsets of STRESS_CORE_v1. They do not replace the locked 18-gene stress mean.",
    }
    (OUT / "upr_arm_notes.json").write_text(json.dumps(notes, indent=2))
    log(coverage.to_string(index=False))
    log(compact.to_string(index=False))
    log(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
