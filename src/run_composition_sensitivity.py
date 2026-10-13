#!/usr/bin/env python3
"""Additional bulk-composition covariate checks around the locked identity score.

Does not replace the primary four tests.
Does not use IDENTITY_CORE_v1 genes as the abundance covariate.
Marker means are not cell-type deconvolution.
"""
from __future__ import annotations

import json

import pandas as pd
from scipy import stats

from locked_scoring import (
    GSE50244_ID_MAP,
    load_locked_sets,
    load_matrices,
    log,
    ols_fit,
    prepare_discovery,
    prepare_validation,
    score_sets,
)
from paths import DATA, EXPECTED, RESULTS, SETS

COMP = DATA / "locked_gene_sets" / "composition_proxy_sets_v1.json"
OUT = RESULTS / "composition_sensitivity_v1"


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


def spearman(df: pd.DataFrame, a: str, b: str) -> dict:
    work = df[[a, b]].dropna()
    if len(work) < 8:
        return {"n": len(work), "error": "too_few"}
    rho, p = stats.spearmanr(work[a], work[b])
    return {"n": int(len(work)), "rho": float(rho), "p": float(p)}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    locked = load_locked_sets()
    proxies = json.loads(COMP.read_text())
    identity = locked["IDENTITY_CORE_v1"]
    stress = locked["STRESS_CORE_v1"]
    proxy_sets = {k: v["genes"] for k, v in proxies["sets"].items()}
    overlap_id = sorted(set().union(*proxy_sets.values()) & set(identity))
    overlap_st = sorted(set().union(*proxy_sets.values()) & set(stress))
    if overlap_id or overlap_st:
        raise SystemExit(f"circular overlap identity={overlap_id} stress={overlap_st}")

    gene_sets = {
        "IDENTITY_CORE_v1": identity,
        "STRESS_CORE_v1": stress,
        **proxy_sets,
    }
    data = load_matrices(gene_sets)
    s502, cov502 = score_sets(data["e502"], gene_sets, GSE50244_ID_MAP)
    s768, cov768 = score_sets(data["e768"], gene_sets)
    cov502.insert(0, "series", "GSE50244")
    cov768.insert(0, "series", "GSE76894")
    coverage = pd.concat([cov502, cov768], ignore_index=True)
    coverage.to_csv(OUT / "gene_coverage.csv", index=False)

    d502 = prepare_discovery(data["p502"], s502)
    d768, batch_cols = prepare_validation(data["p768"], s768)
    d502.to_csv(OUT / "GSE50244_composition_scores.csv", index=False)
    d768.to_csv(OUT / "GSE76894_composition_scores.csv", index=False)

    proxy_names = list(proxy_sets)
    rows = []
    for proxy in proxy_names:
        sp = spearman(d502, proxy, "hba1c_n")
        rows.append({"analysis": "proxy_vs_endpoint", "series": "GSE50244", "outcome": proxy, "model": "spearman_hba1c", **{f"stat_{k}": v for k, v in sp.items()}})
        fit = ols_fit(d502, proxy, ["hba1c_n", "age_n", "sex_m", "bmi_n"])
        rows.append({"analysis": "proxy_vs_endpoint", "series": "GSE50244", "outcome": proxy, "model": "ols_hba1c_age_sex_bmi", **focus(fit, "hba1c_n")})
        fitv = ols_fit(d768, proxy, ["t2d", "age_n", "sex_m", "bmi_n"] + batch_cols)
        rows.append({"analysis": "proxy_vs_endpoint", "series": "GSE76894", "outcome": proxy, "model": "ols_t2d_age_sex_bmi_batch", **focus(fitv, "t2d")})
        sp_id = spearman(d502, "IDENTITY_CORE_v1", proxy)
        rows.append({"analysis": "identity_vs_proxy", "series": "GSE50244", "outcome": "IDENTITY_CORE_v1", "model": f"spearman_{proxy}", **{f"stat_{k}": v for k, v in sp_id.items()}})
        sp_idv = spearman(d768, "IDENTITY_CORE_v1", proxy)
        rows.append({"analysis": "identity_vs_proxy", "series": "GSE76894", "outcome": "IDENTITY_CORE_v1", "model": f"spearman_{proxy}", **{f"stat_{k}": v for k, v in sp_idv.items()}})

    base502 = ["hba1c_n", "age_n", "sex_m", "bmi_n"]
    base768 = ["t2d", "age_n", "sex_m", "bmi_n"] + batch_cols
    rows.append({"analysis": "identity_adjusted", "series": "GSE50244", "outcome": "IDENTITY_CORE_v1", "model": "ols_hba1c_age_sex_bmi", "covariate": "none", **focus(ols_fit(d502, "IDENTITY_CORE_v1", base502), "hba1c_n")})
    rows.append({"analysis": "identity_adjusted", "series": "GSE76894", "outcome": "IDENTITY_CORE_v1", "model": "ols_t2d_age_sex_bmi_batch", "covariate": "none", **focus(ols_fit(d768, "IDENTITY_CORE_v1", base768), "t2d")})
    for proxy in proxy_names:
        rows.append({"analysis": "identity_adjusted", "series": "GSE50244", "outcome": "IDENTITY_CORE_v1", "model": f"ols_hba1c_age_sex_bmi_{proxy}", "covariate": proxy, **focus(ols_fit(d502, "IDENTITY_CORE_v1", base502 + [proxy]), "hba1c_n")})
        rows.append({"analysis": "identity_adjusted", "series": "GSE76894", "outcome": "IDENTITY_CORE_v1", "model": f"ols_t2d_age_sex_bmi_batch_{proxy}", "covariate": proxy, **focus(ols_fit(d768, "IDENTITY_CORE_v1", base768 + [proxy]), "t2d")})
        rows.append({"analysis": "stress_adjusted", "series": "GSE50244", "outcome": "STRESS_CORE_v1", "model": f"ols_hba1c_age_sex_bmi_{proxy}", "covariate": proxy, **focus(ols_fit(d502, "STRESS_CORE_v1", base502 + [proxy]), "hba1c_n")})
        rows.append({"analysis": "stress_adjusted", "series": "GSE76894", "outcome": "STRESS_CORE_v1", "model": f"ols_t2d_age_sex_bmi_batch_{proxy}", "covariate": proxy, **focus(ols_fit(d768, "STRESS_CORE_v1", base768 + [proxy]), "t2d")})
    rows.append({"analysis": "identity_adjusted", "series": "GSE50244", "outcome": "IDENTITY_CORE_v1", "model": "ols_hba1c_age_sex_bmi_all_proxies", "covariate": "+".join(proxy_names), **focus(ols_fit(d502, "IDENTITY_CORE_v1", base502 + proxy_names), "hba1c_n")})
    rows.append({"analysis": "identity_adjusted", "series": "GSE76894", "outcome": "IDENTITY_CORE_v1", "model": "ols_t2d_age_sex_bmi_batch_all_proxies", "covariate": "+".join(proxy_names), **focus(ols_fit(d768, "IDENTITY_CORE_v1", base768 + proxy_names), "t2d")})

    rdf = pd.DataFrame(rows)
    rdf.to_csv(OUT / "composition_models.csv", index=False)
    compact = rdf.loc[rdf["analysis"] == "identity_adjusted", ["series", "covariate", "n", "beta", "ci95_low", "ci95_high", "p", "r2"]].copy()
    compact.to_csv(OUT / "composition_identity_compact.csv", index=False)
    notes = {
        "role": "additional_composition_covariates",
        "replaces_primary_four_tests": False,
        "circular_overlap_identity": overlap_id,
        "circular_overlap_stress": overlap_st,
        "beta_proxy_genes": proxy_sets["BETA_ABUNDANCE_PROXY_v1"],
        "exocrine_proxy_genes": proxy_sets["EXOCRINE_CONTAMINATION_PROXY_v1"],
        "endocrine_nonbeta_genes": proxy_sets["ENDOCRINE_NONBETA_PROXY_v1"],
        "note": "Marker-mean proxies, not deconvolution. IDENTITY_CORE_v1 genes were not used as abundance covariates.",
        "sets_file": str(SETS),
        "expected_snapshot": str(EXPECTED / "composition_sensitivity_v1"),
    }
    (OUT / "composition_notes.json").write_text(json.dumps(notes, indent=2))
    log(coverage.to_string(index=False))
    log(compact.to_string(index=False))
    log(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
