#!/usr/bin/env python3
"""Additional mean-versus-median scoring around the locked four tests.

Does not replace the primary four tests (mean of gene-wise z-scores).
Does not add genes from T2D differential expression.
"""
from __future__ import annotations

import json

import numpy as np
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
from paths import EXPECTED, RESULTS

OUT = RESULTS / "median_scoring_v1"


def focus(fit: dict, term: str) -> dict:
    return {
        "n": fit.get("n"),
        "beta": fit.get(f"beta_{term}"),
        "se": fit.get(f"se_{term}"),
        "p": fit.get(f"p_{term}"),
        "ci95_low": fit.get(f"ci95_low_{term}"),
        "ci95_high": fit.get(f"ci95_high_{term}"),
        "r2": fit.get("r2"),
        "error": fit.get("error"),
    }


def spearman(a: pd.Series, b: pd.Series) -> dict:
    work = pd.concat([a, b], axis=1).dropna()
    if len(work) < 8:
        return {"n": len(work), "error": "too_few"}
    rho, p = stats.spearmanr(work.iloc[:, 0], work.iloc[:, 1])
    return {"n": int(len(work)), "rho": float(rho), "p": float(p)}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    locked = load_locked_sets()
    identity = list(locked["IDENTITY_CORE_v1"])
    stress = list(locked["STRESS_CORE_v1"])
    gene_sets = {"IDENTITY_CORE_v1": identity, "STRESS_CORE_v1": stress}
    data = load_matrices(gene_sets)

    mean502, _ = score_sets(data["e502"], gene_sets, GSE50244_ID_MAP, reduce="mean")
    med502, _ = score_sets(data["e502"], gene_sets, GSE50244_ID_MAP, reduce="median")
    mean768, _ = score_sets(data["e768"], gene_sets, reduce="mean")
    med768, _ = score_sets(data["e768"], gene_sets, reduce="median")

    d502_mean = prepare_discovery(
        data["p502"], mean502.rename(columns={"IDENTITY_CORE_v1": "IDENTITY_mean", "STRESS_CORE_v1": "STRESS_mean"})
    )
    d502_med = prepare_discovery(
        data["p502"], med502.rename(columns={"IDENTITY_CORE_v1": "IDENTITY_median", "STRESS_CORE_v1": "STRESS_median"})
    )
    d502 = d502_mean.merge(d502_med[["geo_accession", "IDENTITY_median", "STRESS_median"]], on="geo_accession", how="left")
    d768_mean, batch_cols = prepare_validation(
        data["p768"], mean768.rename(columns={"IDENTITY_CORE_v1": "IDENTITY_mean", "STRESS_CORE_v1": "STRESS_mean"})
    )
    d768_med, _ = prepare_validation(
        data["p768"], med768.rename(columns={"IDENTITY_CORE_v1": "IDENTITY_median", "STRESS_CORE_v1": "STRESS_median"})
    )
    d768 = d768_mean.merge(d768_med[["geo_accession", "IDENTITY_median", "STRESS_median"]], on="geo_accession", how="left")

    rows = []
    corrs = []
    base502 = ["hba1c_n", "age_n", "sex_m", "bmi_n"]
    base768 = ["t2d", "age_n", "sex_m", "bmi_n"] + batch_cols
    for series, df, y_mean, y_med, term, model in [
        ("GSE50244", d502, "IDENTITY_mean", "IDENTITY_median", "hba1c_n", "ols_hba1c_age_sex_bmi"),
        ("GSE50244", d502, "STRESS_mean", "STRESS_median", "hba1c_n", "ols_hba1c_age_sex_bmi"),
        ("GSE76894", d768, "IDENTITY_mean", "IDENTITY_median", "t2d", "ols_t2d_age_sex_bmi_batch"),
        ("GSE76894", d768, "STRESS_mean", "STRESS_median", "t2d", "ols_t2d_age_sex_bmi_batch"),
    ]:
        set_name = "IDENTITY_CORE_v1" if "IDENTITY" in y_mean else "STRESS_CORE_v1"
        x = base502 if series == "GSE50244" else base768
        fit_mean = ols_fit(df.rename(columns={y_mean: "y"}), "y", x)
        fit_med = ols_fit(df.rename(columns={y_med: "y"}), "y", x)
        sp = spearman(df[y_mean], df[y_med])
        corrs.append({"series": series, "set": set_name, **sp})
        rows.append({"analysis": "mean_zscore", "series": series, "set": set_name, "reduce": "mean", "model": model, **focus(fit_mean, term), "spearman_mean_vs_median": sp.get("rho")})
        rows.append({"analysis": "median_zscore", "series": series, "set": set_name, "reduce": "median", "model": model, **focus(fit_med, term), "spearman_mean_vs_median": sp.get("rho")})

    no_iapp = [g for g in identity if g != "IAPP"]
    s_iapp, _ = score_sets(data["e502"], {"IDENTITY_CORE_v1": no_iapp}, GSE50244_ID_MAP, reduce="mean")
    d_iapp = prepare_discovery(data["p502"], s_iapp)
    fit_iapp = ols_fit(d_iapp, "IDENTITY_CORE_v1", base502)
    rows.append({"analysis": "drop_IAPP", "series": "GSE50244", "set": "IDENTITY_CORE_v1_minus_IAPP", "reduce": "mean", "model": "ols_hba1c_age_sex_bmi", **focus(fit_iapp, "hba1c_n"), "spearman_mean_vs_median": np.nan})

    no_ins = [g for g in identity if g != "INS"]
    s768_noins, _ = score_sets(data["e768"], {"IDENTITY_CORE_v1": no_ins}, reduce="mean")
    d768_noins, batch_cols2 = prepare_validation(data["p768"], s768_noins)
    fit768_noins = ols_fit(d768_noins, "IDENTITY_CORE_v1", ["t2d", "age_n", "sex_m", "bmi_n"] + batch_cols2)
    rows.append({"analysis": "drop_INS", "series": "GSE76894", "set": "IDENTITY_CORE_v1_minus_INS", "reduce": "mean", "model": "ols_t2d_age_sex_bmi_batch", **focus(fit768_noins, "t2d"), "spearman_mean_vs_median": np.nan})

    rdf = pd.DataFrame(rows)
    cdf = pd.DataFrame(corrs)
    rdf.to_csv(OUT / "median_scoring_models.csv", index=False)
    cdf.to_csv(OUT / "mean_vs_median_spearman.csv", index=False)
    compact = rdf.loc[rdf["analysis"].isin(["median_zscore", "drop_IAPP", "drop_INS"]), ["analysis", "series", "set", "n", "beta", "ci95_low", "ci95_high", "p", "spearman_mean_vs_median"]].copy()
    compact.to_csv(OUT / "median_scoring_compact.csv", index=False)
    notes = {
        "role": "additional_mean_vs_median_scoring",
        "replaces_primary_four_tests": False,
        "primary_reduce": "mean of gene-wise z-scores",
        "additional_reduce": "median of gene-wise z-scores",
        "also": ["drop_IAPP on GSE50244 identity", "drop_INS on GSE76894 identity"],
        "note": "Median scoring down-weights extreme genes such as INS. It does not replace the locked mean score.",
        "expected_snapshot": str(EXPECTED / "median_scoring_v1"),
    }
    (OUT / "median_scoring_notes.json").write_text(json.dumps(notes, indent=2))
    log(rdf.to_string(index=False))
    log(cdf.to_string(index=False))
    log(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
