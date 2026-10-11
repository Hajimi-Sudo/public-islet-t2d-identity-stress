#!/usr/bin/env python3
"""Primary locked-score reconstruction (Gate-2/3)."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from locked_scoring import (
    ALL_SETS,
    GSE50244_ID_MAP,
    family_of,
    load_locked_sets,
    load_matrices,
    log,
    ols_fit,
    prepare_companion,
    prepare_discovery,
    prepare_orthogonal,
    prepare_validation,
    score_sets,
    spearman,
    welch,
)
from paths import LOCKED_OUT


def main() -> None:
    gene_sets = load_locked_sets()
    data = load_matrices(gene_sets)
    LOCKED_OUT.mkdir(parents=True, exist_ok=True)
    Path(LOCKED_OUT / "probe_map_sizes.json").write_text(json.dumps(data["mapping_provenance"], indent=2))
    log(json.dumps(data["mapping_provenance"], indent=2))

    log("Scoring")
    s502, c502 = score_sets(data["e502"], gene_sets, GSE50244_ID_MAP)
    s768, c768 = score_sets(data["e768"], gene_sets)
    cov_frames = []
    for name, cov in [("GSE50244", c502), ("GSE76894", c768)]:
        cov.insert(0, "series", name)
        cov_frames.append(cov)
    s386 = s209 = None
    if data["e386"] is not None:
        s386, c386 = score_sets(data["e386"], gene_sets)
        c386.insert(0, "series", "GSE38642")
        cov_frames.append(c386)
    if data["e209"] is not None:
        s209, c209 = score_sets(data["e209"], gene_sets)
        c209.insert(0, "series", "GSE20966")
        cov_frames.append(c209)
    coverage = pd.concat(cov_frames, ignore_index=True)
    coverage.to_csv(LOCKED_OUT / "gene_coverage.csv", index=False)

    d502 = prepare_discovery(data["p502"], s502)
    d768, batch_cols = prepare_validation(data["p768"], s768)
    d386 = prepare_companion(data["p386"], s386) if s386 is not None else None
    d209 = prepare_orthogonal(data["p209"], s209) if s209 is not None else None

    for name, df in [("GSE50244", d502), ("GSE38642", d386), ("GSE76894", d768), ("GSE20966", d209)]:
        if df is not None:
            df.to_csv(LOCKED_OUT / f"{name}_scores.csv", index=False)

    models = []
    for sname in ALL_SETS:
        family = family_of(sname)
        models.append({"series": "GSE50244", "set": sname, "family": family, "model": "spearman_hba1c", **{f"stat_{k}": v for k, v in spearman(d502, sname, "hba1c_n").items()}})
        models.append({"series": "GSE50244", "set": sname, "family": family, "model": "ols_hba1c_age_sex_bmi", "focus": "hba1c_n", **{f"stat_{k}": v for k, v in ols_fit(d502, sname, ["hba1c_n", "age_n", "sex_m", "bmi_n"]).items()}})
        models.append({"series": "GSE50244", "set": sname, "family": family, "model": "ols_hba1c65_age_sex_bmi", "focus": "hba1c65", **{f"stat_{k}": v for k, v in ols_fit(d502, sname, ["hba1c65", "age_n", "sex_m", "bmi_n"]).items()}})

    for sname in ALL_SETS:
        family = family_of(sname)
        models.append({"series": "GSE76894", "set": sname, "family": family, "model": "welch_t2d_unadjusted", **{f"stat_{k}": v for k, v in welch(d768, sname, "t2d_label", "T2D", "ND").items()}})
        models.append({"series": "GSE76894", "set": sname, "family": family, "model": "ols_t2d_age_sex_bmi_batch", "focus": "t2d", **{f"stat_{k}": v for k, v in ols_fit(d768, sname, ["t2d", "age_n", "sex_m", "bmi_n"] + batch_cols).items()}})
        models.append({"series": "GSE76894", "set": sname, "family": family, "model": "ols_t2d_age_sex_bmi", "focus": "t2d", **{f"stat_{k}": v for k, v in ols_fit(d768, sname, ["t2d", "age_n", "sex_m", "bmi_n"]).items()}})

    if d386 is not None:
        for sname in ALL_SETS:
            family = family_of(sname)
            models.append({"series": "GSE38642", "set": sname, "family": family, "model": "welch_t2d_unadjusted_companion", **{f"stat_{k}": v for k, v in welch(d386, sname, "t2d_label", "T2D", "ND").items()}})
            models.append({"series": "GSE38642", "set": sname, "family": family, "model": "ols_t2d_age_sex_bmi_companion", "focus": "t2d", **{f"stat_{k}": v for k, v in ols_fit(d386, sname, ["t2d", "age_n", "sex_m", "bmi_n"]).items()}})
    if d209 is not None:
        for sname in ALL_SETS:
            family = family_of(sname)
            models.append({"series": "GSE20966", "set": sname, "family": family, "model": "welch_t2d_unadjusted_orthogonal", **{f"stat_{k}": v for k, v in welch(d209, sname, "t2d_label", "T2D", "ND").items()}})
            models.append({"series": "GSE20966", "set": sname, "family": family, "model": "ols_t2d_age_sex_bmi_orthogonal", "focus": "t2d", **{f"stat_{k}": v for k, v in ols_fit(d209, sname, ["t2d", "age_n", "sex_m", "bmi_n"]).items()}})

    mdf = pd.DataFrame(models)
    mdf.to_csv(LOCKED_OUT / "model_results.csv", index=False)
    mdf.to_json(LOCKED_OUT / "model_results.json", orient="records", indent=2)
    primary_rows = mdf[
        (mdf["family"] == "primary")
        & (
            ((mdf["series"] == "GSE50244") & (mdf["model"] == "ols_hba1c_age_sex_bmi"))
            | ((mdf["series"] == "GSE76894") & (mdf["model"] == "ols_t2d_age_sex_bmi_batch"))
        )
    ].copy()
    primary_rows.to_csv(LOCKED_OUT / "primary_four_tests.csv", index=False)
    log(f"Wrote {LOCKED_OUT}")
    log(coverage.to_string(index=False))
    show_cols = [c for c in mdf.columns if c in {"series", "set", "family", "model", "stat_n", "stat_rho", "stat_p", "stat_beta_hba1c_n", "stat_p_hba1c_n", "stat_beta_t2d", "stat_p_t2d", "stat_diff"}]
    log(mdf[show_cols].to_string(index=False))


if __name__ == "__main__":
    main()
