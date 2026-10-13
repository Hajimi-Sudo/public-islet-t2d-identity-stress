#!/usr/bin/env python3
"""Prespecified robustness around the primary four tests.

Does not add or drop genes from T2D differential expression.
Does not reverse the locked dual-program family.
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd

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
from paths import IDENTITY_MISSING_ON_GSE76894, SENS_OUT

# Default 500 permutations. PUBLIC_ISLET_N_PERM can shorten a smoke test.
N_PERM = int(os.environ.get("PUBLIC_ISLET_N_PERM", "500"))
RNG_SEED = 66


def _focus_beta(fit: dict, term: str) -> dict:
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


def score_identity(expr, genes, id_map=None):
    scores, cov = score_sets(expr, {"IDENTITY_CORE_v1": genes}, id_map)
    return scores.rename(columns={"IDENTITY_CORE_v1": "score"}), cov.iloc[0].to_dict()


def permutation_p(observed: float, null: np.ndarray) -> float:
    if observed is None or np.isnan(observed) or len(null) == 0:
        return float("nan")
    return float((np.sum(np.abs(null) >= abs(observed)) + 1) / (len(null) + 1))


def main() -> None:
    SENS_OUT.mkdir(parents=True, exist_ok=True)
    gene_sets = load_locked_sets()
    identity = list(gene_sets["IDENTITY_CORE_v1"])
    data = load_matrices(gene_sets)

    s502, _ = score_sets(data["e502"], gene_sets, GSE50244_ID_MAP)
    s768, _ = score_sets(data["e768"], gene_sets)
    d502 = prepare_discovery(data["p502"], s502)
    d768, batch_cols = prepare_validation(data["p768"], s768)

    rows = []

    # 1. Drop the five identity genes missing on GSE76894, then rescore GSE50244.
    drop5 = [g for g in identity if g not in IDENTITY_MISSING_ON_GSE76894]
    s502_drop5, cov5 = score_identity(data["e502"], drop5, GSE50244_ID_MAP)
    d502_drop5 = prepare_discovery(data["p502"], s502_drop5.rename(columns={"score": "IDENTITY_CORE_v1"}))
    fit5 = ols_fit(d502_drop5, "IDENTITY_CORE_v1", ["hba1c_n", "age_n", "sex_m", "bmi_n"])
    rows.append(
        {
            "analysis": "drop5_identity_missing_on_GSE76894",
            "series": "GSE50244",
            "set": "IDENTITY_CORE_v1_minus_PDX1_MAFA_KCNJ11_UCN3_GJD2",
            "n_genes": cov5["n_present"],
            "missing": cov5["missing"],
            **_focus_beta(fit5, "hba1c_n"),
        }
    )

    # 2. INS-removed identity on GSE50244.
    no_ins = [g for g in identity if g != "INS"]
    s502_noins, cov_ins = score_identity(data["e502"], no_ins, GSE50244_ID_MAP)
    d502_noins = prepare_discovery(data["p502"], s502_noins.rename(columns={"score": "IDENTITY_CORE_v1"}))
    fit_ins = ols_fit(d502_noins, "IDENTITY_CORE_v1", ["hba1c_n", "age_n", "sex_m", "bmi_n"])
    rows.append(
        {
            "analysis": "drop_INS",
            "series": "GSE50244",
            "set": "IDENTITY_CORE_v1_minus_INS",
            "n_genes": cov_ins["n_present"],
            "missing": cov_ins["missing"],
            **_focus_beta(fit_ins, "hba1c_n"),
        }
    )

    # 3. Leave-one-gene-out identity on GSE50244.
    loo_rows = []
    for g in identity:
        keep = [x for x in identity if x != g]
        s_loo, cov_loo = score_identity(data["e502"], keep, GSE50244_ID_MAP)
        d_loo = prepare_discovery(data["p502"], s_loo.rename(columns={"score": "IDENTITY_CORE_v1"}))
        fit_loo = ols_fit(d_loo, "IDENTITY_CORE_v1", ["hba1c_n", "age_n", "sex_m", "bmi_n"])
        rec = {
            "dropped_gene": g,
            "n_genes": cov_loo["n_present"],
            **_focus_beta(fit_loo, "hba1c_n"),
        }
        loo_rows.append(rec)
        rows.append(
            {
                "analysis": "leave_one_gene_out",
                "series": "GSE50244",
                "set": f"IDENTITY_CORE_v1_minus_{g}",
                "n_genes": cov_loo["n_present"],
                "missing": g,
                **_focus_beta(fit_loo, "hba1c_n"),
            }
        )
    pd.DataFrame(loo_rows).to_csv(SENS_OUT / "leave_one_gene_out_GSE50244.csv", index=False)

    # 4. GSE76894 identity without batch 6 (batch 6 is entirely T2D).
    d768_nob6 = d768.loc[d768["batch"].astype(str) != "6"].copy()
    # rebuild batch dummies on the reduced table
    batch_dummies = pd.get_dummies(d768_nob6["batch"].astype(str), prefix="batch", drop_first=True).astype(float)
    keep_cols = [c for c in d768_nob6.columns if not str(c).startswith("batch_")]
    d768_nob6 = pd.concat([d768_nob6[keep_cols], batch_dummies], axis=1)
    fit_nob6 = ols_fit(d768_nob6, "IDENTITY_CORE_v1", ["t2d", "age_n", "sex_m", "bmi_n"] + list(batch_dummies.columns))
    rows.append(
        {
            "analysis": "drop_batch6",
            "series": "GSE76894",
            "set": "IDENTITY_CORE_v1",
            "n_genes": 16,
            "missing": "batch 6 samples excluded",
            **_focus_beta(fit_nob6, "t2d"),
        }
    )
    fit_nob6_stress = ols_fit(d768_nob6, "STRESS_CORE_v1", ["t2d", "age_n", "sex_m", "bmi_n"] + list(batch_dummies.columns))
    rows.append(
        {
            "analysis": "drop_batch6",
            "series": "GSE76894",
            "set": "STRESS_CORE_v1",
            "n_genes": 18,
            "missing": "batch 6 samples excluded",
            **_focus_beta(fit_nob6_stress, "t2d"),
        }
    )

    # 5. Outcome permutation of the four primary models (phenotype labels shuffled).
    rng = np.random.default_rng(RNG_SEED)
    primary = {
        "GSE50244_IDENTITY": ols_fit(d502, "IDENTITY_CORE_v1", ["hba1c_n", "age_n", "sex_m", "bmi_n"]).get("beta_hba1c_n"),
        "GSE50244_STRESS": ols_fit(d502, "STRESS_CORE_v1", ["hba1c_n", "age_n", "sex_m", "bmi_n"]).get("beta_hba1c_n"),
        "GSE76894_IDENTITY": ols_fit(d768, "IDENTITY_CORE_v1", ["t2d", "age_n", "sex_m", "bmi_n"] + batch_cols).get("beta_t2d"),
        "GSE76894_STRESS": ols_fit(d768, "STRESS_CORE_v1", ["t2d", "age_n", "sex_m", "bmi_n"] + batch_cols).get("beta_t2d"),
    }
    nulls = {k: [] for k in primary}
    work502 = d502.copy()
    work768 = d768.copy()
    for _ in range(N_PERM):
        work502["hba1c_n"] = rng.permutation(work502["hba1c_n"].to_numpy())
        fit_i = ols_fit(work502, "IDENTITY_CORE_v1", ["hba1c_n", "age_n", "sex_m", "bmi_n"])
        fit_s = ols_fit(work502, "STRESS_CORE_v1", ["hba1c_n", "age_n", "sex_m", "bmi_n"])
        nulls["GSE50244_IDENTITY"].append(fit_i.get("beta_hba1c_n", np.nan))
        nulls["GSE50244_STRESS"].append(fit_s.get("beta_hba1c_n", np.nan))
        work768["t2d"] = rng.permutation(work768["t2d"].to_numpy())
        fit_vi = ols_fit(work768, "IDENTITY_CORE_v1", ["t2d", "age_n", "sex_m", "bmi_n"] + batch_cols)
        fit_vs = ols_fit(work768, "STRESS_CORE_v1", ["t2d", "age_n", "sex_m", "bmi_n"] + batch_cols)
        nulls["GSE76894_IDENTITY"].append(fit_vi.get("beta_t2d", np.nan))
        nulls["GSE76894_STRESS"].append(fit_vs.get("beta_t2d", np.nan))
        # restore original outcomes for next independent shuffle
        work502["hba1c_n"] = d502["hba1c_n"].to_numpy()
        work768["t2d"] = d768["t2d"].to_numpy()

    perm_rows = []
    for key, obs in primary.items():
        arr = np.array(nulls[key], dtype=float)
        rec = {
            "analysis": "outcome_permutation",
            "test": key,
            "n_perm": N_PERM,
            "observed_beta": obs,
            "null_mean": float(np.nanmean(arr)),
            "null_sd": float(np.nanstd(arr)),
            "permutation_p": permutation_p(obs, arr),
        }
        perm_rows.append(rec)
        rows.append(
            {
                "analysis": "outcome_permutation",
                "series": key.split("_")[0],
                "set": "IDENTITY_CORE_v1" if "IDENTITY" in key else "STRESS_CORE_v1",
                "n_genes": None,
                "missing": f"{N_PERM} permutations",
                "n": None,
                "beta": obs,
                "se": rec["null_sd"],
                "p": rec["permutation_p"],
                "ci95_low": None,
                "ci95_high": None,
                "r2": None,
                "error": None,
            }
        )
    pd.DataFrame(perm_rows).to_csv(SENS_OUT / "outcome_permutation.csv", index=False)
    pd.DataFrame(rows).to_csv(SENS_OUT / "sensitivity_summary.csv", index=False)
    (SENS_OUT / "sensitivity_notes.json").write_text(
        json.dumps(
            {
                "n_perm": N_PERM,
                "rng_seed": RNG_SEED,
                "drop5_genes": IDENTITY_MISSING_ON_GSE76894,
                "note": "Sensitivity only. Primary four tests remain the locked family.",
            },
            indent=2,
        )
    )
    log(pd.DataFrame(rows).to_string(index=False))
    log(f"Wrote {SENS_OUT}")


if __name__ == "__main__":
    main()
