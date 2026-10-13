#!/usr/bin/env python3
"""Detectable-effect checks around the locked four tests.

Not observed-effect post-hoc power. Reports two-sided alpha=0.05, 80% power
minimum detectable effects from fitted OLS coefficient SEs.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from scipy.optimize import brentq
from statsmodels.stats.power import TTestIndPower

from locked_scoring import load_locked_sets, load_matrices, log, prepare_discovery, prepare_validation, score_sets
from paths import GSE50244_ID_MAP, RESULTS

OUT = RESULTS / "power_analysis_v1"


def ols_full(df, y, x_terms):
    work = df[[y] + x_terms].dropna().copy()
    Y = work[y].astype(float)
    X = sm.add_constant(work[x_terms].astype(float), has_constant="add")
    return work, sm.OLS(Y, X).fit()


def ncp_for_power(df: float, alpha: float = 0.05, power: float = 0.80) -> float:
    tcrit = stats.t.ppf(1 - alpha / 2, df)

    def f(ncp: float) -> float:
        return float(stats.nct.sf(tcrit, df, ncp) - power)

    return float(brentq(f, 0.05, 8.0))


def power_from_se(beta, se, df, alpha=0.05):
    tcrit = stats.t.ppf(1 - alpha / 2, df)
    ncp = abs(beta) / se
    return float(stats.nct.sf(tcrit, df, ncp) + stats.nct.cdf(-tcrit, df, ncp))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    gene_sets = load_locked_sets()
    data = load_matrices(gene_sets)
    s502, _ = score_sets(data["e502"], gene_sets, GSE50244_ID_MAP)
    s768, _ = score_sets(data["e768"], gene_sets)
    d502 = prepare_discovery(data["p502"], s502)
    d768, batch_cols = prepare_validation(data["p768"], s768)
    base502 = ["hba1c_n", "age_n", "sex_m", "bmi_n"]
    base768 = ["t2d", "age_n", "sex_m", "bmi_n"] + batch_cols
    rows = []
    ncp_cache = {}
    for series, df, terms, focus in [
        ("GSE50244", d502, base502, "hba1c_n"),
        ("GSE76894", d768, base768, "t2d"),
    ]:
        for y in ["IDENTITY_CORE_v1", "STRESS_CORE_v1"]:
            work, model = ols_full(df, y, terms)
            df_resid = float(model.df_resid)
            se = float(model.bse[focus])
            beta = float(model.params[focus])
            if df_resid not in ncp_cache:
                ncp_cache[df_resid] = ncp_for_power(df_resid)
            mde = ncp_cache[df_resid] * se
            id_beta = abs(float(ols_full(df, "IDENTITY_CORE_v1", terms)[1].params[focus]))
            rows.append(
                {
                    "analysis": "ols_mde_80pct",
                    "series": series,
                    "set": y,
                    "focus": focus,
                    "n": int(len(work)),
                    "df_resid": df_resid,
                    "beta": beta,
                    "se": se,
                    "p": float(model.pvalues[focus]),
                    "sigma_resid": float(np.sqrt(model.scale)),
                    "mde_80pct_alpha05": mde,
                    "power_at_identity_sized_effect": power_from_se(id_beta, se, df_resid),
                }
            )
    n_t2d = int((d768["t2d_label"] == "T2D").sum())
    n_nd = int((d768["t2d_label"] == "ND").sum())
    ratio = n_nd / n_t2d
    ind = TTestIndPower()
    for d in [0.2, 0.5, 0.8]:
        rows.append(
            {
                "analysis": "unadjusted_two_sample_power",
                "series": "GSE76894",
                "n_t2d": n_t2d,
                "n_nd": n_nd,
                "cohens_d": d,
                "power": float(ind.power(effect_size=d, nobs1=n_t2d, alpha=0.05, ratio=ratio, alternative="two-sided")),
            }
        )
    d_mde = float(ind.solve_power(effect_size=None, nobs1=n_t2d, alpha=0.05, power=0.80, ratio=ratio, alternative="two-sided"))
    rows.append({"analysis": "unadjusted_two_sample_mde", "series": "GSE76894", "n_t2d": n_t2d, "n_nd": n_nd, "mde_cohens_d_80pct": d_mde})
    rdf = pd.DataFrame(rows)
    rdf.to_csv(OUT / "power_models.csv", index=False)
    compact = rdf.loc[rdf["analysis"] == "ols_mde_80pct"].copy()
    compact.to_csv(OUT / "power_compact.csv", index=False)
    (OUT / "power_notes.json").write_text(
        json.dumps(
            {
                "role": "detectable_effect_not_observed_power",
                "alpha": 0.05,
                "target_power": 0.80,
                "gse76894_t2d_n": n_t2d,
                "gse76894_nd_n": n_nd,
                "unadjusted_mde_cohens_d_80pct": d_mde,
            },
            indent=2,
        )
    )
    log(compact.to_string(index=False))
    log(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
