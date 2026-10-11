#!/usr/bin/env python3
"""Compare results/locked_scores_v1 to expected_results/locked_scores_v1.

Checks the four primary coefficients and GSE76894 identity coverage (16/21).
"""
from __future__ import annotations

import sys

import pandas as pd

from paths import EXPECTED, LOCKED_OUT

ATOL_BETA = 1e-6
ATOL_P = 1e-8


def _fail(msg: str) -> None:
    print(f"FAIL  {msg}", flush=True)
    raise SystemExit(1)


def _ok(msg: str) -> None:
    print(f"OK    {msg}", flush=True)


def _coef(row: pd.Series) -> tuple[str, float, float]:
    if row["series"] == "GSE50244":
        return "hba1c_n", float(row["stat_beta_hba1c_n"]), float(row["stat_p_hba1c_n"])
    return "t2d", float(row["stat_beta_t2d"]), float(row["stat_p_t2d"])


def main() -> None:
    got_path = LOCKED_OUT / "primary_four_tests.csv"
    exp_path = EXPECTED / "locked_scores_v1" / "primary_four_tests.csv"
    if not got_path.exists():
        _fail(f"missing {got_path}; run src/run_locked_scores.py first")
    got = pd.read_csv(got_path)
    exp = pd.read_csv(exp_path)
    keys = ["series", "set", "model"]
    keep = keys + ["stat_n", "stat_beta_hba1c_n", "stat_p_hba1c_n", "stat_beta_t2d", "stat_p_t2d"]
    merged = exp[keep].merge(got[keep], on=keys, suffixes=("_exp", "_got"), how="outer", indicator=True)
    if (merged["_merge"] != "both").any():
        _fail("primary four test rows do not match expected keys:\n" + merged.to_string(index=False))

    for _, row in merged.iterrows():
        term = "hba1c_n" if row["series"] == "GSE50244" else "t2d"
        b_exp = float(row[f"stat_beta_{term}_exp"])
        b_got = float(row[f"stat_beta_{term}_got"])
        p_exp = float(row[f"stat_p_{term}_exp"])
        p_got = float(row[f"stat_p_{term}_got"])
        n_exp = float(row["stat_n_exp"])
        n_got = float(row["stat_n_got"])
        label = f"{row['series']} {row['set']}"
        if abs(n_got - n_exp) > 0:
            _fail(f"{label} n {n_got} != {n_exp}")
        if abs(b_got - b_exp) > ATOL_BETA:
            _fail(f"{label} beta {b_got:.8f} != {b_exp:.8f}")
        if abs(p_got - p_exp) > ATOL_P:
            _fail(f"{label} p {p_got:.8g} != {p_exp:.8g}")
        _ok(f"{label}: β={b_got:.6f} p={p_got:.4g} n={int(n_got)}")

    cov_got = pd.read_csv(LOCKED_OUT / "gene_coverage.csv")
    cov_exp = pd.read_csv(EXPECTED / "locked_scores_v1" / "gene_coverage.csv")
    row = cov_got[(cov_got["series"] == "GSE76894") & (cov_got["set"] == "IDENTITY_CORE_v1")]
    if row.empty:
        _fail("missing GSE76894 IDENTITY coverage")
    n_present = int(row.iloc[0]["n_present"])
    if n_present != 16:
        _fail(f"GSE76894 identity coverage {n_present}/21 (expected 16/21)")
    missing = str(row.iloc[0]["missing"])
    for gene in ["PDX1", "MAFA", "KCNJ11", "UCN3", "GJD2"]:
        if gene not in missing:
            _fail(f"GSE76894 identity missing-list should include {gene}: {missing}")
    _ok(f"GSE76894 IDENTITY coverage 16/21 (missing {missing})")

    # Coverage n_present should match the frozen snapshot for primary sets.
    for series, set_name in [
        ("GSE50244", "IDENTITY_CORE_v1"),
        ("GSE50244", "STRESS_CORE_v1"),
        ("GSE76894", "IDENTITY_CORE_v1"),
        ("GSE76894", "STRESS_CORE_v1"),
    ]:
        a = cov_exp[(cov_exp["series"] == series) & (cov_exp["set"] == set_name)].iloc[0]
        b = cov_got[(cov_got["series"] == series) & (cov_got["set"] == set_name)].iloc[0]
        if int(a["n_present"]) != int(b["n_present"]):
            _fail(f"{series} {set_name} n_present {b['n_present']} != {a['n_present']}")

    print("Primary four tests match expected_results/locked_scores_v1.")
    return None


if __name__ == "__main__":
    main()
