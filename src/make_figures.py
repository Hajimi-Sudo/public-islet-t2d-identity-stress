#!/usr/bin/env python3
"""Publication figures from locked-score tables. No gene-set peeking."""
from __future__ import annotations

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

from paths import FIG_OUT, GTEX_OUT, LOCKED_OUT, TAB_OUT

NAVY = "#1B4F72"
RUST = "#A04000"
SLATE = "#2C3E50"
GRAY = "#7F8C8D"
ND = "#5D6D7E"
T2D = "#922B21"
PANEL = "#F7F9F9"

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def save(fig: plt.Figure, stem: str) -> None:
    FIG_OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT / f"{stem}.png", bbox_inches="tight", facecolor="white")
    fig.savefig(FIG_OUT / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def fig1_design() -> None:
    fig, ax = plt.subplots(figsize=(9.4, 5.6))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6.2)
    ax.axis("off")
    ax.set_title("Study design: freeze both programs, then test the same scores", loc="left", fontweight="bold", color=SLATE)

    def box(x, y, w, h, text, fc="#EAF2F8", ec=NAVY, fs=8):
        r = mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04,rounding_size=0.08", linewidth=1.1, facecolor=fc, edgecolor=ec)
        ax.add_patch(r)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=SLATE)

    box(0.25, 4.85, 4.6, 1.1, "IDENTITY_CORE_v1  (21 genes)\nINS, PDX1, MAFA, NKX6-1, SLC30A8…\nLocked before any T2D contrast", fc="#D4E6F1", ec=NAVY)
    box(5.15, 4.85, 4.6, 1.1, "STRESS_CORE_v1  (18 genes)\nDDIT3, ATF4, XBP1, HSPA5, TXNIP…\nPrespecified direction: up in T2D / hyperglycemia", fc="#F5CBA7", ec=RUST)
    ax.text(5.0, 4.62, "Both scores applied to every cohort below", ha="center", va="center", fontsize=8, color=GRAY, style="italic")
    box(0.25, 2.95, 2.25, 1.4, "Discovery\nGSE50244 RNA-seq\nn=89 islets\nHbA1c endpoint\n(n=77 complete)", fc="#E8F6F3")
    box(2.65, 2.95, 2.25, 1.4, "External\nvalidation\nGSE76894 array\n19 T2D / 84 ND\nage/sex/BMI/batch", fc="#D5F5E3", ec="#196F3D")
    box(5.05, 2.95, 2.25, 1.4, "Same-biobank\ncompanion\nGSE38642 array\n53/89 donor overlap\nnot independent", fc="#F4F6F7", ec=GRAY)
    box(7.45, 2.95, 2.3, 1.4, "Orthogonal LCM\nGSE20966\n10 vs 10\nunderpowered check", fc="#F4F6F7", ec=GRAY)
    box(0.25, 0.25, 4.6, 2.35, "Primary family (4 tests)\nIdentity and stress on:\n• GSE50244  score ~ HbA1c + age + sex + BMI\n• GSE76894  score ~ T2D + age + sex + BMI + batch\nNo gene added from differential expression.\nStress-up is allowed to fail.", fc="#FCF3CF", ec="#B7950B")
    box(5.15, 0.25, 4.6, 2.35, "Overlays, not extra endpoints\nGTEx v8 public median TPM:\ntissue selectivity (bulk pancreas, not islets)\nOpen Targets T2D overlay (not MR)\nNHANES / MIMIC / WGCNA / two-sample MR excluded", fc="#EBF5FB", ec=NAVY)
    save(fig, "Figure1_study_design")


def fig2_primary_panels() -> None:
    d502 = pd.read_csv(LOCKED_OUT / "GSE50244_scores.csv")
    d768 = pd.read_csv(LOCKED_OUT / "GSE76894_scores.csv")
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 7.2))
    fig.subplots_adjust(wspace=0.32, hspace=0.38)

    def scatter(ax, ycol, color, ylabel, subtitle, beta, p):
        work = d502[["hba1c", ycol]].dropna()
        ax.scatter(work["hba1c"], work[ycol], s=22, c=color, alpha=0.75, edgecolors="white", linewidths=0.3, zorder=3)
        x = work["hba1c"].to_numpy()
        y = work[ycol].to_numpy()
        slope, intercept = np.polyfit(x, y, 1)
        xs = np.linspace(x.min(), x.max(), 50)
        ax.plot(xs, slope * xs + intercept, color=color, lw=1.6, zorder=2)
        ax.axvline(6.5, color=GRAY, ls="--", lw=0.8, alpha=0.8)
        ax.set_xlabel("Donor HbA1c (%)")
        ax.set_ylabel(ylabel)
        ax.set_title(subtitle, loc="left", fontweight="bold")
        ax.text(0.04, 0.04, f"unadjusted slope shown\nadjusted β={beta:.3f}, p={p:.3g}", transform=ax.transAxes, fontsize=7.5, color=SLATE, va="bottom")
        ax.set_facecolor(PANEL)

    scatter(axes[0, 0], "IDENTITY_CORE_v1", NAVY, "Identity score", "A  GSE50244 identity vs HbA1c", -0.200, 0.035)
    scatter(axes[0, 1], "STRESS_CORE_v1", RUST, "Stress/UPR score", "B  GSE50244 stress vs HbA1c", 0.031, 0.71)

    def boxes(ax, ycol, ylabel, subtitle, note):
        nd = d768.loc[d768["t2d_label"] == "ND", ycol].dropna()
        t2 = d768.loc[d768["t2d_label"] == "T2D", ycol].dropna()
        bp = ax.boxplot([nd, t2], positions=[1, 2], widths=0.55, patch_artist=True, showfliers=False, medianprops=dict(color="black", lw=1.2), whiskerprops=dict(color=SLATE, lw=1), capprops=dict(color=SLATE, lw=1))
        bp["boxes"][0].set_facecolor("#D5D8DC")
        bp["boxes"][1].set_facecolor("#F5B7B1")
        rng = np.random.default_rng(66)
        ax.scatter(1 + rng.uniform(-0.08, 0.08, len(nd)), nd, s=12, c=ND, alpha=0.55, zorder=3)
        ax.scatter(2 + rng.uniform(-0.08, 0.08, len(t2)), t2, s=12, c=T2D, alpha=0.7, zorder=3)
        ax.set_xticks([1, 2], [f"ND\nn={len(nd)}", f"T2D\nn={len(t2)}"])
        ax.set_ylabel(ylabel)
        ax.set_title(subtitle, loc="left", fontweight="bold")
        ax.text(0.04, 0.04, note, transform=ax.transAxes, fontsize=7.5, color=SLATE, va="bottom")
        ax.set_facecolor(PANEL)

    boxes(axes[1, 0], "IDENTITY_CORE_v1", "Identity score", "C  GSE76894 identity", "adjusted β_T2D=−0.569, p=0.015")
    boxes(axes[1, 1], "STRESS_CORE_v1", "Stress/UPR score", "D  GSE76894 stress", "adjusted β_T2D=−0.181, p=0.19\n(not up)")
    save(fig, "Figure2_discovery_validation")


def fig3_forest() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.6), sharey=False)
    fig.subplots_adjust(wspace=0.45, left=0.18)

    def panel(ax, rows, xlabel, title, xlim):
        ax.axvline(0, color=GRAY, lw=0.9)
        for i, (lab, b, lo, hi, p, c) in enumerate(rows[::-1]):
            ax.plot([lo, hi], [i, i], color=c, lw=2.4, solid_capstyle="round")
            ax.plot(b, i, "o", color=c, ms=8, zorder=3)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([f"{rows[1][0]}\np={rows[1][4]:.3g}", f"{rows[0][0]}\np={rows[0][4]:.3g}"])
        ax.set_xlabel(xlabel)
        ax.set_title(title, loc="left", fontweight="bold")
        ax.set_xlim(*xlim)
        ax.set_ylim(-0.7, 1.7)
        ax.set_facecolor(PANEL)

    left = [("Identity", -0.199851, -0.385076, -0.014625, 0.0348, NAVY), ("Stress/UPR", 0.030901, -0.135330, 0.197133, 0.712, RUST)]
    right = [("Identity", -0.568676, -1.023723, -0.113628, 0.0149, NAVY), ("Stress/UPR", -0.181425, -0.454398, 0.091547, 0.190, RUST)]
    panel(axes[0], left, "β per HbA1c point (age/sex/BMI)", "A  GSE50244 discovery, n=77", (-0.50, 0.28))
    panel(axes[1], right, "β for recorded T2D (age/sex/BMI/batch)", "B  GSE76894 validation, n=102", (-1.20, 0.28))
    save(fig, "Figure3_primary_forest")


def fig4_gtex() -> None:
    g = pd.read_csv(GTEX_OUT / "gtex_set_mean_tpm_focus_tissues.csv")
    tissues = ["Pancreas", "Liver", "Whole Blood", "Muscle - Skeletal", "Adipose - Subcutaneous", "Thyroid", "Adrenal Gland", "Brain - Cortex"]
    short = ["Pancreas", "Liver", "Blood", "Muscle", "Adipose SC", "Thyroid", "Adrenal", "Cortex"]
    iden = g.loc[g["set"] == "IDENTITY_CORE_v1", tissues].iloc[0].to_numpy(dtype=float)
    stre = g.loc[g["set"] == "STRESS_CORE_v1", tissues].iloc[0].to_numpy(dtype=float)
    x = np.arange(len(tissues))
    w = 0.38
    fig, ax = plt.subplots(figsize=(8.8, 4.2))
    ax.bar(x - w / 2, np.log10(iden + 0.1), width=w, color=NAVY, label="Identity core (mean TPM)")
    ax.bar(x + w / 2, np.log10(stre + 0.1), width=w, color=RUST, label="Stress core (mean TPM)")
    ax.set_xticks(x, short, rotation=25, ha="right")
    ax.set_ylabel("log10(mean TPM + 0.1)")
    ax.set_title("GTEx v8 public medians: identity high in bulk pancreas; stress genes show limited tissue specificity", loc="left", fontweight="bold")
    ax.legend(frameon=False, loc="upper left")
    ax.set_facecolor(PANEL)
    ax.text(0.0, -0.32, "Bulk pancreas, not purified islets. INS dominates the identity mean; selectivity remains after removing INS.", transform=ax.transAxes, fontsize=7.5, color=GRAY)
    save(fig, "Figure4_gtex_overlay")


def write_tables() -> None:
    TAB_OUT.mkdir(parents=True, exist_ok=True)
    t1 = pd.DataFrame(
        [
            {"Series": "GSE50244", "Role": "Discovery", "Platform": "RNA-seq (GPL11154)", "n": 89, "T2D / ND": "no GEO T2D field; HbA1c on 77/89 (11 ≥ 6.5%)", "Age, years": "58.3 ± 10.7", "BMI": "26.0 ± 3.1", "Sex (M/F)": "54/35"},
            {"Series": "GSE76894", "Role": "External validation", "Platform": "HG-U133 Plus 2 (GPL570)", "n": 103, "T2D / ND": "19 / 84", "Age, years": "ND 60.0 ± 16.1; T2D 72.2 ± 7.5", "BMI": "ND 25.8 ± 4.2; T2D 26.5 ± 3.6", "Sex (M/F)": "ND 38/46; T2D 13/6"},
            {"Series": "GSE38642", "Role": "Same-biobank companion", "Platform": "HuGene 1.0 ST (GPL6244)", "n": 63, "T2D / ND": "9 / 54", "Age, years": "ND 56.6 ± 9.8; T2D 57.0 ± 13.1", "BMI": "ND 26.0; T2D 28.5", "Sex (M/F)": "ND 31/23; T2D 5/4"},
            {"Series": "GSE20966", "Role": "Orthogonal LCM check", "Platform": "X3P (GPL1352)", "n": 20, "T2D / ND": "10 / 10", "Age, years": "ND 60.3; T2D 67.3", "BMI": "ND 30.6; T2D 30.9", "Sex (M/F)": "ND 6/4; T2D 7/3"},
        ]
    )
    t1.to_csv(TAB_OUT / "Table1_cohorts.csv", index=False)
    t2 = pd.DataFrame(
        [
            {"Series": "GSE50244", "Set": "IDENTITY_CORE_v1", "Model": "score ~ HbA1c + age + sex + BMI", "n": 77, "beta": -0.200, "SE": 0.093, "CI95": "−0.385 to −0.015", "p": 0.035, "Protocol": "identity down", "Verdict": "supported"},
            {"Series": "GSE50244", "Set": "STRESS_CORE_v1", "Model": "score ~ HbA1c + age + sex + BMI", "n": 77, "beta": 0.031, "SE": 0.083, "CI95": "−0.135 to 0.197", "p": 0.71, "Protocol": "stress up", "Verdict": "not supported"},
            {"Series": "GSE76894", "Set": "IDENTITY_CORE_v1", "Model": "score ~ T2D + age + sex + BMI + batch", "n": 102, "beta": -0.569, "SE": 0.229, "CI95": "−1.024 to −0.114", "p": 0.015, "Protocol": "identity down", "Verdict": "supported"},
            {"Series": "GSE76894", "Set": "STRESS_CORE_v1", "Model": "score ~ T2D + age + sex + BMI + batch", "n": 102, "beta": -0.181, "SE": 0.137, "CI95": "−0.454 to 0.092", "p": 0.19, "Protocol": "stress up", "Verdict": "not supported"},
        ]
    )
    t2.to_csv(TAB_OUT / "Table2_primary_four_tests.csv", index=False)
    pd.read_csv(LOCKED_OUT / "gene_coverage.csv").to_csv(TAB_OUT / "Table3_gene_coverage.csv", index=False)


def main() -> None:
    fig1_design()
    fig2_primary_panels()
    fig3_forest()
    fig4_gtex()
    write_tables()
    print("wrote figures to", FIG_OUT)
    print("wrote tables to", TAB_OUT)


if __name__ == "__main__":
    main()
