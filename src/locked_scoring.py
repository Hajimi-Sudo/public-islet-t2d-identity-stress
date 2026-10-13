#!/usr/bin/env python3
"""Locked identity/stress scores. Frozen gene sets; no DE-driven gene adding."""
from __future__ import annotations

import gzip
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

from paths import (
    ALL_SETS,
    ANNOT,
    GEO,
    GSE50244_ID_MAP,
    PHENO,
    PRIMARY_SETS,
    SETS,
)

EXPECTED_ANNOT_BYTES = {
    "GPL570.annot.gz": 8471521,
    "GPL6244.annot.gz": 7193437,
    "GPL1352.annot.gz": 8765182,
}

SYM_RE = re.compile(r"^[A-Za-z][A-Za-z0-9.\-]*$")


def log(msg: str) -> None:
    print(msg, flush=True)


def parse_gpl_annot(path: Path) -> pd.DataFrame:
    rows = []
    header = None
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("#") or line.startswith("^") or line.startswith("!"):
                continue
            if header is None:
                header = line.rstrip("\n").split("\t")
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            rec = {header[i]: parts[i] if i < len(parts) else "" for i in range(len(header))}
            rows.append(rec)
    df = pd.DataFrame(rows)
    id_col = None
    for c in df.columns:
        if c.lower() in {"id", "id_ref"}:
            id_col = c
            break
    if id_col is None:
        id_col = df.columns[0]
    gene_col = None
    for cand in ["Gene symbol", "Gene Symbol", "gene_assignment", "Gene title"]:
        if cand in df.columns:
            gene_col = cand
            break
    if gene_col is None:
        raise ValueError(f"No gene symbol column in {path}: {list(df.columns)[:12]}")
    out = df[[id_col, gene_col]].rename(columns={id_col: "probe", gene_col: "gene_raw"})
    out["probe"] = out["probe"].astype(str).str.strip().str.strip('"')
    return out


def explode_symbols(raw: str) -> list[str]:
    if raw is None or (isinstance(raw, float) and math.isnan(raw)):
        return []
    s = str(raw).strip().strip('"')
    if s in {"", "---", "NA", "null", "NONE"}:
        return []
    symbols: list[str] = []
    if "///" in s:
        for tok in s.split("///"):
            tok = tok.strip()
            if SYM_RE.fullmatch(tok):
                symbols.append(tok)
            elif "//" in tok:
                parts = [t.strip() for t in tok.split("//")]
                if len(parts) >= 2 and SYM_RE.fullmatch(parts[1]):
                    symbols.append(parts[1])
        return list(dict.fromkeys(symbols))
    if "//" in s:
        toks = [t.strip() for t in s.split("//")]
        if len(toks) >= 2:
            for i in range(1, len(toks), 3):
                sym = toks[i].strip()
                if SYM_RE.fullmatch(sym):
                    symbols.append(sym)
        return list(dict.fromkeys(symbols))
    parts = re.split(r"[;,/| ]+", s)
    return [p for p in parts if SYM_RE.fullmatch(p)]


def probe_to_genes(annot_df: pd.DataFrame) -> dict[str, list[str]]:
    m = {}
    for probe, raw in zip(annot_df["probe"], annot_df["gene_raw"]):
        m[probe] = explode_symbols(raw)
    return m


def annot_complete(path: Path) -> bool:
    if not path.exists() or path.stat().st_size < 100_000:
        return False
    expected = EXPECTED_ANNOT_BYTES.get(path.name)
    if expected is not None and path.stat().st_size != expected:
        return False
    try:
        with gzip.open(path, "rb") as f:
            f.read(64)
        return True
    except OSError:
        return False


def load_probe_map(platform: str) -> tuple[dict[str, list[str]] | None, str]:
    annot = ANNOT / f"{platform}.annot.gz"
    if annot_complete(annot):
        return probe_to_genes(parse_gpl_annot(annot)), f"geo_annot:{annot.name}"
    return None, "missing"


def read_series_matrix_expr(path: Path) -> pd.DataFrame:
    header = None
    rows = []
    with gzip.open(path, "rt", errors="replace") as f:
        in_table = False
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                header_line = next(f)
                header = [x.strip().strip('"') for x in header_line.rstrip("\n").split("\t")]
                continue
            if not in_table:
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            parts = line.rstrip("\n").split("\t")
            probe = parts[0].strip().strip('"')
            vals = []
            for x in parts[1:]:
                x = x.strip().strip('"')
                if x in {"", "NA", "NaN", "null"}:
                    vals.append(np.nan)
                else:
                    vals.append(float(x))
            rows.append((probe, vals))
    cols = header[1:]
    data = np.array([r[1] for r in rows], dtype=float)
    idx = [r[0] for r in rows]
    return pd.DataFrame(data, index=idx, columns=cols)


def read_gse50244_gene_matrix(path: Path, pheno: pd.DataFrame) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", compression="gzip")
    idcol = df.columns[0]
    df = df.set_index(idcol)
    title_to_gsm = dict(zip(pheno["title"].astype(str), pheno["geo_accession"]))
    newcols = []
    for c in df.columns:
        gsm = title_to_gsm.get(str(c))
        if gsm is None:
            raise KeyError(f"GSE50244 matrix column {c} not in phenotype titles")
        newcols.append(gsm)
    df.columns = newcols
    return df.astype(float)


def collapse_probes_to_genes(
    expr: pd.DataFrame, probe_map: dict[str, list[str]], gene_sets: dict[str, list[str]]
) -> pd.DataFrame:
    needed = set()
    for genes in gene_sets.values():
        needed.update(genes)
    gene_to_probes: dict[str, list[str]] = {g: [] for g in needed}
    for probe, symbols in probe_map.items():
        if probe not in expr.index:
            continue
        for g in symbols:
            if g in gene_to_probes:
                gene_to_probes[g].append(probe)
    rows = {}
    for g, probes in gene_to_probes.items():
        if not probes:
            continue
        sub = expr.loc[probes]
        rows[g] = sub.median(axis=0)
    if not rows:
        return pd.DataFrame(index=[], columns=expr.columns)
    return pd.DataFrame(rows).T


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    arr = df.to_numpy(dtype=float)
    mu = np.nanmean(arr, axis=1, keepdims=True)
    sd = np.nanstd(arr, axis=1, keepdims=True, ddof=0)
    sd = np.where(sd == 0, np.nan, sd)
    z = (arr - mu) / sd
    return pd.DataFrame(z, index=df.index, columns=df.columns)


def score_sets(
    gene_expr: pd.DataFrame, gene_sets: dict[str, list[str]], id_map: dict[str, str] | None = None, reduce: str = "mean"
) -> tuple[pd.DataFrame, pd.DataFrame]:
    id_map = id_map or {}
    reduce = str(reduce).lower()
    present_idx = set(gene_expr.index)
    z = zscore_rows(gene_expr)
    score_rows = {}
    coverage = []
    for set_name, genes in gene_sets.items():
        mapped = []
        missing = []
        used_ids = []
        for g in genes:
            gid = id_map.get(g, g)
            if gid in present_idx:
                mapped.append(g)
                used_ids.append(gid)
            elif g in present_idx:
                mapped.append(g)
                used_ids.append(g)
            else:
                missing.append(g)
        if used_ids:
            block = z.loc[used_ids]
            if reduce == "mean":
                score_rows[set_name] = block.mean(axis=0, skipna=True)
            elif reduce == "median":
                score_rows[set_name] = block.median(axis=0, skipna=True)
            else:
                raise ValueError(f"unknown reduce={reduce!r}")
        else:
            score_rows[set_name] = pd.Series(np.nan, index=gene_expr.columns)
        coverage.append(
            {
                "set": set_name,
                "n_locked": len(genes),
                "n_present": len(mapped),
                "frac_present": len(mapped) / len(genes) if genes else np.nan,
                "present": ",".join(mapped),
                "missing": ",".join(missing),
            }
        )
    scores = pd.DataFrame(score_rows)
    scores.index.name = "geo_accession"
    return scores, pd.DataFrame(coverage)


def ols_fit(df: pd.DataFrame, y: str, x_terms: list[str]) -> dict:
    work = df[[y] + x_terms].dropna().copy()
    if len(work) < 8:
        return {"n": len(work), "error": "too_few"}
    Y = work[y].astype(float)
    X = work[x_terms].astype(float)
    X = sm.add_constant(X, has_constant="add")
    model = sm.OLS(Y, X).fit()
    out = {
        "n": int(len(work)),
        "r2": float(model.rsquared),
        "adj_r2": float(model.rsquared_adj),
    }
    for term in ["const"] + x_terms:
        if term not in model.params.index:
            continue
        out[f"beta_{term}"] = float(model.params[term])
        out[f"se_{term}"] = float(model.bse[term])
        out[f"t_{term}"] = float(model.tvalues[term])
        out[f"p_{term}"] = float(model.pvalues[term])
        ci = model.conf_int().loc[term]
        out[f"ci95_low_{term}"] = float(ci.iloc[0])
        out[f"ci95_high_{term}"] = float(ci.iloc[1])
    return out


def spearman(df: pd.DataFrame, a: str, b: str) -> dict:
    work = df[[a, b]].dropna()
    if len(work) < 8:
        return {"n": len(work), "error": "too_few"}
    rho, p = stats.spearmanr(work[a], work[b])
    return {"n": int(len(work)), "rho": float(rho), "p": float(p)}


def welch(df: pd.DataFrame, y: str, group: str, pos: str, neg: str) -> dict:
    a = df.loc[df[group] == pos, y].dropna()
    b = df.loc[df[group] == neg, y].dropna()
    if len(a) < 3 or len(b) < 3:
        return {"n_pos": int(len(a)), "n_neg": int(len(b)), "error": "too_few"}
    t, p = stats.ttest_ind(a, b, equal_var=False)
    return {
        "n_pos": int(len(a)),
        "n_neg": int(len(b)),
        "mean_pos": float(a.mean()),
        "mean_neg": float(b.mean()),
        "diff": float(a.mean() - b.mean()),
        "t": float(t),
        "p": float(p),
    }


def load_locked_sets() -> dict[str, list[str]]:
    locked = json.loads(SETS.read_text())
    return {k: v["genes"] for k, v in locked["sets"].items() if k in ALL_SETS}


def load_matrices(gene_sets: dict[str, list[str]]) -> dict:
    log("Loading phenotypes")
    p502 = pd.read_csv(PHENO / "GSE50244_phenotype.csv")
    p386 = pd.read_csv(PHENO / "GSE38642_phenotype.csv")
    p768 = pd.read_csv(PHENO / "GSE76894_phenotype.csv")
    p209 = pd.read_csv(PHENO / "GSE20966_phenotype.csv")

    log("Loading GPL annotations")
    gpl570, src570 = load_probe_map("GPL570")
    gpl6244, src6244 = load_probe_map("GPL6244")
    gpl1352, src1352 = load_probe_map("GPL1352")
    if gpl570 is None:
        raise FileNotFoundError("GPL570 probe map required for GSE76894. Run src/download_public_data.py")
    mapping_provenance = {
        "GPL570": {"source": src570, "n_probes": len(gpl570), "n_with_symbol": sum(1 for v in gpl570.values() if v)},
        "GPL6244": {
            "source": src6244,
            "n_probes": 0 if gpl6244 is None else len(gpl6244),
            "n_with_symbol": 0 if gpl6244 is None else sum(1 for v in gpl6244.values() if v),
        },
        "GPL1352": {
            "source": src1352,
            "n_probes": 0 if gpl1352 is None else len(gpl1352),
            "n_with_symbol": 0 if gpl1352 is None else sum(1 for v in gpl1352.values() if v),
        },
    }

    log("Loading expression")
    e502 = read_gse50244_gene_matrix(
        GEO / "GSE50244/GSE50244_Genes_counts_TMM_NormLength_atLeastMAF5_expressed.txt.gz", p502
    )
    e768_probe = read_series_matrix_expr(GEO / "GSE76894/GSE76894_series_matrix.txt.gz")
    e386_probe = None
    e209_probe = None
    if gpl6244 is not None:
        e386_probe = read_series_matrix_expr(GEO / "GSE38642/GSE38642_series_matrix.txt.gz")
    else:
        log("SKIP GSE38642: GPL6244 map incomplete")
    if gpl1352 is not None and src1352.startswith("geo_annot"):
        e209_probe = read_series_matrix_expr(GEO / "GSE20966/GSE20966_series_matrix.txt.gz")
    else:
        log(f"SKIP GSE20966: GPL1352 official annot not ready (source={src1352})")

    log("Collapsing probes")
    e768 = collapse_probes_to_genes(e768_probe, gpl570, gene_sets)
    e386 = None
    e209 = None
    if e386_probe is not None and gpl6244 is not None:
        e386 = collapse_probes_to_genes(e386_probe, gpl6244, gene_sets)
    if e209_probe is not None and gpl1352 is not None:
        e209 = collapse_probes_to_genes(e209_probe, gpl1352, gene_sets)

    return {
        "p502": p502,
        "p386": p386,
        "p768": p768,
        "p209": p209,
        "e502": e502,
        "e768": e768,
        "e386": e386,
        "e209": e209,
        "mapping_provenance": mapping_provenance,
    }


def attach(pheno: pd.DataFrame, scores: pd.DataFrame, gsm_col: str = "geo_accession") -> pd.DataFrame:
    return pheno.copy().merge(
        scores.reset_index(), left_on=gsm_col, right_on="geo_accession", how="left", suffixes=("", "_scoreidx")
    )


def prepare_discovery(p502: pd.DataFrame, s502: pd.DataFrame) -> pd.DataFrame:
    d502 = attach(p502, s502)
    d502["hba1c_n"] = pd.to_numeric(d502["hba1c"], errors="coerce")
    d502["age_n"] = pd.to_numeric(d502["age"], errors="coerce")
    d502["bmi_n"] = pd.to_numeric(d502["bmi"], errors="coerce")
    d502["sex_m"] = (d502["gender"].astype(str).str.lower().str[0] == "m").astype(float)
    d502["hba1c65"] = np.where(d502["hba1c_n"].isna(), np.nan, (d502["hba1c_n"] >= 6.5).astype(float))
    return d502


def prepare_validation(p768: pd.DataFrame, s768: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    d768 = attach(p768, s768)
    d768["age_n"] = pd.to_numeric(d768["age (years)"], errors="coerce")
    d768["bmi_n"] = pd.to_numeric(d768["body mass index (bmi)"], errors="coerce")
    d768["sex_m"] = (d768["gender (m/f)"].astype(str).str.lower().str[0] == "m").astype(float)
    d768["t2d"] = (d768["t2d_label"] == "T2D").astype(float)
    batch_dummies = pd.get_dummies(d768["batch"].astype(str), prefix="batch", drop_first=True).astype(float)
    d768 = pd.concat([d768, batch_dummies], axis=1)
    return d768, list(batch_dummies.columns)


def prepare_companion(p386: pd.DataFrame, s386: pd.DataFrame) -> pd.DataFrame:
    d386 = attach(p386, s386)
    d386["hba1c_n"] = pd.to_numeric(d386["hba1c"].replace("--", np.nan), errors="coerce")
    d386["age_n"] = pd.to_numeric(d386["age"], errors="coerce")
    d386["bmi_n"] = pd.to_numeric(d386["bmi"], errors="coerce")
    d386["sex_m"] = (d386["sex"].astype(str).str.lower().str[0] == "m").astype(float)
    d386["t2d"] = (d386["t2d_label"] == "T2D").astype(float)
    return d386


def prepare_orthogonal(p209: pd.DataFrame, s209: pd.DataFrame) -> pd.DataFrame:
    d209 = attach(p209, s209)
    d209["age_n"] = pd.to_numeric(d209["age (years)"], errors="coerce")
    d209["bmi_n"] = pd.to_numeric(d209["bmi (body mass index, kg/m2)"], errors="coerce")
    d209["sex_m"] = (d209["gender"].astype(str).str.lower().str[0] == "m").astype(float)
    d209["t2d"] = (d209["t2d_label"] == "T2D").astype(float)
    return d209


def family_of(set_name: str) -> str:
    return "primary" if set_name in PRIMARY_SETS else "secondary"
