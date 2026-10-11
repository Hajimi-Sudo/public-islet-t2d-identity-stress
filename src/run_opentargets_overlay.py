#!/usr/bin/env python3
"""Open Targets MONDO_0005148 overlay. Not GWAS Catalog. Not Mendelian randomization."""
from __future__ import annotations

import json

import pandas as pd

from locked_scoring import load_locked_sets, log
from paths import OT_OUT, OVERLAYS


def datatype_score(row: dict, key: str) -> float | None:
    for item in row.get("datatypeScores", []) or []:
        if item.get("id") == key:
            return float(item.get("score"))
    return None


def main() -> None:
    OT_OUT.mkdir(parents=True, exist_ok=True)
    snap = json.loads((OVERLAYS / "opentargets_mondo0005148_top500.json").read_text())
    rows = snap["data"]["disease"]["associatedTargets"]["rows"]
    by_symbol = {}
    for rec in rows:
        symbol = rec["target"]["approvedSymbol"]
        by_symbol[symbol] = {
            "overall_score": float(rec.get("score")),
            "genetic_association_score": datatype_score(rec, "genetic_association"),
            "known_drug_score": datatype_score(rec, "known_drug"),
        }

    gene_sets = load_locked_sets()
    gene_rows = []
    summaries = []
    for set_name, genes in gene_sets.items():
        in_top = []
        with_ga = []
        for g in genes:
            hit = by_symbol.get(g)
            gene_rows.append(
                {
                    "set": set_name,
                    "gene": g,
                    "in_opentargets_t2d_top500": hit is not None,
                    "overall_score": None if hit is None else hit["overall_score"],
                    "genetic_association_score": None if hit is None else hit["genetic_association_score"],
                    "known_drug_score": None if hit is None else hit["known_drug_score"],
                }
            )
            if hit is not None:
                in_top.append(g)
                if hit["genetic_association_score"] not in (None, 0.0):
                    with_ga.append(g)
        summaries.append(
            {
                "set": set_name,
                "n_locked": len(genes),
                "n_in_ot_t2d_top500": len(in_top),
                "frac_in_ot_t2d_top500": len(in_top) / len(genes) if genes else None,
                "n_with_genetic_association_score": len(with_ga),
                "genes_in_top500": ",".join(in_top),
                "genes_with_genetic_association": ",".join(with_ga),
            }
        )

    pd.DataFrame(gene_rows).to_csv(OT_OUT / "opentargets_t2d_top500_overlap.csv", index=False)
    pd.DataFrame(summaries).to_csv(OT_OUT / "opentargets_t2d_top500_set_summary.csv", index=False)
    source = json.loads((OVERLAYS / "overlay_source.json").read_text())
    source["not_mr"] = True
    source["not_gwas_catalog_full_dump"] = True
    (OT_OUT / "overlay_source.json").write_text(json.dumps(source, indent=2))
    log(pd.DataFrame(summaries).to_string(index=False))
    log(f"Wrote {OT_OUT}")


if __name__ == "__main__":
    main()
