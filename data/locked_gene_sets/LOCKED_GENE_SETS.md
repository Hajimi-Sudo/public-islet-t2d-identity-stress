# Locked gene sets v1

Frozen: 2026-10-07  
Rule: frozen **before** any T2D vs ND expression contrast on GSE50244 / GSE38642 / GSE76894 / GSE20966.

Machine copy: `locked_gene_sets_v1.json`

## IDENTITY_CORE_v1 (n=21)

INS, IAPP, PDX1, MAFA, MAFB, NKX6-1, NKX2-2, NEUROD1, PAX6, ISL1, SLC2A2, SLC30A8, GCK, PCSK1, PCSK2, ABCC8, KCNJ11, UCN3, FFAR1, GLP1R, GJD2

GSE50244 identifier exception: `INS` is stored as `IGF2+INS-IGF2+INS` in the processed Fadista gene matrix. Map to that fused ID; do not drop insulin from the score.

Missing in GSE50244 processed matrix after mapping: none expected except fused INS. `MAFA` is present as `MAFA`.

## STRESS_CORE_v1 (n=18)

DDIT3, ATF3, ATF4, ATF6, XBP1, HSPA5, EIF2AK3, ERN1, TRIB3, DNAJC3, PDIA4, TXNIP, PPP1R15A, HERPUD1, DERL3, WFS1, HSP90B1, HYOU1

All 15 original UPR core genes plus WFS1/HSP90B1/HYOU1 are literature UPR/ER-stress genes; HSP90B1 is absent from the Fadista processed matrix and will be skipped only in that matrix (documented missingness, not post-hoc gene dropping from a T2D result).

## HALLMARK_PANCREAS_BETA_CELLS (n=40, MSigDB, PMID 26771021)

Secondary identity program. Includes INS, PDX1, NKX6-1, IAPP, GCG, SST (whole-islet endocrine, not beta-only). Use as a prespecified sensitivity set, not as a beta-cell-pure score.

## HALLMARK_UNFOLDED_PROTEIN_RESPONSE (n=113, MSigDB, PMID 26771021)

Secondary stress program.

## Explicitly not locked from these GEO series

No gene is added or removed because it was DE in T2D in these four series. WGCNA hub genes are forbidden.
