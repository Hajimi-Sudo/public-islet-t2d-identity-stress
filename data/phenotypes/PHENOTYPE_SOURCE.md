# Phenotype tables shipped with the code

These CSVs are **donor-level GEO sample characteristics**, parsed from each
series matrix / SOFT metadata. They are small enough to keep in git so that
the locked models can be rerun without re-parsing GEO headers.

| File | n | Label used in the analysis |
|---|---:|---|
| `GSE50244_phenotype.csv` | 89 | continuous `hba1c` (77 non-missing). `t2d_label` is empty / not a GEO field and must not be used as the primary endpoint |
| `GSE76894_phenotype.csv` | 103 | `t2d_label` from `diabetes status (nd …, t2d …)`; models also use `age (years)`, `gender (m/f)`, `body mass index (bmi)`, `batch` |
| `GSE38642_phenotype.csv` | 63 | `t2d_label` from `status`; same Nordic biobank as GSE50244 (53/89 donor overlap) |
| `GSE20966_phenotype.csv` | 20 | `t2d_label` from `disease` |

## Intentionally not shipped

- `GSE50244_phenotype_with_transferred_t2d.csv` — transferring GSE38642 clinical T2D onto GSE50244 is forbidden. Among overlapping donors, HbA1c ≥ 6.5 and GSE38642 “non-diabetic” are discordant for at least three donors.

## Rules encoded in the scripts

1. GSE50244 discovery = `IDENTITY` / `STRESS` ~ HbA1c + age + sex + BMI.
2. GSE76894 validation = same scores ~ T2D + age + sex + BMI + batch.
3. GSE38642 is a companion, not an independent cohort.
4. GSE20966 is an underpowered orthogonal check.
5. ICU glycemia on GSE76894 is **not** used as a diagnostic endpoint.
