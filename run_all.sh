#!/usr/bin/env bash
# Reproduce the locked analysis. Large GEO/GTEx files go to cache/ or $PUBLIC_ISLET_CACHE.
set -euo pipefail
cd "$(dirname "$0")"

export PYTHONPATH="${PWD}/src${PYTHONPATH:+:${PYTHONPATH}}"

echo "==> public matrices -> ${PUBLIC_ISLET_CACHE:-$PWD/cache}"
python3 src/download_public_data.py

echo "==> locked scores (primary four tests)"
python3 src/run_locked_scores.py

echo "==> GTEx overlay"
python3 src/run_gtex_overlay.py

echo "==> Open Targets overlay"
python3 src/run_opentargets_overlay.py

echo "==> sensitivity (does not replace the primary four tests)"
python3 src/run_sensitivity.py

echo "==> additional mean-versus-median scoring (does not replace the primary four tests)"
python3 src/run_median_scoring.py

echo "==> additional composition covariates (not deconvolution; does not replace the primary four tests)"
python3 src/run_composition_sensitivity.py

echo "==> additional UPR-arm split (does not replace STRESS_CORE_v1)"
python3 src/run_upr_arm_split.py

echo "==> detectable-effect checks (does not replace the primary four tests)"
python3 src/run_power_analysis.py

echo "==> figures"
python3 src/make_figures.py

echo "==> compare to expected_results"
python3 src/compare_to_expected.py

echo "Done. Outputs in $PWD/results"
