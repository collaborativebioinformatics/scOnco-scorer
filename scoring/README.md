# Scoring metrics

Group-provided scoring utilities copied from:
`Group9_2026:/code/scripts/`

Files:
- `scoring_metrics.R` — AUROC, AUPRC, Brier score, ECE, calibration and classification metrics.
- `scoring_metrics_test.R` — smoke tests for the scoring functions.

## Smoke test

Validated on the DNAnexus Cloud Workstation on 2026-08-27.

All 16 checks passed.

The supplied test script currently contains:

    source("03_metrics.R")

while the DNAnexus scoring file is named:

    scoring_metrics.R

For validation, an uncommitted local symlink
`03_metrics.R -> scoring_metrics.R` was used. The group-provided scripts
were not modified.

## Caller scoring note

For categorical caller outputs, use `binary_only = TRUE`.

Continuous caller scores should be retained separately for ranking-based
evaluation and downstream calibration.
