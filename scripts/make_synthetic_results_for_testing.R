#!/usr/bin/env Rscript
# =============================================================================
# Synthetic results_summary.tsv, so the figure scripts can be written and tested
# before any caller output exists. Numbers here are INVENTED - the point is to
# exercise the schema and the plotting, never to stand in for a result.
# =============================================================================
source("results_summary_schema.R")

FRACTIONS <- c(0.01, 0.05, 0.10, 0.20, 0.40, 0.80)
set.seed(20260827)

rows <- list()
add <- function(...) rows[[length(rows) + 1]] <<- results_rows(...)

for (site in c("LLU", "NCI")) {
  off <- if (site == "NCI") -0.045 else 0    # sites differ slightly, as they would
  for (f in FRACTIONS) {
    # discrimination improves with purity; intervals widen as it falls
    sens <- plogis(3.2 * log10(f / 0.01) - 1.1) + off
    spec <- 0.93 + 0.05 * f + off / 3
    prec <- plogis(4.0 * f - 0.4) + off
    f1   <- 2 * prec * sens / (prec + sens)
    w    <- 0.16 * (1 - f) + 0.02
    for (nm in c("sensitivity", "specificity", "precision", "f1")) {
      est <- switch(nm, sensitivity = sens, specificity = spec, precision = prec, f1 = f1)
      add("copykat", site, "10x", f, NA_integer_, nm, est,
          max(0, est - w), min(1, est + w),
          n_cells = 3000L, n_source_cells = 1600L,
          ci_method = "cluster_bootstrap", n_boot = 2000L, seed = 100L,
          preliminary = TRUE, note = "single replicate; preliminary")
    }
    add("copykat", site, "10x", f, NA_integer_, "frac_aneuploid",
        min(1, max(0, f * 0.86 + 0.02)), max(0, f * 0.86 - 0.04), min(1, f * 0.86 + 0.08),
        n_cells = 3000L, n_source_cells = 1600L,
        ci_method = "cluster_bootstrap", n_boot = 2000L, seed = 100L,
        preliminary = TRUE, note = "single replicate; preliminary")
    add("copykat", site, "10x", f, NA_integer_, "frac_filtered",
        0.06 + 0.03 * (1 - f), NA_real_, NA_real_, n_cells = 3000L,
        preliminary = TRUE, note = "single replicate; preliminary")
    # continuous burden score from the masked CNAmat
    auroc <- plogis(2.6 * log10(f / 0.01) - 0.2) * 0.45 + 0.55 + off / 2
    add("copykat", site, "10x", f, NA_integer_, "auroc", auroc,
        max(0.5, auroc - w), min(1, auroc + w), n_cells = 3000L, n_source_cells = 1600L,
        ci_method = "cluster_bootstrap", n_boot = 2000L, seed = 100L,
        preliminary = TRUE, note = "burden score; single replicate")
    add("copykat", site, "10x", f, NA_integer_, "auprc",
        min(1, auroc * f + 0.05), NA_real_, NA_real_, n_cells = 3000L,
        preliminary = TRUE, note = "burden score; single replicate")
  }

  # --- controls -------------------------------------------------------------
  add("copykat", site, "10x", NA_real_, NA_integer_, "false_aneuploid_rate",
      0.031, 0.019, 0.047, n_cells = 1400L, n_source_cells = 1400L,
      ci_method = "cluster_bootstrap", n_boot = 2000L, seed = 200L,
      preliminary = TRUE, note = "B-vs-B: held-out scored-B vs disjoint reference-B")
  add("copykat", site, "10x", NA_real_, NA_integer_, "auroc",
      0.503, 0.472, 0.535, n_cells = 3000L, n_source_cells = 1600L,
      ci_method = "cluster_bootstrap", n_boot = 2000L, seed = 300L,
      preliminary = TRUE, note = "label shuffle")
  add("copykat", site, "10x", NA_real_, NA_integer_, "barcode_match_rate",
      0.998, NA_real_, NA_real_, n_cells = 3000L, note = "truth join on source_cell")
  add("copykat", site, "10x", NA_real_, NA_integer_, "reference_misclass_rate",
      0.012, NA_real_, NA_real_, n_cells = 1400L,
      note = "reference-B cells called aneuploid")
}

df <- do.call(rbind, rows)
df$note[df$metric == "auroc" & is.na(df$fraction)] <- "label shuffle"
write_results(df, "results_summary_SYNTHETIC.tsv")
cat("\nSYNTHETIC DATA - invented numbers for testing the figures only.\n")
