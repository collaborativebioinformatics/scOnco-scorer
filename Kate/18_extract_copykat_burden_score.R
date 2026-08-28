# ============================================================
# 18_extract_copykat_burden_score.R
#
# Extracts a continuous per-cell CNA-burden score from each
# CopyKAT run's full R object ($CNAmat -- the smoothed relative
# copy-number matrix, genomic bins x cells), joins it against
# the manifest for truth/tumor_fraction/replicate, and writes it
# in exactly the schema 16_apply_calibration_to_expression_scores.R
# expects, at exactly the path its loader looks for -- so once
# this runs, 16's CopyKAT section is unblocked with no further
# edits needed.
#
# Prerequisite: 07b_pull_copykat_full_objects.sh has been run,
# since 07_pull_predictions_dnanexus.sh deliberately skips these
# (large, not needed for the earlier run-level calibration work).
#
# Burden metric: CopyKAT's CNAmat is a log2-ratio-like relative
# copy-number estimate per genomic bin. A cell with real copy-
# number alterations deviates from ~0 across many bins; a normal
# diploid cell stays close to ~0 throughout. Three standard
# summary statistics are computed per cell so downstream
# calibration can pick whichever correlates best with truth via
# CV, rather than assuming one upfront:
#   - burden_mean_abs : mean(|value|) across all bins
#   - burden_rms       : sqrt(mean(value^2)) across all bins
#   - burden_var        : var(value) across all bins
# ============================================================

PROJECT_ROOT  <- "/pi/thomas.fazzio-umw/Kate/work"
copykat_root  <- file.path(PROJECT_ROOT, "results/copykat")
manifest_path <- file.path(PROJECT_ROOT, "manifests/LLU_all_mixtures_with_p00.tsv")
ref_cells_path <- file.path(PROJECT_ROOT, "metadata/LLU_B_reference_cells.txt")

# Output path matches exactly what 16's load_copykat_scores()
# looks for -- results/expression_callers/copykat/... under
# PROJECT_ROOT, not results/calibration/...
out_dir <- file.path(PROJECT_ROOT, "results/expression_callers/copykat")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(out_dir, "copykat_continuous_scores.tsv")

master_manifest <- read.delim(manifest_path, stringsAsFactors = FALSE)

ref_cells <- if (file.exists(ref_cells_path)) readLines(ref_cells_path) else character(0)

dataset_ids <- list.dirs(copykat_root, full.names = FALSE, recursive = FALSE)

if (length(dataset_ids) == 0) {
  stop("No dataset directories found under ", copykat_root)
}

cat("Found", length(dataset_ids), "dataset directories.\n")

# ------------------------------------------------------------
# Per-dataset extraction
# ------------------------------------------------------------

extract_dataset <- function(dataset_id) {

  run_dir <- file.path(copykat_root, dataset_id)
  full_rds <- file.path(run_dir, paste0(dataset_id, "_copykat_full.rds"))
  done_file <- file.path(run_dir, "DONE")

  if (!file.exists(done_file)) {
    cat("SKIP (no DONE marker):", dataset_id, "\n")
    return(NULL)
  }
  if (!file.exists(full_rds)) {
    cat(
      "SKIP (no full R object locally -- run",
      "07b_pull_copykat_full_objects.sh first):", dataset_id, "\n"
    )
    return(NULL)
  }

  ck <- readRDS(full_rds)

  if (!"CNAmat" %in% names(ck)) {
    cat(
      "SKIP (no CNAmat in R object -- unexpected CopyKAT output",
      "structure):", dataset_id, "\n"
    )
    return(NULL)
  }

  cna <- ck$CNAmat

  meta_cols <- intersect(c("chrom", "chrompos", "abspos"), colnames(cna))

  if (length(meta_cols) == 0) {
    cat(
      "WARNING:", dataset_id, "-- expected metadata columns",
      "(chrom/chrompos/abspos) not found in CNAmat. Check CopyKAT",
      "version/output structure. Skipping this dataset.\n"
    )
    return(NULL)
  }

  cell_cols <- setdiff(colnames(cna), meta_cols)

  if (length(cell_cols) == 0) {
    cat("SKIP (no cell columns found in CNAmat):", dataset_id, "\n")
    return(NULL)
  }

  cna_matrix <- as.matrix(cna[, cell_cols, drop = FALSE])
  mode(cna_matrix) <- "numeric"

  burden_mean_abs <- colMeans(abs(cna_matrix), na.rm = TRUE)
  burden_rms <- sqrt(colMeans(cna_matrix^2, na.rm = TRUE))
  burden_var <- apply(cna_matrix, 2, var, na.rm = TRUE)

  burden_df <- data.frame(
    dataset_id = dataset_id,
    cell_id = cell_cols,
    burden_mean_abs = burden_mean_abs,
    burden_rms = burden_rms,
    burden_var = burden_var,
    n_bins = nrow(cna_matrix),
    stringsAsFactors = FALSE
  )

  # ------------------------------------------------------------
  # Join truth/tumor_fraction/replicate. Two sources of truth:
  # mixture cells (from the per-dataset manifest slice) and
  # reference cells (fixed list, always truth = "normal", always
  # supplied to every run as norm.cell.names but NOT part of the
  # scored 1000-cell mixture in the manifest).
  # ------------------------------------------------------------

  dataset_manifest <- master_manifest[
    master_manifest$dataset_id == dataset_id, , drop = FALSE
  ]

  mixture_join <- merge(
    burden_df,
    dataset_manifest[, c("cell_id", "truth", "tumor_fraction", "replicate", "source")],
    by = "cell_id",
    all.x = FALSE
  )

  ref_in_this_run <- burden_df[burden_df$cell_id %in% ref_cells, , drop = FALSE]

  if (nrow(ref_in_this_run) > 0) {
    ref_join <- ref_in_this_run
    ref_join$truth <- "normal"
    ref_join$tumor_fraction <- unique(dataset_manifest$tumor_fraction)[1]
    ref_join$replicate <- unique(dataset_manifest$replicate)[1]
    ref_join$source <- "LLU_B_reference"

    ref_join <- ref_join[, names(mixture_join)]
  } else {
    ref_join <- mixture_join[0, ]
  }

  n_unmatched <- nrow(burden_df) - nrow(mixture_join) - nrow(ref_join)

  if (n_unmatched > 0) {
    cat(
      "  NOTE:", dataset_id, "--", n_unmatched,
      "cells in CNAmat had no matching manifest row or reference",
      "entry (not mixture, not reference -- likely filtered before",
      "prediction or a barcode-format mismatch worth checking).\n"
    )
  }

  rbind(mixture_join, ref_join)
}

all_results <- lapply(dataset_ids, extract_dataset)
all_results <- all_results[!vapply(all_results, is.null, logical(1))]

if (length(all_results) == 0) {
  stop(
    "No datasets successfully extracted. Most likely cause: full",
    " R objects not pulled locally yet -- run",
    " 07b_pull_copykat_full_objects.sh first."
  )
}

combined <- do.call(rbind, all_results)

cat(
  "\nExtracted burden scores for", nrow(combined), "cells across",
  length(unique(combined$dataset_id)), "datasets.\n"
)

# ------------------------------------------------------------
# Which burden metric best separates truth? Quick AUC check per
# metric -- informational only, does not choose a "score" column
# for downstream calibration (16 will fit/evaluate whichever is
# passed in; that decision should go through the same CV
# machinery, not be hardcoded here).
# ------------------------------------------------------------

quick_auc <- function(x, truth) {
  labels <- as.integer(truth == "tumor")
  n1 <- sum(labels == 1); n0 <- sum(labels == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(x)
  (sum(r[labels == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

cat("\nQuick informational AUC per candidate burden metric (pooled, all cells):\n")
cat("  burden_mean_abs:", round(quick_auc(combined$burden_mean_abs, combined$truth), 3), "\n")
cat("  burden_rms:", round(quick_auc(combined$burden_rms, combined$truth), 3), "\n")
cat("  burden_var:", round(quick_auc(combined$burden_var, combined$truth), 3), "\n")
cat(
  "(This is a naive pooled check, not grouped/CV'd -- use",
  "16_apply_calibration_to_expression_scores.R's grouped CV for",
  "the real per-caller evaluation, not this number.)\n\n"
)

# ------------------------------------------------------------
# Write in the exact schema 16's load_copykat_scores() expects.
# `score` defaults to burden_rms; edit this line if the AUC check
# above suggests a different metric is more informative, or leave
# all three columns in and let 16's loader pick.
# ------------------------------------------------------------

combined$population <- "LLU"
combined$score <- combined$burden_rms
combined$score_type <- "continuous"

output_table <- combined[, c(
  "dataset_id", "cell_id", "replicate", "tumor_fraction", "population",
  "score", "score_type", "truth",
  "burden_mean_abs", "burden_rms", "burden_var", "source"
)]

write.table(
  output_table, out_path,
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("Written to", out_path, "\n")
cat("16_apply_calibration_to_expression_scores.R can now be run\n")
cat("as-is for the CopyKAT caller -- no loader edits needed.\n")
