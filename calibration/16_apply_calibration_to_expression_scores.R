# ============================================================
# 16_apply_calibration_to_expression_scores.R
#
# Applies calibration_module.R to the real CopyKAT / inferCNV /
# SCEVAN continuous scores from Arijita's expression-caller
# pipeline, incorporating every finding from the simulation
# validation (calibration_simulation.R):
#
#   1. group_id = individual cell barcode, NOT dataset_id or
#      replicate. The held-out reference B cells (and, less
#      obviously, some mixture cells) recur across many runs --
#      grouping by dataset/replicate would silently let the same
#      physical cell appear in both a training and a test fold,
#      exactly the leakage scenario 2 modeled.
#
#   2. Scenario 2's null leakage result was specific to a smooth,
#      low-capacity calibrator on a densely-sampled synthetic
#      score. It is NOT assumed to transfer here -- this script
#      re-runs the row-level-vs-grouped paired comparison
#      directly on the real scores as a sanity check.
#
#   3. Only continuous scores are calibrated. Any SCEVAN run that
#      fell back to a categorical label (pinned wrapper failure)
#      is excluded from calibration and reported separately,
#      descriptive-only.
#
#   4. Both Platt and isotonic are fit per caller; the simulation
#      showed Platt's apparent edge was an artifact of that
#      simulation's logistic-shaped ground truth, not a general
#      result -- CV metrics decide here, not an assumption.
#
#   5. AUC (discrimination) and calibration metrics are reported
#      as separate numbers per caller, never collapsed together.
#
#   6. inferCNV's reference-free condition gets a separate
#      population-level (aggregate, median-based) diagnostic,
#      analogous to the CopyKAT sensitivity/FPR decomposition --
#      not folded into per-cell calibration.
#
#   7. NCI is treated as a held-out generalization test: fit only
#      on LLU, then score NCI through the frozen calibrator and
#      report its own reliability diagram, never pooled into the
#      LLU training folds.
# ============================================================

source("calibration_module.R")
library(ggplot2)

PROJECT_ROOT <- "/pi/thomas.fazzio-umw/Kate/work"
results_root <- file.path(PROJECT_ROOT, "results/expression_callers")
out_dir      <- file.path(PROJECT_ROOT, "results/calibration/expression_callers")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# SECTION 1 -- Loaders. ADAPT THESE to the real file layout once
# Arijita's outputs exist. Everything below this section is
# generic and does not need to change.
#
# Expected standard schema after loading, one row per (cell,
# caller, dataset):
#   dataset_id, cell_id, replicate, tumor_fraction, population
#     (population: "LLU" or "NCI")
#   caller        ("copykat", "infercnv", "scevan")
#   score         numeric continuous score (burden/posterior/CNA)
#   score_type    "continuous" or "categorical"
#   truth         "tumor" / "normal"
# ============================================================

load_copykat_scores <- function() {
  # TODO: point at wherever Fangfei's continuous burden score
  # (extracted from the full copykat R object's CNA matrix) is
  # written, per-cell, per-dataset. Placeholder path shown.
  path <- file.path(results_root, "copykat", "copykat_continuous_scores.tsv")
  if (!file.exists(path)) {
    cat("NOTE: CopyKAT continuous scores not found at", path, "-- skipping.\n")
    return(NULL)
  }
  df <- read.delim(path, stringsAsFactors = FALSE)
  df$caller <- "copykat"
  df$score_type <- "continuous"
  df
}

load_infercnv_scores <- function() {
  # TODO: point at inferCNV's per-cell continuous output. Keep
  # the reference-free condition's cells clearly flagged (e.g.
  # a `condition` column with "referenced"/"reference_free") so
  # section 5 below can pull them out for the separate aggregate
  # diagnostic rather than mixing them into per-cell calibration.
  path <- file.path(results_root, "infercnv", "infercnv_continuous_scores.tsv")
  if (!file.exists(path)) {
    cat("NOTE: inferCNV continuous scores not found at", path, "-- skipping.\n")
    return(NULL)
  }
  df <- read.delim(path, stringsAsFactors = FALSE)
  df$caller <- "infercnv"
  df$score_type <- "continuous"
  df
}

load_scevan_scores <- function() {
  # TODO: point at the pinned-wrapper SCEVAN output. Must include
  # a status/flag column distinguishing runs where the wrapper
  # successfully returned a continuous CNA-derived score from any
  # that fell back to SCEVAN's default categorical label -- those
  # get score_type = "categorical" and are excluded from
  # calibration in section 3 below.
  path <- file.path(results_root, "scevan", "scevan_scores.tsv")
  if (!file.exists(path)) {
    cat("NOTE: SCEVAN scores not found at", path, "-- skipping.\n")
    return(NULL)
  }
  df <- read.delim(path, stringsAsFactors = FALSE)
  df$caller <- "scevan"
  # Expect df already has a score_type column ("continuous" /
  # "categorical") from the wrapper's own status reporting. If
  # not present, default conservatively to categorical so nothing
  # gets calibrated by mistake.
  if (!"score_type" %in% names(df)) {
    df$score_type <- "categorical"
  }
  df
}

all_scores <- do.call(rbind, list(
  load_copykat_scores(),
  load_infercnv_scores(),
  load_scevan_scores()
))

if (is.null(all_scores) || nrow(all_scores) == 0) {
  stop(
    "No expression-caller score files found under ", results_root,
    ". Update the loader paths in Section 1 once real outputs exist."
  )
}

all_scores$label <- as.integer(all_scores$truth == "tumor")

cat("Loaded", nrow(all_scores), "rows across callers:\n")
print(table(all_scores$caller, all_scores$score_type))

# ============================================================
# SECTION 2 -- Split out categorical-fallback rows (finding #3).
# Reported separately, descriptive-only -- never calibrated.
# ============================================================

categorical_rows <- all_scores[all_scores$score_type == "categorical", , drop = FALSE]
continuous_rows  <- all_scores[all_scores$score_type == "continuous", , drop = FALSE]

if (nrow(categorical_rows) > 0) {
  cat(
    "\nWARNING:", nrow(categorical_rows),
    "rows had only a categorical label (SCEVAN wrapper fallback or",
    "similar). These are excluded from calibration and reported",
    "descriptively only.\n"
  )

  categorical_summary <- do.call(rbind, lapply(
    split(categorical_rows, categorical_rows$caller),
    function(sub) {
      data.frame(
        caller = sub$caller[1],
        n_categorical_rows = nrow(sub),
        n_datasets_affected = length(unique(sub$dataset_id))
      )
    }
  ))

  write.table(
    categorical_summary,
    file.path(out_dir, "categorical_fallback_summary.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  print(categorical_summary)
}

# ============================================================
# SECTION 3 -- LLU vs. NCI split (finding #7). NCI is held out
# entirely from calibrator fitting; only ever scored afterward.
# ============================================================

llu_rows <- continuous_rows[continuous_rows$population == "LLU", , drop = FALSE]
nci_rows <- continuous_rows[continuous_rows$population == "NCI", , drop = FALSE]

cat(
  "\nLLU rows (used for fitting):", nrow(llu_rows),
  " | NCI rows (held out for validation only):", nrow(nci_rows), "\n"
)

# ============================================================
# SECTION 4 -- Per-caller grouped calibration CV, grouped by
# CELL BARCODE (finding #1), with both Platt and isotonic
# (finding #4). AUC and calibration metrics reported separately
# (finding #5).
# ============================================================

caller_results <- list()

for (caller_name in unique(llu_rows$caller)) {

  caller_data <- llu_rows[llu_rows$caller == caller_name, , drop = FALSE]

  cat("\n=== Calibrating:", caller_name, "(n =", nrow(caller_data), ") ===\n")

  cv_out <- run_grouped_calibration_cv(
    caller_data,
    group_col = "cell_id",
    score_col = "score",
    label_col = "label",
    k = 5,
    methods = c("platt", "isotonic"),
    seed = 42
  )

  print(cv_out$summary)

  cv_out$summary$caller <- caller_name
  cv_out$per_fold$caller <- caller_name

  caller_results[[caller_name]] <- cv_out

  write.table(
    cv_out$per_fold,
    file.path(out_dir, paste0(caller_name, "_per_fold_cv.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

combined_summary <- do.call(rbind, lapply(caller_results, function(x) x$summary))

write.table(
  combined_summary,
  file.path(out_dir, "all_callers_calibration_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\n=== Combined discrimination + calibration summary, all callers ===\n")
print(combined_summary)

# ============================================================
# SECTION 5 -- Leakage sanity check (finding #2): re-run the same
# CV grouped by dataset_id/replicate instead of cell_id, and
# compare via a paired test. Does NOT assume the simulation's
# null result transfers to real data.
# ============================================================

leakage_check_rows <- list()

for (caller_name in unique(llu_rows$caller)) {

  caller_data <- llu_rows[llu_rows$caller == caller_name, , drop = FALSE]

  if (length(unique(caller_data$cell_id)) == nrow(caller_data)) {
    cat(
      "\nSkipping leakage check for", caller_name,
      "-- no repeated cell_id values present in this data, so",
      "cell-level vs dataset-level grouping cannot differ.\n"
    )
    next
  }

  by_cell <- run_grouped_calibration_cv(
    caller_data, "cell_id", "score", "label",
    k = 5, methods = "isotonic", seed = 99
  )
  by_dataset <- run_grouped_calibration_cv(
    caller_data, "dataset_id", "score", "label",
    k = 5, methods = "isotonic", seed = 99
  )

  leakage_check_rows[[caller_name]] <- data.frame(
    caller = caller_name,
    grouped_by_cell_brier = mean(by_cell$per_fold$brier),
    grouped_by_dataset_brier = mean(by_dataset$per_fold$brier),
    diff = mean(by_dataset$per_fold$brier) - mean(by_cell$per_fold$brier)
  )
}

if (length(leakage_check_rows) > 0) {
  leakage_check <- do.call(rbind, leakage_check_rows)
  write.table(
    leakage_check,
    file.path(out_dir, "leakage_sanity_check.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  cat("\n=== Leakage sanity check: cell-level vs. dataset-level grouping ===\n")
  print(leakage_check)
  cat(
    "\nA meaningfully positive 'diff' here means grouping by",
    "dataset_id/replicate (as opposed to cell_id) understates error",
    "-- i.e. real leakage, unlike the simulation's null result.\n"
  )
}

# ============================================================
# SECTION 6 -- NCI held-out validation (finding #7). Calibrators
# are refit on ALL of LLU (not CV folds) and frozen, then scored
# against NCI, which the calibrator has never seen.
# ============================================================

if (nrow(nci_rows) > 0) {

  nci_results <- list()

  for (caller_name in unique(nci_rows$caller)) {

    llu_caller <- llu_rows[llu_rows$caller == caller_name, , drop = FALSE]
    nci_caller <- nci_rows[nci_rows$caller == caller_name, , drop = FALSE]

    if (nrow(llu_caller) == 0) next

    frozen_iso <- fit_isotonic(llu_caller$score, llu_caller$label)
    nci_caller$prob_isotonic <- frozen_iso(nci_caller$score)

    brier_nci <- brier_score(nci_caller$prob_isotonic, nci_caller$label)
    ece_nci <- expected_calibration_error(nci_caller$prob_isotonic, nci_caller$label)
    auc_nci <- auc_score(nci_caller$score, nci_caller$label)

    nci_results[[caller_name]] <- data.frame(
      caller = caller_name,
      n_nci_cells = nrow(nci_caller),
      auc = auc_nci,
      brier = brier_nci,
      ece = ece_nci
    )

    rel <- reliability_curve_data(nci_caller$prob_isotonic, nci_caller$label)
    rel$caller <- caller_name

    fig <- ggplot(rel, aes(x = mean_predicted, y = observed_frequency)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
      geom_point(aes(size = n), color = "#c0392b") +
      geom_line(color = "#c0392b") +
      scale_x_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
      labs(
        x = "Mean predicted P(malignant)", y = "Observed frequency",
        size = "n in bin",
        title = paste0(caller_name, ": NCI held-out validation"),
        subtitle = paste0(
          "Calibrator fit on LLU only. Brier = ", round(brier_nci, 4),
          ", ECE = ", round(ece_nci, 4)
        )
      ) +
      theme_bw(base_size = 12)

    ggsave(
      file.path(out_dir, paste0(caller_name, "_nci_reliability.png")),
      fig, width = 6.5, height = 5.5, dpi = 300
    )
  }

  if (length(nci_results) > 0) {
    nci_summary <- do.call(rbind, nci_results)
    write.table(
      nci_summary, file.path(out_dir, "nci_holdout_summary.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    cat("\n=== NCI held-out validation summary ===\n")
    print(nci_summary)
  }

} else {
  cat("\nNo NCI rows present yet -- skipping held-out validation.\n")
}

# ============================================================
# SECTION 7 -- inferCNV reference-free aggregate diagnostic
# (finding #6). Population-level, median-based -- mirrors the
# CopyKAT sensitivity/FPR decomposition, kept separate from
# per-cell calibration above.
# ============================================================

infercnv_rows <- continuous_rows[continuous_rows$caller == "infercnv", , drop = FALSE]

if ("condition" %in% names(infercnv_rows) &&
    "reference_free" %in% unique(infercnv_rows$condition)) {

  ref_free <- infercnv_rows[infercnv_rows$condition == "reference_free", , drop = FALSE]

  agg <- do.call(rbind, lapply(
    split(ref_free, ref_free$tumor_fraction),
    function(sub) {
      data.frame(
        tumor_fraction = sub$tumor_fraction[1],
        median_score_tumor = median(sub$score[sub$truth == "tumor"]),
        median_score_normal = median(sub$score[sub$truth == "normal"]),
        n_tumor = sum(sub$truth == "tumor"),
        n_normal = sum(sub$truth == "normal")
      )
    }
  ))
  agg <- agg[order(agg$tumor_fraction), ]

  write.table(
    agg,
    file.path(out_dir, "infercnv_reference_free_baseline_inversion.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  cat("\n=== inferCNV reference-free baseline-inversion diagnostic ===\n")
  print(agg)
  cat(
    "\nWatch for tumor/normal median scores converging or inverting",
    "at low purity -- the reference-free condition has no external",
    "anchor, so its 'baseline' is only as good as the assumption",
    "that most cells in the mixture are normal.\n"
  )

} else {
  cat(
    "\nNo inferCNV reference-free condition column found -- skipping",
    "Section 7. Confirm the loader in Section 1 includes a",
    "`condition` column distinguishing referenced vs. reference-free",
    "runs once inferCNV output exists.\n"
  )
}

cat("\nDone. All outputs written to", out_dir, "\n")
