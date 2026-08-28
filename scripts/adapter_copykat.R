# =============================================================================
# adapter_copykat.R  -  CopyKAT output -> benchmark_metrics.tsv rows
#
# Per run directory CopyKAT writes four TSVs. Two feed scoring:
#
#   *_predictions.tsv       one row per cell, ALREADY joined to truth
#                           cell_id, original_barcode, source, truth,
#                           tumor_fraction, replicate, copykat.pred,
#                           copykat_status, copykat_binary
#   *_run_summary.tsv       mixture_cells / predicted / filtered / aneuploid /
#                           diploid / not_defined
#   *_reference_summary.tsv reference_supplied / predicted / aneuploid /
#                           diploid / not_defined / absent
#   *_input_summary.tsv     cell and gene counts, reference_mixture_overlap
#
# CopyKAT emits LABELS, not probabilities, so every metric here is binary-only.
# Brier, ECE and log loss are deliberately absent: computing them off a 0/1 call
# returns a real-looking number that measures error rate, not calibration.
# A continuous burden score needs CNAmat out of *_copykat_full.rds; see
# cnv_burden() in scoring_mask.R.
#
# Requires: benchmark_metrics_schema.R, scoring_metrics.R
# =============================================================================

if (!exists("METRIC_COLS")) source("benchmark_metrics_schema.R")
if (!exists("auroc"))        source("scoring_metrics.R")

# CopyKAT returns three states, not two. "not.defined" is an indeterminate call,
# NOT a diploid one - counting it as diploid silently inflates specificity,
# because at low purity almost every cell is truly normal.
#   "filtered" (default) - excluded from scoring, counted in frac_filtered
#   "diploid"            - treated as a normal call (only if the team decides so)
NOT_DEFINED_POLICY <- "filtered"

#' Read one predictions TSV and normalise it.
read_copykat_predictions <- function(path, not_defined = NOT_DEFINED_POLICY) {
  d <- read.delim(path, stringsAsFactors = FALSE)
  need <- c("cell_id", "original_barcode", "source", "truth",
            "tumor_fraction", "replicate", "copykat.pred")
  miss <- setdiff(need, names(d))
  if (length(miss)) stop(basename(path), ": missing ", paste(miss, collapse = ", "))

  d$y <- ifelse(d$truth == "tumor", 1L, ifelse(d$truth == "normal", 0L, NA_integer_))
  if (anyNA(d$y)) stop(basename(path), ": unrecognised truth value(s): ",
                       paste(unique(d$truth[is.na(d$y)]), collapse = ", "))

  pred <- tolower(gsub("[^a-z]", "", tolower(d$copykat.pred)))
  d$call <- ifelse(pred == "aneuploid", 1L,
            ifelse(pred == "diploid",   0L, NA_integer_))   # not.defined -> NA
  if (not_defined == "diploid") d$call[is.na(d$call)] <- 0L

  # source_cell is the CLUSTER KEY for the bootstrap. Use cell_id, not
  # original_barcode: a raw 10x barcode is unique only within a source, so a
  # tumour and a normal cell can share one.
  d$source_cell <- d$cell_id
  d$fraction <- as.numeric(d$tumor_fraction)
  d$fraction[is.na(d$fraction)] <- 0        # p00 leaves the column blank
  d
}

#' Binary metrics for one run. Cells with an indeterminate call are excluded
#' from the confusion matrix and reported as frac_filtered instead.
copykat_metrics <- function(d) {
  scored <- d[!is.na(d$call), ]
  n_in   <- nrow(d)
  n_sc   <- nrow(scored)

  tp <- sum(scored$call == 1 & scored$y == 1); fp <- sum(scored$call == 1 & scored$y == 0)
  fn <- sum(scored$call == 0 & scored$y == 1); tn <- sum(scored$call == 0 & scored$y == 0)

  sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  f1   <- if (!is.na(prec) && !is.na(sens) && (prec + sens) > 0)
            2 * prec * sens / (prec + sens) else NA_real_

  c(sensitivity = sens, specificity = spec, precision = prec, f1 = f1,
    frac_aneuploid = if (n_sc > 0) mean(scored$call == 1) else NA_real_,
    frac_filtered  = if (n_in > 0) 1 - n_sc / n_in else NA_real_,
    # at 0% tumour every aneuploid call is false: this IS the B-vs-B control
    false_aneuploid_rate = if (all(scored$y == 0) && n_sc > 0)
                             mean(scored$call == 1) else NA_real_)
}

#' Cluster-bootstrap the binary metrics over source cells.
copykat_boot <- function(d, n_boot = 2000, seed = 1L) {
  if (!exists("cluster_bootstrap")) source("scoring_bootstrap.R")
  b <- d[!is.na(d$call), c("source_cell", "y")]
  b$p <- d$call[!is.na(d$call)]           # 0/1 treated as a degenerate score
  bs <- cluster_bootstrap(b, n_boot = n_boot, seed = seed, verbose = FALSE,
                          metric_fn = function(y, p)
                            c(sensitivity = if (sum(y == 1) > 0) mean(p[y == 1] == 1) else NA_real_,
                              specificity = if (sum(y == 0) > 0) mean(p[y == 0] == 0) else NA_real_,
                              precision   = if (sum(p == 1) > 0) mean(y[p == 1] == 1) else NA_real_,
                              frac_aneuploid = mean(p == 1)))
  bs
}

#' Everything for one run directory -> schema rows.
score_copykat_run <- function(dir, site = "LLU", platform = "10x",
                              n_boot = 2000, seed = 1L, preliminary = TRUE,
                              note = "") {
  f <- function(suffix) {
    hits <- list.files(dir, pattern = paste0(suffix, "$"), full.names = TRUE)
    if (!length(hits)) NULL else hits[1]
  }
  pf <- f("_predictions.tsv"); rf <- f("_run_summary.tsv"); ff <- f("_reference_summary.tsv")
  if (is.null(pf)) stop("no *_predictions.tsv in ", dir)

  d <- read_copykat_predictions(pf)
  fr  <- d$fraction[1]
  rep <- d$replicate[1]
  pt  <- copykat_metrics(d)

  bs <- copykat_boot(d, n_boot = n_boot, seed = seed)
  lo <- setNames(bs$lower, names(bs$point)); hi <- setNames(bs$upper, names(bs$point))

  rows <- list()
  for (nm in names(pt)) {
    if (is.na(pt[[nm]])) next
    has_ci <- nm %in% names(lo)
    rows[[length(rows) + 1]] <- metric_rows(
      "copykat", site, platform, fr, rep, nm, pt[[nm]],
      lower = if (has_ci) lo[[nm]] else NA_real_,
      upper = if (has_ci) hi[[nm]] else NA_real_,
      n_cells = nrow(d), n_source_cells = length(unique(d$source_cell)),
      ci_method = if (has_ci) "cluster_bootstrap" else "none",
      n_boot = if (has_ci) n_boot else NA_integer_,
      seed = if (has_ci) seed else NA_integer_,
      preliminary = preliminary,
      note = if (nm == "false_aneuploid_rate")
        paste(c("B-vs-B: p00, all-normal mixture", note), collapse = "; ") else note)
  }

  # reference misclassification comes from CopyKAT's own reference accounting
  if (!is.null(ff)) {
    rs <- read.delim(ff, stringsAsFactors = FALSE)
    if (rs$reference_predicted > 0)
      rows[[length(rows) + 1]] <- metric_rows(
        "copykat", site, platform, fr, rep, "reference_misclass_rate",
        rs$reference_aneuploid / rs$reference_predicted,
        n_cells = rs$reference_predicted, ci_method = "none",
        preliminary = preliminary,
        note = sprintf("held-out reference-B called aneuploid at fraction %.2f; %d not_defined, %d absent",
                       fr, rs$reference_not_defined, rs$reference_absent))
  }
  # Truth join: CopyKAT's predictions file already carries truth, so this is a
  # check that no cell was dropped between the mixture and the output, not a
  # fuzzy barcode match.
  rows[[length(rows) + 1]] <- metric_rows(
    "copykat", site, platform, fr, rep, "barcode_match_rate",
    mean(!is.na(d$y)), n_cells = nrow(d), ci_method = "none",
    preliminary = preliminary,
    note = "predictions.tsv ships joined to truth; checks for dropped cells")

  # Label shuffle: permute labels within source cell, rescore. Discrimination
  # must collapse to chance. Skipped where only one class is present (p00).
  if (length(unique(d$y)) == 2) {
    if (!exists("label_shuffle_control")) source("scoring_bootstrap.R")
    b <- d[!is.na(d$call), c("source_cell", "y")]
    b$p <- d$call[!is.na(d$call)]
    ls <- label_shuffle_control(b, n_boot = 500, seed = seed + 5000L)
    rows[[length(rows) + 1]] <- metric_rows(
      "copykat", site, platform, fr, rep, "auroc",
      ls$point[["auroc"]], ls$lower[["auroc"]], ls$upper[["auroc"]],
      n_cells = nrow(b), n_source_cells = length(unique(b$source_cell)),
      ci_method = "cluster_bootstrap", n_boot = 500L, seed = seed + 5000L,
      preliminary = preliminary, note = "label shuffle control")
  }

  do.call(rbind, rows)
}

#' Every run directory under a root, scored one at a time.
#' Use score_copykat_pooled() instead when replicates are available.
score_copykat_all <- function(root, pattern = "_r01$", ...) {
  dirs <- list.dirs(root, recursive = FALSE)
  dirs <- dirs[grepl(pattern, basename(dirs))]
  if (!length(dirs)) stop("no run directories matching ", pattern, " under ", root)
  out <- lapply(seq_along(dirs), function(i) {
    message("scoring ", basename(dirs[i]))
    score_copykat_run(dirs[i], seed = 1000L + i, ...)
  })
  do.call(rbind, out)
}

# =============================================================================
# Pooled scoring across replicates
#
# METHODS.md section 5: uncertainty is "aggregated over mixture replicates".
# Cells are pooled WITHIN a fraction and the cluster bootstrap resamples source
# cells across the pooled set, so a cell recurring in several replicates is one
# cluster rather than several independent observations.
#
# This is also the point at which the cluster bootstrap starts doing real work.
# With a single replicate no cell recurs, so cluster resampling is close to row
# resampling; across ten replicates each cell appears ten times and the
# intervals widen to reflect that.
# =============================================================================

#' Read every replicate at one fraction and stack them.
read_fraction <- function(dirs, not_defined = NOT_DEFINED_POLICY) {
  do.call(rbind, lapply(dirs, function(d) {
    pf <- list.files(d, pattern = "_predictions\\.tsv$", full.names = TRUE)
    if (!length(pf)) stop("no *_predictions.tsv in ", d)
    read_copykat_predictions(pf[1], not_defined)
  }))
}

#' Score one fraction with its replicates pooled.
score_fraction_pooled <- function(dirs, site = "LLU", platform = "10x",
                                  n_boot = 2000, seed = 1L, preliminary = FALSE,
                                  note = "") {
  d <- read_fraction(dirs)
  fr <- d$fraction[1]
  reps <- sort(unique(d$replicate))
  pt <- copykat_metrics(d)

  bs <- copykat_boot(d, n_boot = n_boot, seed = seed)
  lo <- setNames(bs$lower, names(bs$point)); hi <- setNames(bs$upper, names(bs$point))

  n_src <- length(unique(d$source_cell))
  note_full <- paste(c(sprintf("%d replicates pooled (%s); %.1f rows per source cell",
                               length(reps), paste(range(reps), collapse = "-"),
                               nrow(d) / n_src), note), collapse = "; ")

  rows <- list()
  for (nm in names(pt)) {
    if (is.na(pt[[nm]])) next
    has_ci <- nm %in% names(lo)
    rows[[length(rows) + 1]] <- metric_rows(
      "copykat", site, platform, fr, NA_integer_, nm, pt[[nm]],
      lower = if (has_ci) lo[[nm]] else NA_real_,
      upper = if (has_ci) hi[[nm]] else NA_real_,
      n_cells = nrow(d), n_source_cells = n_src,
      ci_method = if (has_ci) "cluster_bootstrap" else "none",
      n_boot = if (has_ci) n_boot else NA_integer_,
      seed = if (has_ci) seed else NA_integer_,
      preliminary = preliminary,
      note = if (nm == "false_aneuploid_rate")
        paste("B-vs-B: p00, all-normal mixture", note_full, sep = "; ") else note_full)
  }

  # Per-replicate values as well, for the replicate-consistency question:
  # "detected in k of n replicates" is a different statement from a pooled
  # estimate with a wide interval, and METHODS.md section 6A needs both.
  for (r in reps) {
    dr <- d[d$replicate == r, ]
    pr <- copykat_metrics(dr)
    for (nm in c("sensitivity", "specificity", "precision", "f1")) {
      if (is.na(pr[[nm]])) next
      rows[[length(rows) + 1]] <- metric_rows(
        "copykat", site, platform, fr, r, nm, pr[[nm]],
        n_cells = nrow(dr), n_source_cells = length(unique(dr$source_cell)),
        ci_method = "none", preliminary = preliminary,
        note = "per-replicate value; not pooled")
    }
  }

  # reference misclassification, averaged over the replicates at this fraction
  rf <- unlist(lapply(dirs, function(d0)
    list.files(d0, pattern = "_reference_summary\\.tsv$", full.names = TRUE)))
  if (length(rf)) {
    rs <- do.call(rbind, lapply(rf, read.delim, stringsAsFactors = FALSE))
    rows[[length(rows) + 1]] <- metric_rows(
      "copykat", site, platform, fr, NA_integer_, "reference_misclass_rate",
      sum(rs$reference_aneuploid) / sum(rs$reference_predicted),
      n_cells = sum(rs$reference_predicted), ci_method = "none",
      preliminary = preliminary,
      note = sprintf("held-out reference-B called aneuploid, %d replicates pooled",
                     nrow(rs)))
  }

  rows[[length(rows) + 1]] <- metric_rows(
    "copykat", site, platform, fr, NA_integer_, "barcode_match_rate",
    mean(!is.na(d$y)), n_cells = nrow(d), ci_method = "none",
    preliminary = preliminary,
    note = "predictions.tsv ships joined to truth; checks for dropped cells")

  if (length(unique(d$y)) == 2) {
    b <- d[!is.na(d$call), c("source_cell", "y")]
    b$p <- d$call[!is.na(d$call)]
    ls <- label_shuffle_control(b, n_boot = 500, seed = seed + 5000L)
    rows[[length(rows) + 1]] <- metric_rows(
      "copykat", site, platform, fr, NA_integer_, "auroc",
      ls$point[["auroc"]], ls$lower[["auroc"]], ls$upper[["auroc"]],
      n_cells = nrow(b), n_source_cells = length(unique(b$source_cell)),
      ci_method = "cluster_bootstrap", n_boot = 500L, seed = seed + 5000L,
      preliminary = preliminary, note = "label shuffle control")
  }
  do.call(rbind, rows)
}

#' All fractions, replicates pooled within each.
#' Run directories are expected to be named <site>_p<NN>_r<NN>.
score_copykat_pooled <- function(root, site = "LLU", platform = "10x",
                                 n_boot = 2000, seed_base = 1000L, ...) {
  dirs <- list.dirs(root, recursive = FALSE)
  dirs <- dirs[grepl("_p\\d+_r\\d+$", basename(dirs))]
  if (!length(dirs)) stop("no <site>_pNN_rNN directories under ", root)

  frac_tag <- sub("^.*_(p\\d+)_r\\d+$", "\\1", basename(dirs))
  by_frac <- split(dirs, frac_tag)
  by_frac <- by_frac[order(as.integer(sub("^p", "", names(by_frac))))]

  out <- lapply(seq_along(by_frac), function(i) {
    tag <- names(by_frac)[i]
    message("scoring ", tag, ": ", length(by_frac[[i]]), " replicate(s) pooled")
    score_fraction_pooled(by_frac[[i]], site = site, platform = platform,
                          n_boot = n_boot, seed = seed_base + i, ...)
  })
  do.call(rbind, out)
}
