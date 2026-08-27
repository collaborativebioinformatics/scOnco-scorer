# =============================================================================
# Stage 3 - Scoring metrics
#
# Functions only - source() this, don't Rscript it.
# Validate with 04_validate_metrics.R BEFORE pointing it at real caller output.
#
# Base R only, no dependencies: this has to run inside four different
# caller containers without installing anything.
# =============================================================================

#' Area under the ROC curve, via the Mann-Whitney identity.
#' Handles ties correctly (average ranks). Returns NA if only one class present.
auroc <- function(y, p) {
  ok <- !is.na(y) & !is.na(p); y <- y[ok]; p <- p[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' Average precision (step-interpolated area under the PR curve).
#' Compare against the prevalence baseline, not against 0.5.
auprc <- function(y, p) {
  ok <- !is.na(y) & !is.na(p); y <- y[ok]; p <- p[ok]
  if (sum(y == 1) == 0 || sum(y == 0) == 0) return(NA_real_)
  o  <- order(p, decreasing = TRUE)
  ys <- y[o]
  tp <- cumsum(ys); fp <- cumsum(1 - ys)
  prec <- tp / (tp + fp)
  rec  <- tp / sum(ys)
  sum(diff(c(0, rec)) * prec)
}

#' Brier score (mean squared error of the probability). Lower is better.
brier <- function(y, p) {
  ok <- !is.na(y) & !is.na(p)
  mean((p[ok] - y[ok])^2)
}

#' Calibration curve: equal-width bins over [0,1].
#' Returns one row per bin with mean predicted prob (conf) and observed rate (acc).
calibration_curve <- function(y, p, n_bins = 10) {
  ok <- !is.na(y) & !is.na(p); y <- y[ok]; p <- p[ok]
  br  <- seq(0, 1, length.out = n_bins + 1)
  bin <- cut(p, breaks = br, include.lowest = TRUE, labels = FALSE)
  out <- data.frame(
    bin = seq_len(n_bins),
    bin_mid = (br[-1] + br[-length(br)]) / 2,
    n = 0L, conf = NA_real_, acc = NA_real_
  )
  for (i in seq_len(n_bins)) {
    idx <- which(bin == i)
    if (length(idx)) {
      out$n[i]    <- length(idx)
      out$conf[i] <- mean(p[idx])
      out$acc[i]  <- mean(y[idx])
    }
  }
  out
}

#' Expected calibration error: n-weighted mean |observed - predicted| across bins.
ece <- function(y, p, n_bins = 10) {
  cc <- calibration_curve(y, p, n_bins)
  n_tot <- sum(cc$n)
  if (n_tot == 0) return(NA_real_)
  sum((cc$n / n_tot) * abs(cc$acc - cc$conf), na.rm = TRUE)
}

#' Maximum calibration error - the worst single bin. Catches localised failures
#' that ECE averages away.
mce <- function(y, p, n_bins = 10) {
  cc <- calibration_curve(y, p, n_bins)
  cc <- cc[cc$n > 0, ]
  if (!nrow(cc)) return(NA_real_)
  max(abs(cc$acc - cc$conf), na.rm = TRUE)
}

#' Score one caller run against truth.
#'
#' @param y          integer 0/1 truth, one per cell
#' @param p          predicted probability of malignant, same length/order as y
#' @param threshold  cut for the hard-label metrics
#' @param n_bins     calibration bins
#' @param binary_only TRUE if the caller emits labels not probabilities
#'                    (CopyKAT, SCEVAN) - suppresses calibration metrics, which
#'                    are meaningless on degenerate 0/1 "probabilities"
#' @return list(metrics = one-row data.frame, calibration = data.frame)
score_calls <- function(y, p, threshold = 0.5, n_bins = 10, binary_only = FALSE) {
  stopifnot(length(y) == length(p))
  ok <- !is.na(y) & !is.na(p)
  n_dropped <- sum(!ok)
  y <- y[ok]; p <- p[ok]

  if (length(y) == 0) stop("No usable cells after dropping NAs - check the join on cell_id.")
  if (any(p < 0 | p > 1)) stop("Probabilities outside [0,1] - fix the adapter, don't clamp here.")

  pred <- as.integer(p >= threshold)
  tp <- sum(pred == 1 & y == 1); fp <- sum(pred == 1 & y == 0)
  fn <- sum(pred == 0 & y == 1); tn <- sum(pred == 0 & y == 0)

  sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  f1   <- if (!is.na(prec) && !is.na(sens) && (prec + sens) > 0) 2 * prec * sens / (prec + sens) else NA_real_

  metrics <- data.frame(
    n_cells      = length(y),
    n_dropped    = n_dropped,
    prevalence   = mean(y),
    auroc        = auroc(y, p),
    auprc        = auprc(y, p),
    brier        = if (binary_only) NA_real_ else brier(y, p),
    ece          = if (binary_only) NA_real_ else ece(y, p, n_bins),
    mce          = if (binary_only) NA_real_ else mce(y, p, n_bins),
    sensitivity  = sens,
    specificity  = spec,
    precision    = prec,
    f1           = f1,
    tp = tp, fp = fp, fn = fn, tn = tn,
    binary_only  = binary_only,
    stringsAsFactors = FALSE
  )

  list(metrics = metrics,
       calibration = if (binary_only) NULL else calibration_curve(y, p, n_bins))
}

#' Join caller output to truth and score. Reports the join loss loudly, because
#' silent barcode mismatch is the failure mode that ruins a benchmark.
#'
#' @param truth_df  data.frame with cell_id, y
#' @param calls_df  data.frame with cell_id, prob_malignant
score_run <- function(truth_df, calls_df, ...) {
  stopifnot(all(c("cell_id", "y") %in% names(truth_df)))
  stopifnot(all(c("cell_id", "prob_malignant") %in% names(calls_df)))

  m <- merge(truth_df[, c("cell_id", "y")], calls_df[, c("cell_id", "prob_malignant")],
             by = "cell_id", all.x = TRUE)

  n_truth   <- nrow(truth_df)
  n_matched <- sum(!is.na(m$prob_malignant))
  pct <- 100 * n_matched / n_truth

  if (pct < 90) {
    warning(sprintf("Only %.1f%% of truth cells matched caller output (%d/%d). Barcode format mismatch?",
                    pct, n_matched, n_truth))
    cat("  truth  e.g.: ", paste(head(truth_df$cell_id, 2), collapse = ", "), "\n")
    cat("  caller e.g.: ", paste(head(calls_df$cell_id, 2), collapse = ", "), "\n")
  }

  res <- score_calls(m$y, m$prob_malignant, ...)
  res$metrics$pct_matched <- round(pct, 2)
  res
}
