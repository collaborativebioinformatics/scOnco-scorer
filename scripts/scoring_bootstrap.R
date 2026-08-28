# =============================================================================
# scoring_bootstrap.R  -  source-cell cluster bootstrap
#
# Implements METHODS.md section 5 (Uncertainty):
#   "Because source cells recur across mixture fractions and replicates,
#    ordinary row-level bootstrap resampling is invalid. Uncertainty is
#    estimated with at least 2,000 source-cell cluster-bootstrap replicates,
#    stratified by true A/B label and aggregated over mixture replicates.
#    Every principal per-fraction metric receives a 95% interval."
#
# WHY CLUSTERS AND NOT ROWS
#
# The in-silico sweep reuses one pool of source cells at every purity and in
# every replicate, and METHODS.md builds fractions from deterministic prefixes
# of a single ordering - so recurrence is by design, not accident. One source
# cell contributes up to (6 fractions x 10 replicates) = 60 rows. Those rows are
# the same cell measured again, not independent observations.
#
# Resampling ROWS treats each appearance as fresh, inflates the effective sample
# size, and produces intervals that are far too narrow. Resampling SOURCE CELLS
# keeps a cell's appearances together: it is either in a replicate with all its
# rows, or absent entirely.
#
# Requires: scoring_metrics.R
# =============================================================================

if (!exists("auroc")) source("scoring_metrics.R")

REQUIRED_COLS <- c("source_cell", "y", "p")

check_input <- function(df) {
  miss <- setdiff(REQUIRED_COLS, names(df))
  if (length(miss))
    stop("Missing column(s): ", paste(miss, collapse = ", "),
         "\n  'source_cell' must be the ORIGINAL pre-mixture barcode. If the",
         "\n  mixture output carries only mixture-local IDs, the cluster",
         "\n  bootstrap cannot run and any CI would be invalid.")
  if (!all(df$y %in% c(0, 1)))
    stop("y must be 0/1 (1 = HCC1395 / A, 0 = HCC1395BL / B).")
  if (any(df$p < 0 | df$p > 1, na.rm = TRUE))
    stop("p outside [0,1] - fix the adapter rather than clamping here.")
  invisible(TRUE)
}

#' Principal metrics in one pass. METHODS.md section 8 requires ROC-AUC, PR-AUC,
#' Brier, log loss, ECE, and sensitivity/specificity at prespecified operating
#' points. Returns a NAMED NUMERIC VECTOR so replicates stack into a matrix.
all_metrics <- function(y, p, threshold = 0.5, n_bins = 10) {
  c(auroc       = auroc(y, p),
    auprc       = auprc(y, p),
    brier       = brier(y, p),
    log_loss    = log_loss(y, p),
    ece         = ece(y, p, n_bins),
    mce         = mce(y, p, n_bins),
    sensitivity = if (sum(y == 1) > 0) mean(p[y == 1] >= threshold) else NA_real_,
    specificity = if (sum(y == 0) > 0) mean(p[y == 0] <  threshold) else NA_real_,
    prevalence  = mean(y))
}

#' Source-cell cluster bootstrap with a percentile interval.
#'
#' @param df       needs source_cell, y, p
#' @param n_boot   METHODS.md requires >= 2000
#' @param seed     recorded in the return value for the seed ledger (section 7)
#' @param stratify resample within true-label strata, per METHODS.md. Holds the
#'                 A:B ratio fixed across replicates, which matters enormously at
#'                 1% purity (30 A cells) where an unstratified draw can land on
#'                 very few - or zero - malignant cells.
cluster_bootstrap <- function(df, n_boot = 2000, seed = 1L, stratify = TRUE,
                              conf = 0.95, metric_fn = all_metrics,
                              verbose = TRUE) {
  check_input(df)
  df <- df[!is.na(df$p) & !is.na(df$y), , drop = FALSE]
  if (!nrow(df)) stop("No usable rows after dropping NA.")

  idx_by_cell <- split(seq_len(nrow(df)), df$source_cell)
  cells <- names(idx_by_cell)

  # A source cell has exactly one true label by construction; verify it.
  cell_label <- vapply(idx_by_cell, function(i) {
    u <- unique(df$y[i]); if (length(u) != 1L) NA_real_ else as.numeric(u)
  }, numeric(1))
  if (anyNA(cell_label))
    stop(sum(is.na(cell_label)), " source cell(s) carry more than one y value; ",
         "source_cell is not identifying a unique biological cell.")

  strata <- if (stratify) split(seq_along(cells), cell_label)
            else list(seq_along(cells))

  point <- metric_fn(df$y, df$p)
  reps <- matrix(NA_real_, n_boot, length(point),
                 dimnames = list(NULL, names(point)))

  set.seed(seed)
  for (b in seq_len(n_boot)) {
    drawn <- unlist(lapply(strata, function(s) sample(s, length(s), replace = TRUE)),
                    use.names = FALSE)
    rows <- unlist(idx_by_cell[cells[drawn]], use.names = FALSE)
    reps[b, ] <- metric_fn(df$y[rows], df$p[rows])
    if (verbose && b %% 500 == 0) message("  ", b, "/", n_boot)
  }

  a <- (1 - conf) / 2
  ci <- apply(reps, 2, quantile, probs = c(a, 1 - a), na.rm = TRUE)
  n_valid <- colSums(!is.na(reps))
  if (any(n_valid < 0.9 * n_boot))
    warning("Some metrics undefined in >10% of replicates (too few cells of one ",
            "class). Those intervals are unreliable.")

  list(point = point, lower = ci[1, ], upper = ci[2, ],
       n_valid = n_valid, replicates = reps, seed = seed, n_boot = n_boot,
       stratified = stratify, n_source_cells = length(cells), n_rows = nrow(df),
       rows_per_cell = round(nrow(df) / length(cells), 2))
}

#' Naive row bootstrap. Provided ONLY to demonstrate the understatement in the
#' test suite and in the writeup. Never report these intervals.
row_bootstrap <- function(df, n_boot = 2000, seed = 1L, conf = 0.95,
                          metric_fn = all_metrics) {
  check_input(df)
  df <- df[!is.na(df$p) & !is.na(df$y), , drop = FALSE]
  point <- metric_fn(df$y, df$p)
  reps <- matrix(NA_real_, n_boot, length(point), dimnames = list(NULL, names(point)))
  set.seed(seed)
  for (b in seq_len(n_boot)) {
    i <- sample(seq_len(nrow(df)), nrow(df), replace = TRUE)
    reps[b, ] <- metric_fn(df$y[i], df$p[i])
  }
  a <- (1 - conf) / 2
  ci <- apply(reps, 2, quantile, probs = c(a, 1 - a), na.rm = TRUE)
  list(point = point, lower = ci[1, ], upper = ci[2, ], replicates = reps)
}

boot_summary <- function(bs, ...) {
  extra <- list(...)
  out <- data.frame(metric = names(bs$point), estimate = as.numeric(bs$point),
                    lower = as.numeric(bs$lower), upper = as.numeric(bs$upper),
                    width = as.numeric(bs$upper - bs$lower),
                    n_valid = as.integer(bs$n_valid),
                    n_source_cells = bs$n_source_cells,
                    n_boot = bs$n_boot, seed = bs$seed, stringsAsFactors = FALSE)
  for (nm in names(extra)) out[[nm]] <- extra[[nm]]
  rownames(out) <- NULL
  out
}

#' Bootstrap separately at each tumour fraction. Replicates are POOLED within a
#' fraction, per METHODS.md ("aggregated over mixture replicates") - a source
#' cell recurring across replicates is one cluster, so pooling is safe at cell
#' level and would not be at row level.
boot_by_fraction <- function(df, n_boot = 2000, seed = 1L, ...) {
  if (!"fraction" %in% names(df)) stop("No 'fraction' column.")
  fr <- sort(unique(df$fraction))
  do.call(rbind, lapply(seq_along(fr), function(i) {
    sub <- df[df$fraction == fr[i], , drop = FALSE]
    message("fraction ", fr[i], ": ", nrow(sub), " rows, ",
            length(unique(sub$source_cell)), " source cells")
    bs <- cluster_bootstrap(sub, n_boot = n_boot, seed = seed + i, verbose = FALSE, ...)
    boot_summary(bs, fraction = fr[i],
                 n_replicates = if ("replicate" %in% names(sub))
                   length(unique(sub$replicate)) else NA_integer_)
  }))
}

#' Lowest tumour fraction where a metric clears a threshold BY ITS LOWER 95%
#' BOUND (METHODS.md section 6A). The point estimate is deliberately not used,
#' and non-overlap of adjacent intervals is NOT computed - METHODS.md rules it
#' out as a surrogate significance test.
detectability_floor <- function(summ, metric = "auroc", threshold = 0.80,
                                min_consecutive = 1L) {
  d <- summ[summ$metric == metric, , drop = FALSE]
  d <- d[order(d$fraction), , drop = FALSE]
  ok <- !is.na(d$lower) & d$lower >= threshold
  if (!any(ok)) return(NA_real_)
  r <- rle(ok)
  starts <- cumsum(c(1, head(r$lengths, -1)))
  good <- which(r$values & r$lengths >= min_consecutive)
  if (!length(good)) return(NA_real_)
  d$fraction[starts[good[1]]]
}

#' Label-shuffle control (METHODS.md section 5, hard controls). Labels are
#' permuted WITHIN source cell so the cluster structure survives; discrimination
#' must collapse to chance, with the interval covering 0.5.
label_shuffle_control <- function(df, n_boot = 500, seed = 99L) {
  check_input(df)
  cl <- unique(df[, c("source_cell", "y")])
  set.seed(seed)
  cl$y <- sample(cl$y)                       # permute at CELL level, not row
  sh <- merge(df[, setdiff(names(df), "y")], cl, by = "source_cell")
  cluster_bootstrap(sh, n_boot = n_boot, seed = seed, verbose = FALSE)
}
