# =============================================================================
# benchmark_metrics_schema.R  -  the contract for benchmark_metrics.tsv
#
# Sourced by adapter_copykat.R (which WRITES the file) and by both figure
# scripts (which READ it), so the three can never disagree about its shape.
#
# Both figures read this file, so its shape is fixed here and validated before
# anything is written. Long format: ONE ROW PER METRIC per (caller, site,
# fraction, replicate). Wide format would need a new column every time a metric
# is added, and the figures would break each time.
#
# Every row carries its own interval and its own provenance, so a figure can
# never plot an estimate whose uncertainty came from somewhere else.
# =============================================================================

METRIC_COLS <- c(
  "caller",         # copykat | infercnv | scevan | numbat
  "site",           # LLU | NCI
  "platform",       # 10x | C1
  "fraction",       # nominal tumour fraction, 0-1; NA for controls
  "replicate",      # replicate id, or NA when replicates are pooled
  "metric",         # see KNOWN_METRICS
  "estimate",
  "lower", "upper", # 95% interval; NA where none was computed
  "n_cells",        # scored cells behind the estimate
  "n_source_cells", # distinct source cells; NA when not applicable
  "ci_method",      # cluster_bootstrap | none
  "n_boot", "seed", # NA when ci_method is none
  "preliminary",    # TRUE where the caveat must travel with the number
  "note",
  "figure"          # which figure plots this row, or "none"
)

# Splitting these matters: METHODS.md forbids calibration metrics on a caller
# that emits labels, so the schema records which family a metric belongs to.
BINARY_METRICS <- c("sensitivity", "specificity", "precision", "f1",
                    "frac_aneuploid", "frac_filtered", "false_aneuploid_rate",
                    "barcode_match_rate", "reference_misclass_rate")
CONTINUOUS_METRICS <- c("auroc", "auprc", "brier", "log_loss", "ece", "mce")
KNOWN_METRICS <- c(BINARY_METRICS, CONTINUOUS_METRICS)

# Metrics bounded on [0,1]. log_loss is not, so it is excluded.
UNIT_METRICS <- setdiff(KNOWN_METRICS, "log_loss")

# Which figure each row belongs to. Defined ONCE here rather than in each figure
# script: two copies of this rule would eventually drift, and a row would end up
# in both figures or in neither without anyone noticing.
CONTROL_METRICS <- c("false_aneuploid_rate", "barcode_match_rate",
                     "reference_misclass_rate")

#' Label every row with the figure that plots it, or "none".
#' Written into benchmark_metrics.tsv so the source-data files are transparently
#' a filtered view of it rather than a second, competing table.
assign_figure <- function(df) {
  is_control <- df$metric %in% CONTROL_METRICS |
                grepl("shuffle|B-vs-B", df$note, ignore.case = TRUE)
  # Per-replicate rows are kept for the replicate-consistency question
  # (METHODS.md 6A) but are not plotted: the figures show pooled estimates with
  # their intervals, and overlaying ten points per fraction would obscure them.
  per_rep <- !is.na(df$replicate)
  fig <- ifelse(per_rep & !is_control, "none",
         ifelse(is_control, "benchmark_controls",
         # p00 is the all-normal control run, not a purity point: sensitivity and
         # precision are undefined there and 0 has no place on a log axis
         ifelse(!is.na(df$fraction) & df$fraction > 0, "purity_performance", "none")))
  fig
}

#' Build rows for one scored run. Vectorised over metric.
metric_rows <- function(caller, site, platform, fraction, replicate,
                         metric, estimate, lower = NA_real_, upper = NA_real_,
                         n_cells = NA_integer_, n_source_cells = NA_integer_,
                         ci_method = "none", n_boot = NA_integer_, seed = NA_integer_,
                         preliminary = FALSE, note = "") {
  data.frame(caller = caller, site = site, platform = platform,
             fraction = fraction, replicate = replicate,
             metric = metric, estimate = estimate, lower = lower, upper = upper,
             n_cells = n_cells, n_source_cells = n_source_cells,
             ci_method = ci_method, n_boot = n_boot, seed = seed,
             preliminary = preliminary, note = note,
             stringsAsFactors = FALSE)
}

#' Convert a boot_summary() data frame into schema rows.
boot_summary_to_metric_rows <- function(bs, caller, site, platform, fraction,
                              replicate = NA_integer_, preliminary = FALSE,
                              note = "") {
  metric_rows(caller, site, platform, fraction, replicate,
               metric = bs$metric, estimate = bs$estimate,
               lower = bs$lower, upper = bs$upper,
               n_cells = NA_integer_, n_source_cells = bs$n_source_cells,
               ci_method = "cluster_bootstrap", n_boot = bs$n_boot, seed = bs$seed,
               preliminary = preliminary, note = note)
}

#' Check a results table before it is written or plotted.
#' Errors are structural; warnings are things a reader should know about.
validate_metrics <- function(df, strict = TRUE) {
  problems <- character(0)
  p <- function(...) problems <<- c(problems, paste0(...))

  miss <- setdiff(setdiff(METRIC_COLS, "figure"), names(df))
  if (length(miss)) p("missing column(s): ", paste(miss, collapse = ", "))
  if (length(problems)) stop(paste(problems, collapse = "\n"))

  unknown <- setdiff(unique(df$metric), KNOWN_METRICS)
  if (length(unknown))
    p("unknown metric(s): ", paste(unknown, collapse = ", "),
      "\n  add to KNOWN_METRICS or fix the writer")

  fr <- df$fraction[!is.na(df$fraction)]
  if (any(fr < 0 | fr > 1))
    p("fraction outside [0,1] - use 0.05, not 5")

  u <- df[df$metric %in% UNIT_METRICS & !is.na(df$estimate), ]
  if (any(u$estimate < 0 | u$estimate > 1))
    p("estimate outside [0,1] for a unit-interval metric")

  has_ci <- !is.na(df$lower) & !is.na(df$upper)
  if (any(df$lower[has_ci] > df$estimate[has_ci] + 1e-9) ||
      any(df$upper[has_ci] < df$estimate[has_ci] - 1e-9))
    p("interval does not contain its estimate on some rows")

  cb <- df$ci_method == "cluster_bootstrap"
  if (any(cb & (is.na(df$n_boot) | is.na(df$seed))))
    p("cluster_bootstrap rows must carry n_boot and seed (the seed ledger)")

  # calibration on a binary caller is exactly the mistake METHODS.md warns about
  bad_cal <- df$metric %in% c("brier", "ece", "mce", "log_loss") &
             df$caller %in% c("copykat", "scevan") & !is.na(df$estimate)
  if (any(bad_cal))
    p("calibration metric present for a label-only caller (",
      paste(unique(df$caller[bad_cal]), collapse = ", "),
      "). Pass binary_only = TRUE, or supply a continuous burden score.")

  key <- paste(df$caller, df$site, df$fraction, df$replicate, df$metric)
  if (any(duplicated(key)))
    p(sum(duplicated(key)), " duplicated (caller, site, fraction, replicate, metric) row(s)")

  if (length(problems)) {
    msg <- paste(problems, collapse = "\n")
    if (strict) stop(msg) else warning(msg)
  }

  # Advisory, not fatal: a single replicate cannot support a final interval,
  # so those rows should carry preliminary = TRUE.
  reps <- unique(df$replicate[!is.na(df$replicate)])
  if (length(reps) == 1 && any(!df$preliminary))
    message("NOTE: only replicate ", reps, " is present, but ",
            sum(!df$preliminary), " row(s) have preliminary = FALSE. ",
            "Single-replicate uncertainty should be flagged preliminary.")
  invisible(TRUE)
}

write_metrics <- function(df, path = "benchmark_metrics.tsv", strict = TRUE) {
  validate_metrics(df, strict = strict)
  df$figure <- assign_figure(df)
  df <- df[order(df$caller, df$site, df$fraction, df$replicate, df$metric), ]
  write.table(df[, METRIC_COLS], path, sep = "\t", row.names = FALSE,
              quote = FALSE, na = "NA")
  cat("wrote", path, "-", nrow(df), "rows,",
      length(unique(df$metric)), "metrics,",
      length(unique(df$caller)), "caller(s)\n")
  invisible(df)
}

read_metrics <- function(path = "benchmark_metrics.tsv", validate = TRUE) {
  df <- read.delim(path, stringsAsFactors = FALSE, na.strings = "NA")
  if (validate) validate_metrics(df, strict = FALSE)
  df
}
