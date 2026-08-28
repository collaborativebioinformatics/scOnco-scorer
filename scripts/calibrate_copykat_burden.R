#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (!length(hit)) return(default)
  i <- hit[length(hit)]
  if (i == length(args)) stop(flag, " requires a value", call. = FALSE)
  args[i + 1L]
}

input <- get_arg("--input", "/w/copykat_burden_percell.tsv")
outdir <- get_arg("--outdir", "/w")
group_cols_arg <- get_arg("--group_cols", NA_character_)
seed <- as.integer(get_arg("--seed", "100"))
outer_k <- as.integer(get_arg("--outer_k", "5"))
inner_k <- as.integer(get_arg("--inner_k", "5"))
ece_bins <- as.integer(get_arg("--ece_bins", "10"))
target_npv <- as.numeric(get_arg("--target_npv", "0.95"))
target_ppv <- as.numeric(get_arg("--target_ppv", "0.95"))

if (!file.exists(input)) stop("input does not exist: ", input, call. = FALSE)
if (is.na(seed)) stop("--seed must be an integer", call. = FALSE)
if (is.na(outer_k) || outer_k < 2L) stop("--outer_k must be >= 2", call. = FALSE)
if (is.na(inner_k) || inner_k < 2L) stop("--inner_k must be >= 2", call. = FALSE)
if (is.na(ece_bins) || ece_bins < 2L) stop("--ece_bins must be >= 2", call. = FALSE)
if (!is.finite(target_npv) || target_npv <= 0 || target_npv > 1) {
  stop("--target_npv must be in (0,1]", call. = FALSE)
}
if (!is.finite(target_ppv) || target_ppv <= 0 || target_ppv > 1) {
  stop("--target_ppv must be in (0,1]", call. = FALSE)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

d <- fread(input, showProgress = FALSE)

required_base <- c("burden", "copykat_scored", "tumor_fraction", "replicate")
missing_base <- setdiff(required_base, names(d))
if (length(missing_base)) {
  stop(
    "input missing required columns: ",
    paste(missing_base, collapse = ", "),
    call. = FALSE
  )
}

if (!"y" %in% names(d)) {
  if (!"truth" %in% names(d)) {
    stop("input must contain either y or truth", call. = FALSE)
  }
  z <- tolower(trimws(as.character(d$truth)))
  d[, y := NA_integer_]
  d[z %in% c("tumor", "tumour", "a", "hcc1395"), y := 1L]
  d[z %in% c("normal", "b", "hcc1395bl"), y := 0L]
}

if (is.logical(d$copykat_scored)) {
  scored <- d$copykat_scored
} else {
  s <- tolower(trimws(as.character(d$copykat_scored)))
  scored <- s %in% c("true", "t", "1", "yes", "y")
}

d <- d[scored & as.numeric(tumor_fraction) > 0]
if (!nrow(d)) stop("no scored rows with tumor_fraction > 0", call. = FALSE)

d[, burden := as.numeric(burden)]
d[, y := as.integer(y)]
d[, tumor_fraction := as.numeric(tumor_fraction)]
d[, replicate := as.character(replicate)]
d[, .row_id := .I]

if (anyNA(d$burden) || any(!is.finite(d$burden))) {
  stop("burden contains NA/non-finite values", call. = FALSE)
}
if (anyNA(d$y) || any(!d$y %in% c(0L, 1L))) {
  stop("y must contain only 0/1", call. = FALSE)
}

if (is.na(group_cols_arg) || !nzchar(trimws(group_cols_arg))) {
  stop(
    paste0(
      "--group_cols is required. Supply the stable biological source-cell identifier ",
      "used across mixture contexts, for example ",
      "--group_cols source_barcode,source_replicate. ",
      "Do not use a mixture/context-specific cell ID."
    ),
    call. = FALSE
  )
}

group_cols <- trimws(strsplit(group_cols_arg, ",", fixed = TRUE)[[1L]])
group_cols <- group_cols[nzchar(group_cols)]
if (!length(group_cols)) stop("--group_cols resolved to zero columns", call. = FALSE)

missing_group_cols <- setdiff(group_cols, names(d))
if (length(missing_group_cols)) {
  stop(
    "input is missing requested grouping column(s): ",
    paste(missing_group_cols, collapse = ", "),
    ". Regenerate copykat_burden_percell.tsv preserving original source-cell metadata.",
    call. = FALSE
  )
}

for (cc in group_cols) {
  if (anyNA(d[[cc]]) || any(!nzchar(trimws(as.character(d[[cc]]))))) {
    stop("grouping column contains NA/empty values: ", cc, call. = FALSE)
  }
}

d[, gid := do.call(paste, c(.SD, sep = "||")), .SDcols = group_cols]

group_truth <- d[, .(n_truth = uniqueN(y), y = first(y)), by = gid]
if (any(group_truth$n_truth != 1L)) {
  stop("at least one source-cell group has inconsistent truth labels", call. = FALSE)
}

group_context <- d[
  ,
  .(
    n_rows = .N,
    n_contexts = uniqueN(paste(tumor_fraction, replicate, sep = "|"))
  ),
  by = gid
]

if (uniqueN(d$tumor_fraction) > 1L && !any(group_context$n_contexts > 1L)) {
  stop(
    paste0(
      "the selected source-cell identifier never recurs across mixture contexts. ",
      "This is incompatible with the stated source-cell-grouped nested design. ",
      "Use the original source barcode plus source replicate/capture identifier."
    ),
    call. = FALSE
  )
}

make_group_folds <- function(dt, k, seed_value) {
  g <- unique(dt[, .(gid, y)])
  counts <- g[, .N, by = y][order(y)]

  if (nrow(counts) != 2L) {
    stop("fold construction requires both classes", call. = FALSE)
  }

  if (any(counts$N < k)) {
    stop(
      "not enough source-cell groups per class for ",
      k,
      "-fold CV: ",
      paste0("class ", counts$y, "=", counts$N, collapse = ", "),
      call. = FALSE
    )
  }

  set.seed(seed_value)
  g[, fold := NA_integer_]

  for (cls in sort(unique(g$y))) {
    idx <- which(g$y == cls)
    g$fold[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }

  setNames(g$fold, g$gid)
}

brier <- function(p, y) {
  ok <- is.finite(p) & !is.na(y)
  if (!any(ok)) return(NA_real_)
  mean((p[ok] - y[ok])^2)
}

logloss <- function(p, y) {
  ok <- is.finite(p) & !is.na(y)
  if (!any(ok)) return(NA_real_)
  p <- p[ok]
  y <- y[ok]
  p <- pmin(pmax(p, 1e-15), 1 - 1e-15)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

ece <- function(p, y, B = 10L) {
  ok <- is.finite(p) & !is.na(y)
  p <- p[ok]
  y <- y[ok]

  if (!length(p)) return(NA_real_)

  br <- cut(
    p,
    breaks = seq(0, 1, length.out = B + 1L),
    include.lowest = TRUE,
    right = TRUE
  )

  tab <- data.table(p = p, y = y, br = br)[
    ,
    .(
      n = .N,
      mean_probability = mean(p),
      observed_fraction = mean(y)
    ),
    by = br
  ]

  sum(
    tab$n / length(p) *
      abs(tab$observed_fraction - tab$mean_probability)
  )
}

auroc <- function(score, y) {
  ok <- is.finite(score) & !is.na(y)
  score <- score[ok]
  y <- y[ok]

  if (length(unique(y)) != 2L) return(NA_real_)

  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)

  if (n1 == 0L || n0 == 0L) return(NA_real_)

  r <- rank(score, ties.method = "average")
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

average_precision <- function(score, y) {
  ok <- is.finite(score) & !is.na(y)
  score <- score[ok]
  y <- y[ok]

  if (length(unique(y)) != 2L) return(NA_real_)

  n_positive <- sum(y == 1L)
  if (n_positive == 0L) return(NA_real_)

  z <- data.table(score = score, y = y)[
    ,
    .(
      positives = sum(y == 1L),
      negatives = sum(y == 0L)
    ),
    by = score
  ][order(-score)]

  z[, tp := cumsum(positives)]
  z[, fp := cumsum(negatives)]
  z[, recall := tp / n_positive]
  z[, precision := tp / (tp + fp)]

  sum(diff(c(0, z$recall)) * z$precision)
}

fit_platt <- function(score, y) {
  if (length(unique(y)) != 2L) return(NULL)

  fit <- tryCatch(
    suppressWarnings(
      glm(
        y ~ score,
        data = data.frame(y = y, score = score),
        family = binomial()
      )
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NULL)

  cf <- coef(fit)
  if (length(cf) != 2L || any(!is.finite(cf))) return(NULL)

  list(type = "platt", fit = fit)
}

predict_platt <- function(obj, score) {
  p <- suppressWarnings(
    predict(
      obj$fit,
      newdata = data.frame(score = score),
      type = "response"
    )
  )
  pmin(pmax(as.numeric(p), 0), 1)
}

fit_isotonic <- function(score, y) {
  if (length(unique(y)) != 2L) return(NULL)

  z <- data.table(
    score = as.numeric(score),
    y = as.numeric(y)
  )[
    ,
    .(
      ybar = mean(y),
      weight = .N
    ),
    by = score
  ][order(score)]

  n <- nrow(z)
  if (!n) return(NULL)

  level <- numeric(n)
  weight <- numeric(n)
  left <- integer(n)
  right <- integer(n)
  m <- 0L

  for (i in seq_len(n)) {
    m <- m + 1L
    level[m] <- z$ybar[i]
    weight[m] <- z$weight[i]
    left[m] <- i
    right[m] <- i

    while (m >= 2L && level[m - 1L] > level[m]) {
      new_weight <- weight[m - 1L] + weight[m]
      level[m - 1L] <- (
        level[m - 1L] * weight[m - 1L] +
          level[m] * weight[m]
      ) / new_weight
      weight[m - 1L] <- new_weight
      right[m - 1L] <- right[m]
      m <- m - 1L
    }
  }

  fitted <- numeric(n)

  for (b in seq_len(m)) {
    fitted[left[b]:right[b]] <- level[b]
  }

  list(
    type = "isotonic",
    x = z$score,
    p = pmin(pmax(fitted, 0), 1)
  )
}

predict_isotonic <- function(obj, score) {
  if (length(obj$x) == 1L) {
    return(rep(obj$p[1L], length(score)))
  }

  approx(
    x = obj$x,
    y = obj$p,
    xout = score,
    method = "linear",
    rule = 2,
    ties = "ordered"
  )$y
}

fit_calibrator <- function(method, score, y) {
  if (method == "platt") return(fit_platt(score, y))
  if (method == "isotonic") return(fit_isotonic(score, y))
  stop("unknown calibrator: ", method, call. = FALSE)
}

predict_calibrator <- function(obj, score) {
  if (is.null(obj)) return(rep(NA_real_, length(score)))
  if (obj$type == "platt") return(predict_platt(obj, score))
  if (obj$type == "isotonic") return(predict_isotonic(obj, score))
  stop("unknown calibrator object", call. = FALSE)
}

choose_thresholds <- function(p, y, npv_target, ppv_target) {
  ok <- is.finite(p) & !is.na(y)
  p <- p[ok]
  y <- y[ok]

  if (!length(p) || length(unique(y)) != 2L) {
    return(list(
      lower = NA_real_,
      upper = NA_real_,
      train_npv = NA_real_,
      train_ppv = NA_real_,
      n_normal_calls = 0L,
      n_malignant_calls = 0L
    ))
  }

  thresholds <- sort(unique(p))

  low <- rbindlist(lapply(thresholds, function(t) {
    idx <- p <= t
    if (!any(idx)) return(NULL)
    data.table(
      threshold = t,
      n_normal_calls = sum(idx),
      npv = mean(y[idx] == 0L)
    )
  }))

  high <- rbindlist(lapply(thresholds, function(t) {
    idx <- p >= t
    if (!any(idx)) return(NULL)
    data.table(
      threshold = t,
      n_malignant_calls = sum(idx),
      ppv = mean(y[idx] == 1L)
    )
  }))

  low <- low[npv >= npv_target]
  high <- high[ppv >= ppv_target]

  if (!nrow(low) || !nrow(high)) {
    return(list(
      lower = NA_real_,
      upper = NA_real_,
      train_npv = NA_real_,
      train_ppv = NA_real_,
      n_normal_calls = 0L,
      n_malignant_calls = 0L
    ))
  }

  candidates <- CJ(
    i = seq_len(nrow(low)),
    j = seq_len(nrow(high))
  )

  candidates <- candidates[
    low$threshold[i] < high$threshold[j]
  ]

  if (!nrow(candidates)) {
    return(list(
      lower = NA_real_,
      upper = NA_real_,
      train_npv = NA_real_,
      train_ppv = NA_real_,
      n_normal_calls = 0L,
      n_malignant_calls = 0L
    ))
  }

  candidates[
    ,
    `:=`(
      lower = low$threshold[i],
      upper = high$threshold[j],
      train_npv = low$npv[i],
      train_ppv = high$ppv[j],
      n_normal_calls = low$n_normal_calls[i],
      n_malignant_calls = high$n_malignant_calls[j]
    )
  ]

  candidates[
    ,
    n_determinate := n_normal_calls + n_malignant_calls
  ]

  setorder(
    candidates,
    -n_determinate,
    lower,
    upper
  )

  best <- candidates[1L]

  list(
    lower = best$lower,
    upper = best$upper,
    train_npv = best$train_npv,
    train_ppv = best$train_ppv,
    n_normal_calls = best$n_normal_calls,
    n_malignant_calls = best$n_malignant_calls
  )
}

make_decision <- function(p, lower, upper) {
  out <- rep("indeterminate", length(p))

  if (is.finite(lower)) {
    out[is.finite(p) & p <= lower] <- "normal"
  }

  if (is.finite(upper)) {
    out[is.finite(p) & p >= upper] <- "malignant"
  }

  out
}

outer_map <- make_group_folds(d, outer_k, seed)
d[, outer_fold := unname(outer_map[gid])]

group_fold_assignment <- unique(
  d[
    ,
    c(
      list(gid = gid, y = y, outer_fold = outer_fold),
      .SD
    ),
    .SDcols = group_cols
  ]
)

methods <- c("platt", "isotonic")
outer_predictions <- vector("list", outer_k)
fold_summaries <- vector("list", outer_k)

for (k in seq_len(outer_k)) {
  tr <- d[outer_fold != k]
  te <- d[outer_fold == k]

  if (!nrow(tr) || !nrow(te)) {
    stop("empty outer fold ", k, call. = FALSE)
  }

  if (length(unique(tr$y)) != 2L || length(unique(te$y)) != 2L) {
    stop("outer fold ", k, " lacks one class", call. = FALSE)
  }

  inner_map <- make_group_folds(
    tr,
    inner_k,
    seed + 1000L * k
  )

  tr[, inner_fold := unname(inner_map[gid])]

  inner_predictions <- vector("list", inner_k)

  for (j in seq_len(inner_k)) {
    itr <- tr[inner_fold != j]
    iva <- tr[inner_fold == j]

    if (!nrow(itr) || !nrow(iva)) {
      stop(
        "empty inner fold ",
        j,
        " in outer fold ",
        k,
        call. = FALSE
      )
    }

    if (
      length(unique(itr$y)) != 2L ||
      length(unique(iva$y)) != 2L
    ) {
      stop(
        "inner fold ",
        j,
        " in outer fold ",
        k,
        " lacks one class",
        call. = FALSE
      )
    }

    inner_dt <- data.table(
      row_id = iva$.row_id,
      y = iva$y
    )

    for (method in methods) {
      obj <- fit_calibrator(
        method,
        itr$burden,
        itr$y
      )

      inner_dt[
        ,
        (paste0("p_", method)) :=
          predict_calibrator(obj, iva$burden)
      ]
    }

    inner_predictions[[j]] <- inner_dt
  }

  inner_oof <- rbindlist(
    inner_predictions,
    use.names = TRUE
  )

  setorder(inner_oof, row_id)

  inner_logloss <- sapply(
    methods,
    function(method) {
      p <- inner_oof[[paste0("p_", method)]]
      if (all(!is.finite(p))) return(Inf)
      logloss(p, inner_oof$y)
    }
  )

  if (all(!is.finite(inner_logloss))) {
    stop(
      "all calibrators failed in outer fold ",
      k,
      call. = FALSE
    )
  }

  best_method <- names(
    inner_logloss == min(inner_logloss, na.rm = TRUE)
  )[1L]

  selected_inner_p <- inner_oof[[paste0("p_", best_method)]]

  thresholds <- choose_thresholds(
    p = selected_inner_p,
    y = inner_oof$y,
    npv_target = target_npv,
    ppv_target = target_ppv
  )

  fit_platt_outer <- fit_platt(
    tr$burden,
    tr$y
  )

  fit_iso_outer <- fit_isotonic(
    tr$burden,
    tr$y
  )

  p_platt_test <- predict_calibrator(
    fit_platt_outer,
    te$burden
  )

  p_iso_test <- predict_calibrator(
    fit_iso_outer,
    te$burden
  )

  selected_fit <- if (best_method == "platt") {
    fit_platt_outer
  } else {
    fit_iso_outer
  }

  p_cal_test <- predict_calibrator(
    selected_fit,
    te$burden
  )

  if (any(!is.finite(p_cal_test))) {
    stop(
      "non-finite selected calibrated probabilities in outer fold ",
      k,
      call. = FALSE
    )
  }

  train_prevalence <- mean(tr$y)
  p_null_test <- rep(
    train_prevalence,
    nrow(te)
  )

  pred <- copy(te)

  pred[, calibrator := best_method]
  pred[, p_null := p_null_test]
  pred[, p_platt := p_platt_test]
  pred[, p_isotonic := p_iso_test]
  pred[, p_cal := p_cal_test]
  pred[, threshold_lower := thresholds$lower]
  pred[, threshold_upper := thresholds$upper]
  pred[
    ,
    decision := make_decision(
      p_cal,
      thresholds$lower,
      thresholds$upper
    )
  ]

  outer_predictions[[k]] <- pred

  fold_summaries[[k]] <- data.table(
    outer_fold = k,
    n_train_rows = nrow(tr),
    n_test_rows = nrow(te),
    n_train_groups = uniqueN(tr$gid),
    n_test_groups = uniqueN(te$gid),
    train_prevalence = train_prevalence,
    selected_calibrator = best_method,
    inner_logloss_platt = unname(inner_logloss["platt"]),
    inner_logloss_isotonic = unname(inner_logloss["isotonic"]),
    threshold_lower = thresholds$lower,
    threshold_upper = thresholds$upper,
    training_NPV_at_lower = thresholds$train_npv,
    training_PPV_at_upper = thresholds$train_ppv,
    n_training_normal_calls = thresholds$n_normal_calls,
    n_training_malignant_calls = thresholds$n_malignant_calls
  )
}

pred <- rbindlist(
  outer_predictions,
  use.names = TRUE,
  fill = TRUE
)

fold_summary <- rbindlist(
  fold_summaries,
  use.names = TRUE
)

setorder(pred, .row_id)

if (nrow(pred) != nrow(d)) {
  stop(
    "outer-CV predictions do not cover every eligible row exactly once",
    call. = FALSE
  )
}

if (anyDuplicated(pred$.row_id)) {
  stop(
    "duplicate outer-CV prediction rows detected",
    call. = FALSE
  )
}

metric_block <- function(z) {
  list(
    n = nrow(z),
    n_groups = uniqueN(z$gid),
    n_tumor = sum(z$y == 1L),
    n_normal = sum(z$y == 0L),
    prevalence = mean(z$y),

    raw_AUROC = auroc(z$burden, z$y),
    raw_AUPRC_AP = average_precision(z$burden, z$y),

    calibrated_AUROC = auroc(z$p_cal, z$y),
    calibrated_AUPRC_AP = average_precision(z$p_cal, z$y),

    Brier_null = brier(z$p_null, z$y),
    Brier_platt = brier(z$p_platt, z$y),
    Brier_isotonic = brier(z$p_isotonic, z$y),
    Brier_cal = brier(z$p_cal, z$y),
    Delta_Brier_cal_minus_null =
      brier(z$p_cal, z$y) -
      brier(z$p_null, z$y),

    logloss_null = logloss(z$p_null, z$y),
    logloss_platt = logloss(z$p_platt, z$y),
    logloss_isotonic = logloss(z$p_isotonic, z$y),
    logloss_cal = logloss(z$p_cal, z$y),
    Delta_logloss_cal_minus_null =
      logloss(z$p_cal, z$y) -
      logloss(z$p_null, z$y),

    ECE_null = ece(z$p_null, z$y, ece_bins),
    ECE_platt = ece(z$p_platt, z$y, ece_bins),
    ECE_isotonic = ece(z$p_isotonic, z$y, ece_bins),
    ECE_cal = ece(z$p_cal, z$y, ece_bins),

    indeterminate_fraction =
      mean(z$decision == "indeterminate"),

    normal_call_fraction =
      mean(z$decision == "normal"),

    malignant_call_fraction =
      mean(z$decision == "malignant"),

    sensitivity_determinate = {
      idx <- z$decision != "indeterminate" & z$y == 1L
      if (!any(idx)) NA_real_
      else mean(z$decision[idx] == "malignant")
    },

    specificity_determinate = {
      idx <- z$decision != "indeterminate" & z$y == 0L
      if (!any(idx)) NA_real_
      else mean(z$decision[idx] == "normal")
    },

    PPV = {
      idx <- z$decision == "malignant"
      if (!any(idx)) NA_real_
      else mean(z$y[idx] == 1L)
    },

    NPV = {
      idx <- z$decision == "normal"
      if (!any(idx)) NA_real_
      else mean(z$y[idx] == 0L)
    }
  )
}

summary_by_fraction <- pred[
  ,
  metric_block(.SD),
  by = tumor_fraction
][order(tumor_fraction)]

summary_overall <- as.data.table(
  metric_block(pred)
)

reliability_bins <- pred[
  ,
  {
    br <- cut(
      p_cal,
      breaks = seq(
        0,
        1,
        length.out = ece_bins + 1L
      ),
      include.lowest = TRUE,
      right = TRUE
    )

    data.table(
      bin = br,
      y = y,
      p_cal = p_cal
    )[
      ,
      .(
        n = .N,
        mean_probability = mean(p_cal),
        observed_tumor_fraction = mean(y)
      ),
      by = bin
    ]
  },
  by = tumor_fraction
][order(tumor_fraction, bin)]

cat(
  "eligible scored rows: ", nrow(d), "\n",
  "unique source-cell groups: ", uniqueN(d$gid), "\n",
  "source groups recurring across >1 context: ",
  sum(group_context$n_contexts > 1L), "\n",
  "group columns: ", paste(group_cols, collapse = ", "), "\n",
  "outer folds: ", outer_k, "\n",
  "inner folds: ", inner_k, "\n",
  "candidate calibrators: Platt, isotonic\n",
  "target NPV: ", target_npv, "\n",
  "target PPV: ", target_ppv, "\n\n",
  sep = ""
)

print(fold_summary)
cat("\n")
print(summary_by_fraction)
cat("\n")
print(summary_overall)

fwrite(
  pred,
  file.path(
    outdir,
    "copykat_calibration_oof_percell.tsv"
  ),
  sep = "\t"
)

fwrite(
  fold_summary,
  file.path(
    outdir,
    "copykat_calibration_outerfolds.tsv"
  ),
  sep = "\t"
)

fwrite(
  group_fold_assignment,
  file.path(
    outdir,
    "copykat_calibration_group_folds.tsv"
  ),
  sep = "\t"
)

fwrite(
  summary_by_fraction,
  file.path(
    outdir,
    "copykat_calibration_summary_by_fraction.tsv"
  ),
  sep = "\t"
)

fwrite(
  summary_overall,
  file.path(
    outdir,
    "copykat_calibration_summary_overall.tsv"
  ),
  sep = "\t"
)

fwrite(
  reliability_bins,
  file.path(
    outdir,
    "copykat_calibration_reliability_bins.tsv"
  ),
  sep = "\t"
)

cat(
  "\nwrote:\n",
  file.path(outdir, "copykat_calibration_oof_percell.tsv"), "\n",
  file.path(outdir, "copykat_calibration_outerfolds.tsv"), "\n",
  file.path(outdir, "copykat_calibration_group_folds.tsv"), "\n",
  file.path(outdir, "copykat_calibration_summary_by_fraction.tsv"), "\n",
  file.path(outdir, "copykat_calibration_summary_overall.tsv"), "\n",
  file.path(outdir, "copykat_calibration_reliability_bins.tsv"), "\n",
  sep = ""
)
