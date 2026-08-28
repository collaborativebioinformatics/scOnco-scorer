# ============================================================
# calibration_module.R
#
# Generic, reusable probability calibration toolkit. Not tied to
# any single caller -- applies to any (group_id, score, label)
# table. Intended application order once real continuous scores
# exist: Numbat joint posterior -> expression-only -> allele-only
# posterior -> CopyKAT continuous burden (once available from the
# copy-number matrix) -> inferCNV/SCEVAN continuous scores, if
# they materialize.
#
# Explicitly NOT for CopyKAT's binary aneuploid/diploid call or
# SCEVAN's categorical label -- those are not probabilities and
# calibrating them as such is out of scope by design (see
# scenario 3 in calibration_simulation.R for why).
#
# All fitting happens on a TRAINING fold only; every calibrator,
# threshold, and metric returned by the CV driver is evaluated on
# a held-out fold the calibrator never saw.
# ============================================================

# ------------------------------------------------------------
# Grouped K-fold assignment
#
# Splits by GROUP, not by row -- every row belonging to the same
# group_id lands in exactly one fold. Prevents leakage when the
# same underlying unit (e.g. a source cell reused across mixtures
# or replicates) appears more than once in the data.
# ------------------------------------------------------------

grouped_kfold <- function(groups, k = 5, seed = 1) {
  set.seed(seed)
  unique_groups <- unique(groups)
  n_groups <- length(unique_groups)
  
  if (n_groups < k) {
    stop("Fewer unique groups (", n_groups, ") than folds (", k, ").")
  }
  
  fold_assignment <- sample(rep(seq_len(k), length.out = n_groups))
  names(fold_assignment) <- unique_groups
  
  fold_assignment[as.character(groups)]
}

# ------------------------------------------------------------
# Calibrators. Each fit_* function is fit on TRAINING data only
# and returns a predict function closing over the fitted model --
# the caller never has access to test-fold data at fit time.
# ------------------------------------------------------------

fit_platt <- function(train_scores, train_labels) {
  fit <- glm(
    train_labels ~ train_scores,
    family = binomial(link = "logit")
  )
  function(new_scores) {
    predict(
      fit,
      newdata = data.frame(train_scores = new_scores),
      type = "response"
    )
  }
}

fit_isotonic <- function(train_scores, train_labels) {
  ord <- order(train_scores)
  iso <- isoreg(train_scores[ord], train_labels[ord])
  # isoreg has no predict() method -- build a monotone step
  # function from its fitted (x, yf) pairs. rule=2 clamps
  # extrapolation to the boundary fitted values.
  step_fun <- approxfun(
    x = iso$x, y = iso$yf, method = "constant", rule = 2, ties = mean
  )
  function(new_scores) {
    pmin(pmax(step_fun(new_scores), 0), 1)
  }
}

# Beta calibration only if the package is actually available --
# per instruction, no manual reimplementation as a substitute.
fit_beta_calibration <- function(train_scores, train_labels) {
  if (!requireNamespace("betacal", quietly = TRUE)) {
    message("betacal package not available -- skipping beta calibration.")
    return(NULL)
  }
  fit <- betacal::beta_calibration(train_scores, train_labels, parameters = "abm")
  function(new_scores) {
    betacal::beta_predict(new_scores, fit)
  }
}

# ------------------------------------------------------------
# Metrics. All take (probs, labels) with labels in {0, 1}.
# ------------------------------------------------------------

brier_score <- function(probs, labels) {
  mean((probs - labels)^2)
}

log_loss <- function(probs, labels, eps = 1e-15) {
  probs <- pmin(pmax(probs, eps), 1 - eps)
  -mean(labels * log(probs) + (1 - labels) * log(1 - probs))
}

expected_calibration_error <- function(probs, labels, n_bins = 10) {
  bins <- cut(
    probs, breaks = seq(0, 1, length.out = n_bins + 1), include.lowest = TRUE
  )
  n <- length(probs)
  
  bin_stats <- do.call(rbind, lapply(
    split(data.frame(probs, labels), bins),
    function(sub) {
      if (nrow(sub) == 0) return(NULL)
      data.frame(
        n = nrow(sub),
        mean_pred = mean(sub$probs),
        mean_obs = mean(sub$labels)
      )
    }
  ))
  
  sum(bin_stats$n / n * abs(bin_stats$mean_pred - bin_stats$mean_obs))
}

reliability_curve_data <- function(probs, labels, n_bins = 10) {
  bins <- cut(
    probs, breaks = seq(0, 1, length.out = n_bins + 1), include.lowest = TRUE
  )
  
  do.call(rbind, lapply(
    split(data.frame(probs, labels), bins),
    function(sub) {
      if (nrow(sub) == 0) return(NULL)
      data.frame(
        mean_predicted = mean(sub$probs),
        observed_frequency = mean(sub$labels),
        n = nrow(sub)
      )
    }
  ))
}

# Rank-based AUC (Mann-Whitney U) -- measures discrimination
# (ranking quality), independent of calibration. A score can have
# high AUC and still be badly calibrated; that distinction is the
# whole point of scenario 1 in calibration_simulation.R.
auc_score <- function(scores, labels) {
  n1 <- sum(labels == 1)
  n0 <- sum(labels == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(scores)
  (sum(r[labels == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# ------------------------------------------------------------
# Threshold selection -- TRAINING DATA ONLY. The CV driver below
# calls this on the training fold and applies the resulting
# threshold to the held-out fold; it is never refit on test data.
# ------------------------------------------------------------

select_threshold_from_training <- function(
    train_probs, train_labels, metric = c("youden", "f1")
) {
  metric <- match.arg(metric)
  candidates <- sort(unique(train_probs))
  
  score_at <- function(t) {
    pred <- as.integer(train_probs >= t)
    tp <- sum(pred == 1 & train_labels == 1)
    fp <- sum(pred == 1 & train_labels == 0)
    fn <- sum(pred == 0 & train_labels == 1)
    tn <- sum(pred == 0 & train_labels == 0)
    
    sens <- if ((tp + fn) > 0) tp / (tp + fn) else 0
    spec <- if ((tn + fp) > 0) tn / (tn + fp) else 0
    prec <- if ((tp + fp) > 0) tp / (tp + fp) else 0
    
    if (metric == "youden") {
      sens + spec - 1
    } else {
      if ((prec + sens) > 0) 2 * prec * sens / (prec + sens) else 0
    }
  }
  
  scores <- vapply(candidates, score_at, numeric(1))
  candidates[which.max(scores)]
}

# ------------------------------------------------------------
# Grouped outer CV driver.
#
# data: data.frame with group_col, score_col, label_col
# methods: any of "platt", "isotonic", "beta"
#
# For each fold: fit every requested calibrator + a training-only
# threshold using ONLY the training rows, then evaluate every
# metric on the held-out rows using that fold's fitted calibrator
# and threshold. Never refits or re-selects anything on test data.
# ------------------------------------------------------------

run_grouped_calibration_cv <- function(
    data, group_col, score_col, label_col,
    k = 5, methods = c("platt", "isotonic"), seed = 1
) {
  groups <- data[[group_col]]
  scores <- data[[score_col]]
  labels <- data[[label_col]]
  
  stopifnot(all(labels %in% c(0, 1)))
  
  folds <- grouped_kfold(groups, k = k, seed = seed)
  
  fitters <- list(
    platt = fit_platt,
    isotonic = fit_isotonic,
    beta = fit_beta_calibration
  )
  
  results <- list()
  reliability_rows <- list()
  result_i <- 1
  
  for (method in methods) {
    if (!method %in% names(fitters)) {
      stop("Unknown calibration method: ", method)
    }
    
    for (fold_id in sort(unique(folds))) {
      
      train_idx <- which(folds != fold_id)
      test_idx <- which(folds == fold_id)
      
      train_scores <- scores[train_idx]
      train_labels <- labels[train_idx]
      test_scores <- scores[test_idx]
      test_labels <- labels[test_idx]
      
      calibrator <- fitters[[method]](train_scores, train_labels)
      
      if (is.null(calibrator)) {
        # e.g. beta calibration unavailable this fold/environment
        next
      }
      
      test_probs <- calibrator(test_scores)
      
      train_probs_for_threshold <- calibrator(train_scores)
      threshold <- select_threshold_from_training(
        train_probs_for_threshold, train_labels
      )
      
      raw_auc <- auc_score(test_scores, test_labels)
      cal_brier <- brier_score(test_probs, test_labels)
      cal_logloss <- log_loss(test_probs, test_labels)
      cal_ece <- expected_calibration_error(test_probs, test_labels)
      
      results[[result_i]] <- data.frame(
        method = method,
        fold = fold_id,
        n_train = length(train_idx),
        n_test = length(test_idx),
        auc = raw_auc,
        brier = cal_brier,
        log_loss = cal_logloss,
        ece = cal_ece,
        threshold_from_training = threshold
      )
      
      rel <- reliability_curve_data(test_probs, test_labels)
      rel$method <- method
      rel$fold <- fold_id
      reliability_rows[[result_i]] <- rel
      
      result_i <- result_i + 1
    }
  }
  
  per_fold <- do.call(rbind, results)
  reliability <- do.call(rbind, reliability_rows)
  
  summary_table <- do.call(rbind, lapply(
    split(per_fold, per_fold$method),
    function(sub) {
      data.frame(
        method = sub$method[1],
        mean_auc = mean(sub$auc),
        mean_brier = mean(sub$brier),
        mean_log_loss = mean(sub$log_loss),
        mean_ece = mean(sub$ece)
      )
    }
  ))
  
  list(
    per_fold = per_fold,
    reliability = reliability,
    summary = summary_table
  )
}