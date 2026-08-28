# ============================================================
# calibration_tests.R
#
# Lightweight stopifnot-based sanity tests for calibration_module.R
# (no testthat dependency, consistent with the pragmatic style
# used elsewhere in this project). Source calibration_module.R
# before running this.
#
# Usage: Rscript calibration_tests.R
# ============================================================

source("calibration_module.R")

cat("Running calibration_module.R tests...\n\n")

test_count <- 0
pass_count <- 0

check <- function(desc, expr) {
  test_count <<- test_count + 1
  ok <- tryCatch(isTRUE(expr), error = function(e) FALSE)
  if (ok) {
    pass_count <<- pass_count + 1
    cat("  PASS:", desc, "\n")
  } else {
    cat("  FAIL:", desc, "\n")
  }
}

# ------------------------------------------------------------
# grouped_kfold: no group ever spans more than one fold
# ------------------------------------------------------------

set.seed(42)
groups <- rep(paste0("g", 1:20), each = 5)  # 20 groups, 5 rows each
folds <- grouped_kfold(groups, k = 4, seed = 1)

check(
  "grouped_kfold never splits a group across folds",
  all(tapply(folds, groups, function(f) length(unique(f)) == 1))
)

check(
  "grouped_kfold produces the requested number of folds",
  length(unique(folds)) == 4
)

check(
  "grouped_kfold errors when k exceeds unique group count",
  tryCatch({
    grouped_kfold(rep("only_one_group", 10), k = 5)
    FALSE
  }, error = function(e) TRUE)
)

# ------------------------------------------------------------
# fit_isotonic: predictions are monotone non-decreasing in score
# ------------------------------------------------------------

set.seed(1)
x <- sort(runif(200))
p_true <- x  # perfectly monotone relationship
y <- rbinom(200, 1, p_true)

iso_predict <- fit_isotonic(x, y)
test_x <- seq(0, 1, length.out = 50)
test_y <- iso_predict(test_x)

check(
  "fit_isotonic output is monotone non-decreasing",
  all(diff(test_y) >= -1e-9)
)

check(
  "fit_isotonic output stays within [0, 1]",
  all(test_y >= 0 & test_y <= 1)
)

# ------------------------------------------------------------
# fit_platt: predictions are monotone in score (logistic is
# monotone by construction, assuming a positive coefficient here)
# ------------------------------------------------------------

platt_predict <- fit_platt(x, y)
platt_y <- platt_predict(test_x)

check(
  "fit_platt output is monotone non-decreasing for a positive relationship",
  all(diff(platt_y) >= -1e-9)
)

# ------------------------------------------------------------
# Metric functions against hand-computed toy examples
# ------------------------------------------------------------

toy_probs <- c(0.9, 0.1, 0.5, 0.5)
toy_labels <- c(1, 0, 1, 0)

expected_brier <- mean((toy_probs - toy_labels)^2)
check(
  "brier_score matches hand computation",
  isTRUE(all.equal(brier_score(toy_probs, toy_labels), expected_brier))
)

expected_ll <- -mean(
  toy_labels * log(toy_probs) + (1 - toy_labels) * log(1 - toy_probs)
)
check(
  "log_loss matches hand computation",
  isTRUE(all.equal(log_loss(toy_probs, toy_labels), expected_ll))
)

check(
  "log_loss on perfectly correct extreme probabilities is near zero",
  log_loss(c(0.999999, 0.000001), c(1, 0)) < 0.01
)

check(
  "brier_score is zero for perfect predictions",
  brier_score(c(1, 0, 1, 0), c(1, 0, 1, 0)) == 0
)

check(
  "brier_score is 0.25 for a perfectly uninformative p=0.5 predictor",
  isTRUE(all.equal(brier_score(rep(0.5, 100), rbinom(100, 1, 0.5) * 0 + c(1, 0)), 0.25))
)

# ------------------------------------------------------------
# auc_score: perfect separation gives AUC = 1, random gives ~0.5
# ------------------------------------------------------------

perfect_scores <- c(1, 2, 3, 4, 5, 6)
perfect_labels <- c(0, 0, 0, 1, 1, 1)
check(
  "auc_score is 1 for perfectly separated scores",
  auc_score(perfect_scores, perfect_labels) == 1
)

reversed_labels <- c(1, 1, 1, 0, 0, 0)
check(
  "auc_score is 0 for perfectly inverted scores",
  auc_score(perfect_scores, reversed_labels) == 0
)

# ------------------------------------------------------------
# select_threshold_from_training: only ever touches its inputs,
# never any external test-fold object -- verified by signature
# (function takes only train_probs/train_labels) plus a sanity
# check that it returns a value seen in the training scores.
# ------------------------------------------------------------

set.seed(2)
train_probs <- runif(100)
train_labels <- rbinom(100, 1, train_probs)
thresh <- select_threshold_from_training(train_probs, train_labels)

check(
  "select_threshold_from_training returns a threshold within observed range",
  thresh >= min(train_probs) && thresh <= max(train_probs)
)

# ------------------------------------------------------------
# run_grouped_calibration_cv: end-to-end smoke test
# ------------------------------------------------------------

set.seed(3)
n_groups <- 60
sim_data <- data.frame(
  group_id = rep(paste0("grp", 1:n_groups), each = 5)
)
sim_p <- rep(runif(n_groups), each = 5)
sim_data$score <- qlogis(sim_p) + rnorm(nrow(sim_data), 0, 0.3)
sim_data$label <- rbinom(nrow(sim_data), 1, sim_p)

cv_out <- run_grouped_calibration_cv(
  sim_data, "group_id", "score", "label",
  k = 5, methods = c("platt", "isotonic"), seed = 7
)

check(
  "run_grouped_calibration_cv returns per_fold, reliability, and summary",
  all(c("per_fold", "reliability", "summary") %in% names(cv_out))
)

check(
  "run_grouped_calibration_cv summary has one row per method",
  nrow(cv_out$summary) == 2
)

check(
  "run_grouped_calibration_cv per_fold Brier scores are all in [0, 1]",
  all(cv_out$per_fold$brier >= 0 & cv_out$per_fold$brier <= 1)
)

cat("\n", pass_count, "/", test_count, "tests passed.\n")

if (pass_count != test_count) {
  quit(save = "no", status = 1)
}