library(ggplot2)

source("calibration_module.R")

# ============================================================
# calibration_simulation.R
#
# Validates calibration_module.R on SIMULATED labelled scores,
# per the assignment, before it is trusted on any real caller
# output. Demonstrates three specific things:
#
#   1. Ranking (AUC) can stay good while calibration drifts --
#      a monotone-but-unscaled score has high discrimination and
#      poor calibration; Platt/isotonic fix the calibration while
#      leaving discrimination essentially unchanged.
#
#   2. Row-level CV splitting leaks when observations repeat
#      (clustered/duplicated units); grouped splitting does not.
#
#   3. A raw binary label is not a calibrated probability -- using
#      one realized 0/1 outcome to predict another independent
#      draw from the same underlying probability is provably worse
#      (in expectation) than using the true probability itself.
#
# Outputs:
#   - simulation_reliability_plot.pdf/png  (scenario 1)
#   - simulation_summary_tables.tsv        (all three scenarios)
# ============================================================

out_dir <- "calibration_sim_outputs"
dir.create(out_dir, showWarnings = FALSE)

theme_sim <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

# ============================================================
# Scenario 1: ranking stays good while calibration drifts
#
# True probability p_true drives the label. The "score" available
# to the calibrator is a monotone transform of p_true (so ranking/
# AUC is preserved) but is NOT on a 0-1 probability scale (poorly
# scaled logit + noise) -- so treating the raw score as a
# probability is badly calibrated even though it ranks well.
# ============================================================

set.seed(101)
n1 <- 4000
group1 <- paste0("cell_", seq_len(n1))
p_true1 <- runif(n1, 0.02, 0.98)
label1 <- rbinom(n1, 1, p_true1)

# Poorly-scaled monotone score: same ranking as p_true, wrong scale
score1 <- qlogis(p_true1) * 2.5 + 1.5

sim1 <- data.frame(group_id = group1, score = score1, label = label1)

# Discrimination (AUC) of the raw score vs. the true probability --
# should be nearly identical, since the transform is monotone.
auc_raw <- auc_score(score1, label1)
auc_true_p <- auc_score(p_true1, label1)

# Calibration of the RAW score if naively treated as a probability
# (min-max squashed to [0,1] just to make Brier/ECE computable at
# all -- illustrating that even a "probability-shaped" naive use
# of a monotone score is poorly calibrated).
naive_prob1 <- (score1 - min(score1)) / (max(score1) - min(score1))
naive_brier <- brier_score(naive_prob1, label1)
naive_ece <- expected_calibration_error(naive_prob1, label1)

cv1 <- run_grouped_calibration_cv(
  sim1, "group_id", "score", "label",
  k = 5, methods = c("platt", "isotonic"), seed = 11
)

scenario1_summary <- rbind(
  data.frame(
    method = "raw_score_naively_scaled_as_prob",
    mean_auc = auc_raw,
    mean_brier = naive_brier,
    mean_log_loss = NA_real_,
    mean_ece = naive_ece
  ),
  cv1$summary
)

cat("=== Scenario 1: ranking vs. calibration ===\n")
print(scenario1_summary)
cat(
  "\nAUC of raw score:", round(auc_raw, 3),
  " | AUC of true p:", round(auc_true_p, 3),
  " (should be nearly identical -- monotone transform preserves ranking)\n"
)
cat(
  "Naive-scaled Brier:", round(naive_brier, 4),
  " vs. calibrated Brier (mean across methods):",
  round(mean(cv1$summary$mean_brier), 4),
  " -- calibration should improve substantially despite AUC being unchanged.\n\n"
)

# Reliability plot: naive (uncalibrated) vs. isotonic-calibrated,
# pooled across CV folds.
naive_rel <- reliability_curve_data(naive_prob1, label1)
naive_rel$method <- "raw_score_naive"

iso_rel <- cv1$reliability[cv1$reliability$method == "isotonic", ]
iso_rel_pooled <- do.call(rbind, lapply(
  split(iso_rel, cut(iso_rel$mean_predicted, breaks = seq(0, 1, 0.1), include.lowest = TRUE)),
  function(sub) {
    if (nrow(sub) == 0) return(NULL)
    data.frame(
      mean_predicted = weighted.mean(sub$mean_predicted, sub$n),
      observed_frequency = weighted.mean(sub$observed_frequency, sub$n),
      n = sum(sub$n)
    )
  }
))
iso_rel_pooled$method <- "isotonic_calibrated"

rel_plot_data <- rbind(
  naive_rel[, c("mean_predicted", "observed_frequency", "n", "method")],
  iso_rel_pooled[, c("mean_predicted", "observed_frequency", "n", "method")]
)

fig_reliability <- ggplot(
  rel_plot_data,
  aes(x = mean_predicted, y = observed_frequency, color = method)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 0.9) +
  geom_point(aes(size = n)) +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(
    values = c(raw_score_naive = "#c0392b", isotonic_calibrated = "#1b6ca8"),
    labels = c(raw_score_naive = "Raw score (naive)", isotonic_calibrated = "Isotonic-calibrated")
  ) +
  labs(
    x = "Mean predicted probability",
    y = "Observed frequency",
    color = NULL, size = "n in bin",
    title = "Simulation: calibration fixes reliability without changing ranking",
    subtitle = paste0(
      "AUC unchanged (", round(auc_raw, 3), " raw vs. ", round(auc_true_p, 3),
      " true p); calibration curve moves onto the diagonal"
    )
  ) +
  theme_sim +
  theme(legend.position = "top")

ggsave(file.path(out_dir, "simulation_reliability_plot.pdf"), fig_reliability, width = 6.5, height = 5.5)
ggsave(file.path(out_dir, "simulation_reliability_plot.png"), fig_reliability, width = 6.5, height = 5.5, dpi = 300)

# ============================================================
# Scenario 2: row-level splitting leaks repeated observations;
# grouped splitting does not.
#
# Simulate clustered data -- each "cluster" (e.g. a source cell
# reused across mixtures/replicates) has a fixed difficulty, and
# appears multiple times with near-duplicate scores/labels. Row-
# level random CV lets near-duplicates of the same cluster land in
# both train and test, inflating apparent performance. Grouped CV
# (by cluster) does not.
# ============================================================

set.seed(202)
n_clusters <- 150
reps_per_cluster <- 6

cluster_id <- rep(paste0("clu_", 1:n_clusters), each = reps_per_cluster)
cluster_p <- rep(runif(n_clusters, 0.05, 0.95), each = reps_per_cluster)

# Near-duplicate rows within a cluster: same underlying p, label
# redrawn each rep, small measurement noise on the score.
score2 <- qlogis(cluster_p) + rnorm(length(cluster_id), 0, 0.15)
label2 <- rbinom(length(cluster_id), 1, cluster_p)

sim2 <- data.frame(
  cluster_id = cluster_id,
  row_id = paste0("row_", seq_along(cluster_id)),
  score = score2,
  label = label2
)

n_repeats <- 8
row_level_briers <- numeric(n_repeats)
grouped_briers <- numeric(n_repeats)

for (rep_i in seq_len(n_repeats)) {
  
  cv_row <- run_grouped_calibration_cv(
    sim2, "row_id", "score", "label",
    k = 5, methods = "isotonic", seed = 1000 + rep_i
  )
  row_level_briers[rep_i] <- mean(cv_row$per_fold$brier)
  
  cv_grouped <- run_grouped_calibration_cv(
    sim2, "cluster_id", "score", "label",
    k = 5, methods = "isotonic", seed = 1000 + rep_i
  )
  grouped_briers[rep_i] <- mean(cv_grouped$per_fold$brier)
}

scenario2_summary <- data.frame(
  split_type = c("row_level_leaky", "grouped_by_cluster"),
  mean_brier_across_repeats = c(mean(row_level_briers), mean(grouped_briers)),
  sd_brier_across_repeats = c(sd(row_level_briers), sd(grouped_briers))
)

cat("=== Scenario 2: row-level leakage vs. grouped CV ===\n")
print(scenario2_summary)
cat(
  "\nRow-level splitting should show LOWER (optimistically biased) Brier",
  "than grouped splitting, because near-duplicate cluster members leak",
  "across the train/test boundary.\n\n"
)

# ============================================================
# Scenario 3: a binary label is not a calibrated probability.
#
# For a fixed true probability p, using one realized 0/1 draw to
# predict an independent second draw from the same p has expected
# Brier score 2*p*(1-p) -- exactly double the Bayes-optimal Brier
# of p*(1-p) achieved by using the true probability itself. This
# holds deterministically in expectation and is easy to show by
# simulation across a range of p.
# ============================================================

set.seed(303)
p_grid <- seq(0.05, 0.95, by = 0.05)
n_trials <- 20000

scenario3_rows <- list()
for (i in seq_along(p_grid)) {
  p <- p_grid[i]
  draw_a <- rbinom(n_trials, 1, p)  # "using a binary label as the prediction"
  draw_b <- rbinom(n_trials, 1, p)  # the outcome being predicted
  
  brier_label_as_pred <- mean((draw_a - draw_b)^2)
  brier_true_p_as_pred <- mean((p - draw_b)^2)
  
  scenario3_rows[[i]] <- data.frame(
    p = p,
    brier_binary_label_as_prediction = brier_label_as_pred,
    brier_true_probability_as_prediction = brier_true_p_as_pred,
    theoretical_2p1minusp = 2 * p * (1 - p),
    theoretical_p1minusp = p * (1 - p)
  )
}

scenario3_summary <- do.call(rbind, scenario3_rows)

cat("=== Scenario 3: binary label vs. true probability as prediction ===\n")
print(scenario3_summary)
cat(
  "\nBrier using a binary label as the prediction should track",
  "2*p*(1-p); Brier using the true probability should track p*(1-p) --",
  "exactly half the error, at every p. This is why CopyKAT's binary",
  "call / SCEVAN's categorical label cannot be calibrated as-is.\n\n"
)

# ============================================================
# Write all summary tables
# ============================================================

write.table(
  cv1$per_fold, file.path(out_dir, "scenario1_per_fold.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  scenario1_summary, file.path(out_dir, "scenario1_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  scenario2_summary, file.path(out_dir, "scenario2_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  scenario3_summary, file.path(out_dir, "scenario3_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("Done. Simulation outputs written to", out_dir, "\n")