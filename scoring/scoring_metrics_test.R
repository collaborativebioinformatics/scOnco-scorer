#!/usr/bin/env Rscript
# =============================================================================
# Stage 3b - Validate the metrics on synthetic callers with known properties.
#
# Run this before 03_metrics.R ever sees real caller output. If the
# overconfident case doesn't produce bad ECE with good AUROC, your calibration
# binning is wrong - and you would never catch that on real data.
#
# Usage:  Rscript 04_validate_metrics.R      (exits non-zero on failure)
# =============================================================================

source("03_metrics.R")

set.seed(20260826)
N <- 5000
pass <- TRUE

check <- function(label, ok, detail = "") {
  cat(sprintf("[%s] %-42s %s\n", if (ok) "PASS" else "FAIL", label, detail))
  if (!ok) pass <<- FALSE
}

# --- 1. perfect caller --------------------------------------------------------
y <- rbinom(N, 1, 0.3)
p <- ifelse(y == 1, 0.999, 0.001)
r <- score_calls(y, p)
check("perfect: AUROC ~ 1",   abs(r$metrics$auroc - 1) < 1e-9, sprintf("auroc=%.4f", r$metrics$auroc))
check("perfect: Brier ~ 0",   r$metrics$brier < 1e-3,          sprintf("brier=%.5f", r$metrics$brier))
check("perfect: ECE ~ 0",     r$metrics$ece   < 1e-2,          sprintf("ece=%.5f",   r$metrics$ece))
check("perfect: sens/spec 1", r$metrics$sensitivity == 1 && r$metrics$specificity == 1)

# --- 2. random caller ---------------------------------------------------------
y <- rbinom(N, 1, 0.3)
p <- runif(N)
r <- score_calls(y, p)
check("random: AUROC ~ 0.5",  abs(r$metrics$auroc - 0.5) < 0.05, sprintf("auroc=%.4f", r$metrics$auroc))
check("random: AUPRC ~ prev", abs(r$metrics$auprc - mean(y)) < 0.05,
      sprintf("auprc=%.3f prev=%.3f", r$metrics$auprc, mean(y)))

# --- 3. well-calibrated vs overconfident, SAME ranking ------------------------
# This is the whole point of reporting calibration alongside discrimination:
# both callers rank cells identically, so AUROC is identical - only ECE and
# Brier can tell them apart.
latent  <- rnorm(N)
p_true  <- plogis(latent)
y       <- rbinom(N, 1, p_true)

r_cal  <- score_calls(y, p_true)                # honest probabilities
r_over <- score_calls(y, plogis(3 * latent))    # same order, pushed to extremes

check("calibrated: ECE small", r_cal$metrics$ece < 0.05, sprintf("ece=%.4f", r_cal$metrics$ece))
check("overconfident: ECE large", r_over$metrics$ece > 0.10, sprintf("ece=%.4f", r_over$metrics$ece))
# With latent ~ N(0,1) and y ~ Bernoulli(plogis(latent)), the *theoretical* ceiling
# on AUROC is about 0.76 - the labels are genuinely noisy. The point of this check
# is that discrimination stays well above chance while calibration falls apart.
check("overconfident: AUROC still high", r_over$metrics$auroc > 0.70,
      sprintf("auroc=%.4f (chance = 0.5, ceiling for this generator ~0.76)",
              r_over$metrics$auroc))
check("AUROC identical (ranking preserved)",
      abs(r_cal$metrics$auroc - r_over$metrics$auroc) < 1e-9)
check("Brier worse when overconfident", r_over$metrics$brier > r_cal$metrics$brier,
      sprintf("%.4f vs %.4f", r_over$metrics$brier, r_cal$metrics$brier))

# --- 4. low-prevalence behaviour (this is your 1-5% purity regime) ------------
y <- rbinom(N, 1, 0.02)
p <- rep(0.0, N)                                 # "everything is normal"
r <- score_calls(y, p)
check("all-normal caller: sens = 0", r$metrics$sensitivity == 0)
check("all-normal caller: spec = 1", r$metrics$specificity == 1)
check("all-normal caller: AUROC = 0.5", abs(r$metrics$auroc - 0.5) < 1e-9,
      "a useless caller still looks 98% accurate at 2% purity - never report accuracy")

# --- 5. edge cases ------------------------------------------------------------
check("single-class truth returns NA, not error",
      is.na(auroc(rep(1L, 100), runif(100))))
check("out-of-range probs rejected",
      inherits(try(score_calls(c(0L,1L), c(0.5, 1.5)), silent = TRUE), "try-error"))

cat("\n", if (pass) "ALL CHECKS PASSED - metrics module is safe to use.\n"
         else "SOME CHECKS FAILED - fix before scoring real output.\n", sep = "")
quit(status = if (pass) 0 else 1)
