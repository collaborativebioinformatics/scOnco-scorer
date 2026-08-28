#!/usr/bin/env Rscript
# =============================================================================
# Validation for scoring_bootstrap.R, on simulated data with known structure.
# Exits non-zero on failure.  Usage: Rscript scoring_bootstrap_test.R
# =============================================================================
source("scoring_bootstrap.R")
pass <- TRUE
check <- function(label, ok, detail = "") {
  cat(sprintf("[%s] %-52s %s\n", if (ok) "PASS" else "FAIL", label, detail))
  if (!ok) pass <<- FALSE
}

# --- simulator: METHODS.md geometry, 6 fractions x N replicates --------------
FRACTIONS <- c(0.01, 0.05, 0.10, 0.20, 0.40, 0.80)
make_sweep <- function(n_pool_A = 400, n_pool_B = 1200, n_rep = 10,
                       effect = 1.4, seed = 1) {
  set.seed(seed)
  # each source cell has a fixed latent quality: this is the cluster effect
  latent <- c(rnorm(n_pool_A, effect, 1), rnorm(n_pool_B, 0, 1))
  ids    <- c(paste0("A", seq_len(n_pool_A)), paste0("B", seq_len(n_pool_B)))
  ytrue  <- c(rep(1L, n_pool_A), rep(0L, n_pool_B))
  rows <- list()
  for (f in FRACTIONS) for (r in seq_len(n_rep)) {
    nA <- max(2, round(3000 * f / 10)); nB <- max(2, round(3000 * (1 - f) / 10))
    iA <- sample(which(ytrue == 1), min(nA, n_pool_A))
    iB <- sample(which(ytrue == 0), min(nB, n_pool_B))
    i  <- c(iA, iB)
    rows[[length(rows) + 1]] <- data.frame(
      source_cell = ids[i], y = ytrue[i],
      p = plogis(latent[i] + rnorm(length(i), 0, 0.35)),  # small per-measurement noise
      fraction = f, replicate = r, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

df <- make_sweep()
cat(sprintf("\nsimulated: %d rows, %d source cells, %.1f rows/cell\n\n",
            nrow(df), length(unique(df$source_cell)),
            nrow(df) / length(unique(df$source_cell))))

# --- 1. the central claim: cluster CIs are WIDER than row CIs ----------------
cb <- cluster_bootstrap(df, n_boot = 400, seed = 7, verbose = FALSE)
rb <- row_bootstrap(df,     n_boot = 400, seed = 7)
w_cl <- cb$upper[["auroc"]] - cb$lower[["auroc"]]
w_rw <- rb$upper[["auroc"]] - rb$lower[["auroc"]]
check("cluster CI wider than row CI (auroc)", w_cl > w_rw,
      sprintf("cluster %.4f vs row %.4f (%.1fx)", w_cl, w_rw, w_cl / w_rw))
check("cluster CI wider than row CI (brier)",
      (cb$upper[["brier"]] - cb$lower[["brier"]]) > (rb$upper[["brier"]] - rb$lower[["brier"]]))

# --- 2. point estimate unchanged by resampling scheme ------------------------
check("point estimates identical across schemes",
      isTRUE(all.equal(cb$point[["auroc"]], rb$point[["auroc"]])))

# --- 3. CI covers the point estimate ----------------------------------------
inside <- cb$point >= cb$lower & cb$point <= cb$upper
check("all CIs contain their point estimate", all(inside, na.rm = TRUE))

# --- 4. stratification holds prevalence fixed -------------------------------
cs <- cluster_bootstrap(df, n_boot = 200, seed = 3, stratify = TRUE,  verbose = FALSE)
cu <- cluster_bootstrap(df, n_boot = 200, seed = 3, stratify = FALSE, verbose = FALSE)
sd_s <- sd(cs$replicates[, "prevalence"]); sd_u <- sd(cu$replicates[, "prevalence"])
check("stratified prevalence varies less", sd_s < sd_u,
      sprintf("sd %.5f vs %.5f", sd_s, sd_u))

# --- 5. determinism ----------------------------------------------------------
a <- cluster_bootstrap(df, n_boot = 100, seed = 42, verbose = FALSE)
b <- cluster_bootstrap(df, n_boot = 100, seed = 42, verbose = FALSE)
check("same seed reproduces the interval", isTRUE(all.equal(a$lower, b$lower)))
check("seed recorded for the ledger", a$seed == 42)

# --- 6. label-shuffle control collapses to chance ---------------------------
ls <- label_shuffle_control(df, n_boot = 300, seed = 5)
check("label shuffle: AUROC CI covers 0.5",
      ls$lower[["auroc"]] <= 0.5 && ls$upper[["auroc"]] >= 0.5,
      sprintf("[%.3f, %.3f]", ls$lower[["auroc"]], ls$upper[["auroc"]]))

# --- 7. per-fraction runs and detectability floor ---------------------------
summ <- boot_by_fraction(df, n_boot = 150, seed = 11)
check("one row per metric per fraction",
      nrow(summ) == length(FRACTIONS) * length(cb$point))
check("every metric has a 95% interval",
      all(!is.na(summ$lower[summ$metric == "auroc"])))
fl <- detectability_floor(summ, "auroc", threshold = 0.80)
check("detectability floor returns a tested fraction or NA",
      is.na(fl) || fl %in% FRACTIONS, sprintf("floor = %s", fl))
lo <- summ$lower[summ$metric == "auroc" & summ$fraction == fl]
check("floor uses the LOWER bound, not the estimate",
      is.na(fl) || lo >= 0.80, sprintf("lower bound at floor = %.3f", lo))

# --- 8. input guards ---------------------------------------------------------
check("missing source_cell is rejected",
      inherits(try(cluster_bootstrap(df[, c("y","p")], n_boot = 2), silent = TRUE), "try-error"))
# flip ONE row of a cell that has several, so that cell carries both labels
multi <- names(which(table(df$source_cell) > 1))[1]
bad <- df; bad$y[which(bad$source_cell == multi)[1]] <- 1L - bad$y[which(bad$source_cell == multi)[1]]
check("source cell with two labels is rejected",
      inherits(try(cluster_bootstrap(bad, n_boot = 2, verbose = FALSE), silent = TRUE), "try-error"))

cat("\n--- per-fraction AUROC (150 reps, illustrative) ---\n")
d <- summ[summ$metric == "auroc", c("fraction","estimate","lower","upper","width")]
print(d, row.names = FALSE, digits = 3)

cat("\n", if (pass) "ALL CHECKS PASSED\n" else "SOME CHECKS FAILED\n", sep = "")
quit(status = if (pass) 0 else 1)
