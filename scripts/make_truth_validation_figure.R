#!/usr/bin/env Rscript
# =============================================================================
# truth_validation.png  -  does our arm-level truth table hold up?
#
# Panel A  deposited consensus arm states across the genome
# Panel B  the consensus against ascatNgs, arm by arm
#
# The ascatNgs values for the masked arms are carried in the subtitle rather than
# a third panel: a bar chart against a structural absence reads as "low" when it
# means "not measured".
#
# Base R graphics only, so this runs in any caller container.
# Inputs : truth_arms.csv, ascat_arms.csv
# Outputs: truth_validation.png, truth_validation_source_data.tsv
# =============================================================================

t <- read.csv("truth_arms.csv", stringsAsFactors = FALSE)
a <- read.csv("ascat_arms.csv", stringsAsFactors = FALSE)

m <- merge(t[, c("arm_id","chrom","callable_mb","frac_gain","frac_loss","frac_loh",
                 "mean_cn","state_call","masked")],
           a[, c("arm_id","frac_gain","frac_loss","frac_loh","mean_cn","state_call")],
           by = "arm_id", suffixes = c("_consensus", "_ascat"))

chr_n <- suppressWarnings(as.integer(sub("[pq]$", "", m$arm_id)))
chr_n[is.na(chr_n)] <- 99
m <- m[order(chr_n, m$arm_id), ]

scored <- m[!m$masked & !is.na(m$state_call_consensus), ]
agree  <- scored$state_call_consensus == scored$state_call_ascat
r_p <- cor(scored$mean_cn_consensus, scored$mean_cn_ascat)
r_s <- cor(scored$mean_cn_consensus, scored$mean_cn_ascat, method = "spearman")
PLOIDY <- 3.087
YMAX   <- 6.0

COL <- c(gain = "#C0392B", loss = "#2471A3", neutral = "#7F8C8D", mixed = "#E0A458")

#' One genome-ordered arm track. Shared y-axis so A and B are directly comparable.
arm_track <- function(cn, state, title, mark_nocall = TRUE) {
  cols <- ifelse(is.na(state), "#FFFFFF", COL[state])
  bp <- barplot(ifelse(is.na(cn), 0, cn), col = cols,
                border = ifelse(is.na(state), "grey35", NA),
                ylim = c(0, YMAX), las = 2, names.arg = m$arm_id,
                ylab = "mean copy number", xlab = "", main = title)
  abline(h = 2, lty = 2, col = "grey45")
  abline(h = PLOIDY, lty = 3, col = "#C0392B")
  if (mark_nocall) {
    for (i in which(is.na(state))) {
      rect(bp[i] - 0.42, 0, bp[i] + 0.42, YMAX, col = "#EDEDED",
           border = "grey55", lty = 3)
      text(bp[i], YMAX / 2, "no call", srt = 90, cex = 0.62, col = "grey30")
    }
  }
  # the masked arms get a marker in every track so they are findable in both
  for (i in which(m$masked)) {
    axis(1, at = bp[i], labels = FALSE, tcl = 0.35, lwd = 0, lwd.ticks = 2,
         col.ticks = "#C0392B")
  }
  invisible(bp)
}

png("truth_validation.png", width = 2400, height = 1750, res = 200)
# A spans the full width and is the taller panel; its key sits in the space
# beside the square scatter rather than on top of the bars.
layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE),
       widths = c(0.72, 2), heights = c(1.28, 1))
par(mar = c(4.4, 4.6, 3.0, 1.2), cex.axis = 0.78, cex.lab = 0.95)

# --- A: consensus ------------------------------------------------------------
arm_track(m$mean_cn_consensus, m$state_call_consensus,
          "A   Deposited SEQC2 consensus: arm-level copy-number truth (HCC1395, GRCh38)")

# --- key panel for A ---------------------------------------------------------
par(mar = c(0.4, 1.2, 0.4, 0.4), pty = "m", xpd = NA)
plot.new()
text(0, 1.00, "Panel A key", adj = c(0, 1), font = 2, cex = 1.05)
legend(0, 0.945, bty = "n", cex = 0.95, border = NA, y.intersp = 1.15,
       fill = c(COL[c("gain", "mixed", "neutral", "loss")], "#EDEDED"),
       legend = c("gain", "mixed", "neutral",
                  "loss (no arm reaches 50%)", "no call (masked)"))
legend(0, 0.575, bty = "n", cex = 0.95, lty = c(2, 3), seg.len = 2.2,
       col = c("grey45", "#C0392B"), y.intersp = 1.05,
       legend = c("diploid baseline (CN 2)",
                  sprintf("genome ploidy %.2f", PLOIDY)))
note <- sprintf(paste("6p and 16q carry no consensus calls across %.0f and %.0f Mb",
                      "of callable sequence. ascatNgs reports mean CN %.2f and",
                      "%.2f there, the two highest arms in the genome."),
                m$callable_mb[m$arm_id == "6p"], m$callable_mb[m$arm_id == "16q"],
                m$mean_cn_ascat[m$arm_id == "6p"], m$mean_cn_ascat[m$arm_id == "16q"])
text(0, 0.400, paste(strwrap(note, width = 34), collapse = "\n"),
     adj = c(0, 1), cex = 0.88, col = "grey25")

# --- B: the two against each other -------------------------------------------
par(mar = c(4.4, 4.8, 2.8, 1.0), pty = "s", xpd = FALSE)
lim <- c(1.5, 7.3)
plot(scored$mean_cn_consensus, scored$mean_cn_ascat, xlim = lim, ylim = lim,
     pch = 21, bg = ifelse(agree, "#2E7D57", "#E0A458"), col = "grey25", cex = 1.3,
     xlab = "consensus mean CN", ylab = "ascatNgs mean CN",
     main = "B   Consensus vs ascatNgs (37 scored arms; 6p and 16q excluded, no consensus call)")
abline(0, 1, col = "grey55", lty = 2)
dis <- scored[!agree, ]
if (nrow(dis)) text(dis$mean_cn_consensus, dis$mean_cn_ascat, dis$arm_id,
                    pos = 4, cex = 0.7, col = "#8A5A00")
legend("topleft", bty = "n", cex = 0.78,
       legend = c(sprintf("state agreement  %d / %d  (%.1f%%)",
                          sum(agree), length(agree), 100 * mean(agree)),
                  sprintf("Pearson r = %.3f", r_p),
                  sprintf("Spearman r = %.3f", r_s),
                  sprintf("median |difference| = %.2f copies",
                          median(abs(scored$mean_cn_consensus - scored$mean_cn_ascat)))))
legend("bottomright", bty = "n", cex = 0.74, pt.cex = 1.2, pch = 21, col = "grey25",
       pt.bg = c("#2E7D57", "#E0A458"), legend = c("state agrees", "state differs"))

# --- separators ---------------------------------------------------------------
# par(fig=) does not override an active layout(), so convert normalized device
# coordinates into the current panel's user space and draw with xpd = NA.
Y_SPLIT <- 1 / (1.28 + 1)          # row boundary: A above, key + B below
X_SPLIT <- 0.72 / (0.72 + 2)       # column boundary: key left, B right
par(xpd = NA)
gx <- function(f) grconvertX(f, from = "ndc", to = "user")
gy <- function(f) grconvertY(f, from = "ndc", to = "user")
# Only a vertical rule between the panel A key and panel B. No horizontal rule:
# the key belongs to A, so separating them would imply otherwise.
segments(gx(X_SPLIT), gy(0.03), gx(X_SPLIT), gy(Y_SPLIT - 0.005), lty = 2, col = "grey65")
par(xpd = FALSE)

invisible(dev.off())

src <- m
src$agrees <- ifelse(src$masked | is.na(src$state_call_consensus), NA,
                     src$state_call_consensus == src$state_call_ascat)
write.table(src, "truth_validation_source_data.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE, na = "NA")

cat(sprintf("scored arms %d | agreement %d/%d (%.1f%%) | Pearson %.3f | Spearman %.3f\n",
            nrow(scored), sum(agree), length(agree), 100 * mean(agree), r_p, r_s))
cat("no-call arms:", paste(m$arm_id[is.na(m$state_call_consensus)], collapse = ", "), "\n")
