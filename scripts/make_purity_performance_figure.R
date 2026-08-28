#!/usr/bin/env Rscript
# =============================================================================
# purity_performance.png  -  caller performance against tumour fraction
#
# Reads benchmark_metrics.tsv (see benchmark_metrics_schema.R). Every point carries its own
# interval from that file; nothing is recomputed here.
#
# Usage: Rscript make_purity_performance_figure.R [benchmark_metrics.tsv] [out.png]
# =============================================================================
source("benchmark_metrics_schema.R")

args <- commandArgs(trailingOnly = TRUE)
IN   <- if (length(args) >= 1) args[1] else "benchmark_metrics.tsv"
OUT  <- if (length(args) >= 2) args[2] else "purity_performance.png"

df <- read_metrics(IN, validate = TRUE)
# Controls now carry a fraction too (p00 is the B-vs-B run), so exclude them by
# metric and by note rather than assuming a control has no fraction.
# Row routing lives in benchmark_metrics_schema.R, so the two figure scripts
# cannot disagree about which rows belong where.
perf <- df[assign_figure(df) == "purity_performance", ]
if (!nrow(perf)) stop("No rows with a tumour fraction in ", IN)

SITES <- sort(unique(perf$site))
SITE_COL <- setNames(c("#C0392B", "#2471A3", "#2E7D57")[seq_along(SITES)], SITES)
PRELIM <- any(perf$preliminary)

#' One metric against tumour fraction, one line per site, interval as a ribbon.
panel <- function(metric, title, ylim = c(0, 1), ref = NULL, ref_lab = NULL,
                  identity = FALSE) {
  d <- perf[perf$metric == metric, ]
  if (!nrow(d)) {
    plot.new(); box(col = "grey80")
    text(0.5, 0.55, title, font = 2, cex = 0.95)
    text(0.5, 0.42, "awaiting data", col = "grey45", cex = 0.9)
    return(invisible(NULL))
  }
  fr <- sort(unique(d$fraction))
  plot(NA, xlim = range(fr), ylim = ylim, log = "x", xaxt = "n",
       xlab = "tumour fraction", ylab = metric, main = title)
  axis(1, at = fr, labels = paste0(fr * 100, "%"))
  if (!is.null(ref)) {
    abline(h = ref, lty = 3, col = "grey45")
    if (!is.null(ref_lab)) text(min(fr), ref + 0.035, ref_lab, adj = 0,
                                cex = 0.7, col = "grey40")
  }
  if (identity) abline(0, 1, lty = 2, col = "grey60")   # only sane on a linear axis

  for (s in SITES) {
    ds <- d[d$site == s, ]
    ds <- ds[order(ds$fraction), ]
    if (!nrow(ds)) next
    ok <- !is.na(ds$lower) & !is.na(ds$upper)
    if (any(ok)) {
      lo <- pmax(ds$lower[ok], ylim[1]); hi <- pmin(ds$upper[ok], ylim[2])
      polygon(c(ds$fraction[ok], rev(ds$fraction[ok])), c(lo, rev(hi)),
              col = adjustcolor(SITE_COL[[s]], 0.18), border = NA)
    }
    lines(ds$fraction, ds$estimate, col = SITE_COL[[s]], lwd = 2)
    points(ds$fraction, ds$estimate, pch = 21, bg = SITE_COL[[s]], col = "white", cex = 1.1)
  }
}

png(OUT, width = 2400, height = 1500, res = 200)
layout(matrix(1:6, 2, 3, byrow = TRUE))
par(mar = c(4.2, 4.4, 3.0, 1.0), oma = c(2.4, 0, 0, 0),
    cex.axis = 0.8, cex.lab = 0.92, cex.main = 0.98)

panel("sensitivity", "A   Sensitivity")
legend("topleft", bty = "n", cex = 0.82, lwd = 2, col = SITE_COL[SITES], legend = SITES)
panel("specificity", "B   Specificity", ylim = c(0.75, 1))
panel("precision",   "C   Precision")
panel("f1",          "D   F1")
panel("auroc",       "E   AUROC (needs continuous burden score)", ylim = c(0.4, 1),
      ref = 0.5, ref_lab = "chance")

# recovered vs nominal purity: the sanity check that a caller is in the right
# ballpark at all. Linear axes so the identity line is meaningful.
d <- perf[perf$metric == "frac_aneuploid", ]
if (nrow(d)) {
  plot(NA, xlim = c(0, 0.85), ylim = c(0, 0.85),
       xlab = "nominal tumour fraction", ylab = "fraction called aneuploid",
       main = "F   Purity recovery")
  abline(0, 1, lty = 2, col = "grey55")
  for (s in SITES) {
    ds <- d[d$site == s, ]; ds <- ds[order(ds$fraction), ]
    ok <- !is.na(ds$lower)
    if (any(ok)) arrows(ds$fraction[ok], ds$lower[ok], ds$fraction[ok], ds$upper[ok],
                        length = 0.03, angle = 90, code = 3, col = SITE_COL[[s]])
    lines(ds$fraction, ds$estimate, col = SITE_COL[[s]], lwd = 2)
    points(ds$fraction, ds$estimate, pch = 21, bg = SITE_COL[[s]], col = "white", cex = 1.1)
  }
  text(0.44, 0.08, "dashed: perfect recovery", cex = 0.72, col = "grey40")
} else {
  plot.new(); box(col = "grey80")
  text(0.5, 0.5, "F   Purity recovery\nawaiting data", col = "grey45", cex = 0.9)
}

mtext(sprintf("CopyKAT, 10x expression arm%s",
              if (PRELIM) "   |   PRELIMINARY: single replicate, uncertainty not final" else ""),
      side = 1, line = 0.9, outer = TRUE, cex = 0.8,
      col = if (PRELIM) "#B03A2E" else "grey30")
invisible(dev.off())

write.table(perf, sub("\\.png$", "_source_data.tsv", OUT), sep = "\t",
            row.names = FALSE, quote = FALSE, na = "NA")
cat("wrote", OUT, "and", sub("\\.png$", "_source_data.tsv", OUT), "\n")
