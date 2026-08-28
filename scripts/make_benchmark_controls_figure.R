#!/usr/bin/env Rscript
# =============================================================================
# benchmark_controls.png  -  the benchmark gates from METHODS.md section 5
#
#   B-vs-B                 held-out scored-B against disjoint reference-B
#   label shuffle          discrimination must collapse to chance
#   barcode match rate     truth joined to caller output on source_cell
#   reference misclass     reference-B cells called aneuploid
#
# Panels with no data render as "awaiting data" rather than failing, so the
# figure can be built before every control has run.
#
# Note: chr6p/16q/chrX are NOT a control here. METHODS.md section 5 treats them
# as a masked-region diagnostic, not false-positive ground truth, and does not
# require them to be silent.
#
# Usage: Rscript make_benchmark_controls_figure.R [benchmark_metrics.tsv] [out.png]
# =============================================================================
source("benchmark_metrics_schema.R")

args <- commandArgs(trailingOnly = TRUE)
IN  <- if (length(args) >= 1) args[1] else "benchmark_metrics.tsv"
OUT <- if (length(args) >= 2) args[2] else "benchmark_controls.png"

df <- read_metrics(IN, validate = TRUE)
# Row routing lives in benchmark_metrics_schema.R (assign_figure), so the two
# figure scripts cannot disagree about which rows belong where.
ctl <- df[assign_figure(df) == "benchmark_controls", ]
PRELIM <- any(ctl$preliminary)

SITES <- sort(unique(df$site))
SITE_COL <- setNames(c("#C0392B", "#2471A3", "#2E7D57")[seq_along(SITES)], SITES)

awaiting <- function(title, why) {
  plot.new(); box(col = "grey80")
  text(0.5, 0.60, title, font = 2, cex = 1.0)
  text(0.5, 0.44, "awaiting data", col = "grey45", cex = 0.92)
  text(0.5, 0.30, why, col = "grey55", cex = 0.72)
}

#' Bar per site with an interval, a pass threshold, and an explicit verdict.
control_panel <- function(metric, title, note_match = NULL, ylim = c(0, 1),
                          threshold = NULL, direction = c("below", "near"),
                          target = NA, why = "") {
  direction <- match.arg(direction)
  d <- ctl[ctl$metric == metric, ]
  if (!is.null(note_match)) d <- d[grepl(note_match, d$note, ignore.case = TRUE), ]
  if (!nrow(d)) return(awaiting(title, why))

  # one bar per fraction when the control varies across fractions, else per site
  by_fraction <- length(unique(d$fraction)) > 1
  d <- if (by_fraction) d[order(d$fraction), ] else d[order(d$site), ]
  labs <- if (by_fraction) paste0(d$fraction * 100, "%") else d$site
  bp <- barplot(d$estimate, names.arg = labs, ylim = ylim, border = NA,
                col = adjustcolor(SITE_COL[d$site], 0.85),
                ylab = metric, main = title,
                xlab = if (by_fraction) "tumour fraction" else "")
  ok <- !is.na(d$lower)
  if (any(ok)) arrows(bp[ok], d$lower[ok], bp[ok], d$upper[ok],
                      length = 0.04, angle = 90, code = 3, col = "grey25")
  # Value labels: inside the bar when it is tall enough, above the interval when
  # it is short. Placing everything above collides with the tolerance band in the
  # label-shuffle panel and runs off the top where the estimate is 1.000.
  span <- diff(ylim)
  tall <- d$estimate > ylim[1] + span * 0.30
  if (any(tall))
    text(bp[tall], ylim[1] + span * 0.06, sprintf("%.3f", d$estimate[tall]),
         cex = 0.82, col = "white", font = 2)
  if (any(!tall)) {
    top <- pmax(d$estimate[!tall], ifelse(is.na(d$upper[!tall]), -Inf, d$upper[!tall]))
    ly  <- top + span * 0.05
    # nudge clear of a gate line the label would otherwise sit on top of
    guide <- if (direction == "below") threshold else NA_real_
    if (!is.null(guide) && !is.na(guide)) {
      hit <- abs(ly - guide) < span * 0.045
      ly[hit] <- guide + span * 0.06
    }
    text(bp[!tall], ly, sprintf("%.3f", d$estimate[!tall]), cex = 0.82, col = "grey20")
  }

  # "below": threshold is a ceiling. "near": threshold is a tolerance around the
  # expected value, so it must be drawn as a band, not as a line at its own value.
  if (direction == "below" && !is.null(threshold)) {
    abline(h = threshold, lty = 2, col = "#B03A2E")
    # Value labels sit above the line, so the gate annotation goes below it, in
    # whichever horizontal gap has no bar under the line. Anchoring it to a fixed
    # edge puts it on top of a bar whenever the last category is non-zero.
    gy <- threshold - diff(ylim) * 0.055
    clear <- which(d$estimate < gy)            # bars that do not reach the label
    gx <- if (length(clear)) bp[clear[which.min(abs(clear - length(bp) / 2))]]
          else par("usr")[2]
    text(gx, gy, sprintf("gate  %.2f", threshold),
         adj = if (length(clear)) 0.5 else 1, cex = 0.72, col = "#B03A2E")
  }
  if (direction == "near" && !is.na(target)) {
    rect(par("usr")[1], target - threshold, par("usr")[2], target + threshold,
         col = adjustcolor("#B03A2E", 0.08), border = NA)
    abline(h = target, lty = 3, col = "grey30")
    abline(h = c(target - threshold, target + threshold), lty = 2, col = "#B03A2E")
    text(par("usr")[2], target + threshold + diff(ylim) * 0.035,
         sprintf("expected %.2f  \u00b1 %.2f", target, threshold),
         adj = 1, cex = 0.72, col = "#B03A2E")
  }

  # verdict from the interval where there is one, else the point estimate
  bound <- if (direction == "below") ifelse(is.na(d$upper), d$estimate, d$upper) else d$estimate
  pass <- if (direction == "below") all(bound <= threshold) else
          all(abs(d$estimate - target) <= threshold)
  mtext(if (pass) "PASS" else "CHECK", side = 3, line = 0.15, adj = 1,
        cex = 0.78, col = if (pass) "#2E7D57" else "#B03A2E", font = 2)
}

png(OUT, width = 2000, height = 1400, res = 200)
layout(matrix(1:4, 2, 2, byrow = TRUE))
par(mar = c(4.0, 4.6, 3.2, 1.4), oma = c(2.4, 0, 1.6, 0),
    cex.axis = 0.82, cex.lab = 0.92, cex.main = 0.98)

control_panel("false_aneuploid_rate",
              "A   B-vs-B: false aneuploid rate", note_match = "B-vs-B",
              ylim = c(0, 0.25), threshold = 0.10, direction = "below",
              why = "needs reference split IDs from the mixture builder")

control_panel("auroc", "B   Label shuffle: discrimination", note_match = "shuffle",
              ylim = c(0, 1), target = 0.5, threshold = 0.05, direction = "near",
              why = "run once caller scores exist")

control_panel("barcode_match_rate", "C   Truth join: barcode match rate",
              ylim = c(0, 1.25), target = 1.0, threshold = 0.10, direction = "near",
              why = "computed when caller output is joined to truth")

control_panel("reference_misclass_rate", "D   Reference-B misclassification",
              ylim = c(0, 0.25), threshold = 0.05, direction = "below",
              why = "needs reference split IDs from the mixture builder")

mtext("Benchmark controls (METHODS.md section 5)", side = 3, line = 0.1,
      outer = TRUE, cex = 0.95, font = 2)
mtext(sprintf("6p / 16q / chrX are a masked-region diagnostic, not a control%s",
              if (PRELIM) "   |   PRELIMINARY: single replicate" else ""),
      side = 1, line = 0.8, outer = TRUE, cex = 0.78,
      col = if (PRELIM) "#B03A2E" else "grey35")
invisible(dev.off())

write.table(ctl, sub("\\.png$", "_source_data.tsv", OUT), sep = "\t",
            row.names = FALSE, quote = FALSE, na = "NA")
cat("wrote", OUT, "and", sub("\\.png$", "_source_data.tsv", OUT), "\n")
