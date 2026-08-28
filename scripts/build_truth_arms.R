#!/usr/bin/env Rscript
# =============================================================================
# Build the GRCh38 arm-level copy-number truth table for the scRNA CNV benchmark.
#
# INPUTS  (all in one directory, passed as the first argument)
#   ngs_benchmark_cnv_gain_loss_loh.bed        Masood consensus, 688 segments
#   exclusion.bed                              regions the benchmark makes no claim on
#   Additional file 3_cnv_gain_cn_median.txt   median CN inside gain regions
#   Additional file 4_cnv_loss_cn_median.txt   median CN inside loss regions
#   WGS_IL_T_1_copynumber_caveman.csv          ascatNgs single replicate (cross-check)
#
# OUTPUTS
#   truth_arms.csv       arm-level truth from the consensus   <- the deliverable
#   ascat_arms.csv       same table built from ascatNgs
#   arm_comparison.csv   the two side by side, with agreement
#
# Usage:  Rscript build_truth_granges.R /path/to/files
# =============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)
DIR  <- if (length(args)) args[1] else "."

# --- knobs -------------------------------------------------------------------
BASELINE_CN   <- 2      # gain/loss are called against a diploid baseline
DOMINANT_FRAC <- 0.50   # fraction of a callable arm needed to call gain or loss
NEUTRAL_MAX   <- 0.20   # below this altered fraction, the arm is neutral
# Aberrant in the HCC1395BL normal line, so tumour-vs-normal interpretation here
# is intrinsically ambiguous. Per METHODS.md section 3 these are a masked-region
# DIAGNOSTIC - they are not false-positive ground truth and are not required to
# be silent.
MASKED_ARMS   <- c("6p", "16q", "Xp", "Xq")

# Published spans from Masood Table S2 (Mb) - used as a hard QC assertion.
# QC targets are the DEPOSIT's own values, as a regression check that the parse
# has not drifted - NOT as agreement with any table in the paper.
#
# Masood et al. gives two figures for the high-confidence set. Page 9 (end of the
# orthogonal-validation section): 346 / 33 / 320 at 1525.6 / 87.9 / 1490.4 Mb.
# Table S2: 340 / 33 / 315 at 1518.5 / 79.7 / 1456.1 Mb. The deposited BED matches
# page 9 on spans and Table S2 on counts, i.e. the full benchmark with contiguous
# segments merged; Table S2 appears to be the strong-evidence subset. See
# truth_arms_PROVENANCE.md.
EXPECT_N  <- c(gain = 340,    loss = 33,   loh = 315)
EXPECT_MB <- c(gain = 1525.5, loss = 87.9, loh = 1490.4)
QC_TOL    <- 0.02       # allow 2% drift on spans

f <- function(x) file.path(DIR, x)

# =============================================================================
# STEP 0 - Reference coordinates
#
# Hardcoded rather than fetched: Seqinfo(genome="hg38") calls out to UCSC and
# the Cloud Workstation may have no egress. These lengths are the standard
# GRCh38 primary assembly values and match exclusion.bed's terminal coordinates.
# =============================================================================
CHROM_LEN <- c(
  chr1  = 248956422, chr2  = 242193529, chr3  = 198295559, chr4  = 190214555,
  chr5  = 181538259, chr6  = 170805979, chr7  = 159345973, chr8  = 145138636,
  chr9  = 138394717, chr10 = 133797422, chr11 = 135086622, chr12 = 133275309,
  chr13 = 114364328, chr14 = 107043718, chr15 = 101991189, chr16 =  90338345,
  chr17 =  83257441, chr18 =  80373285, chr19 =  58617616, chr20 =  64444167,
  chr21 =  46709983, chr22 =  50818468, chrX  = 156040895, chrY  =  57227415
)

# Centromere spans (UCSC hg38 acen, rounded). Only used to split p from q, so
# a few hundred kb of imprecision cannot move an arm-level fraction materially.
CEN <- data.frame(
  chrom = paste0("chr", c(1:22, "X")),
  start = c(121700000, 91800000, 87800000, 48200000, 46100000, 58500000,
             58100000, 43200000, 42200000, 38000000, 51000000, 33200000,
             16500000, 16100000, 17500000, 35300000, 22700000, 15400000,
             24200000, 25700000, 10900000, 13700000, 58100000),
  end   = c(125100000, 96000000, 94000000, 51800000, 50000000, 62600000,
             62100000, 47200000, 45500000, 41600000, 55800000, 37800000,
             18900000, 18200000, 20500000, 38400000, 27400000, 21500000,
             28100000, 30400000, 13000000, 17400000, 63800000),
  stringsAsFactors = FALSE
)
# p-arms of acrocentric chromosomes are satellite repeat, not scoreable.
ACRO <- paste0("chr", c(13, 14, 15, 21, 22))

# =============================================================================
# STEP 1 - Load the consensus and the exclusion set
#
# import() converts BED's 0-based half-open coordinates to R's 1-based inclusive
# convention automatically. Do NOT apply any manual +1/-1 on top of this.
# =============================================================================
cat("== STEP 1: load ==\n")

benchmark <- import(f("ngs_benchmark_cnv_gain_loss_loh.bed"), format = "bed")
exclusion <- import(f("exclusion.bed"), format = "bed")

# The state (gain/loss/LOH) arrives in $name because that is BED column 4.
benchmark$state <- tolower(benchmark$name)
stopifnot(all(benchmark$state %in% c("gain", "loss", "loh")))

# Scope = whatever chromosomes the consensus actually covers. It is autosomes
# only; X and Y are absent, and absent means "no call", not "neutral".
chroms <- seqlevels(benchmark)
stopifnot(all(chroms %in% names(CHROM_LEN)))

# Attach chromosome lengths to every object. Two reasons: GenomicRanges then
# rejects any interval running off the end of a chromosome (a loud failure
# instead of a silent coordinate bug), and operations between objects stop
# emitting Seqinfo-mismatch warnings.
set_lengths <- function(gr) {
  seqlengths(gr) <- CHROM_LEN[seqlevels(gr)]
  gr
}
benchmark <- set_lengths(benchmark)
exclusion <- set_lengths(exclusion)
cat(sprintf("  benchmark: %d segments across %d chromosomes\n",
            length(benchmark), length(chroms)))
missing_chr <- setdiff(names(CHROM_LEN), chroms)
if (length(missing_chr))
  cat("  not covered by the consensus (will be reported as NA):",
      paste(missing_chr, collapse = ", "), "\n")

# =============================================================================
# STEP 2 - QC against the published table, BEFORE anything is transformed
#
# If the parse is wrong, everything downstream is wrong in ways that still look
# plausible. Fail here instead.
# =============================================================================
cat("\n== STEP 2: QC vs Masood Table S2 ==\n")
for (st in names(EXPECT_N)) {
  g   <- benchmark[benchmark$state == st]
  n   <- length(g)
  mb  <- sum(as.numeric(width(g))) / 1e6
  dev <- abs(mb - EXPECT_MB[[st]]) / EXPECT_MB[[st]]
  cat(sprintf("  %-5s %4d segments (expect %3d)  %8.1f Mb (expect %7.1f, %+.2f%%)\n",
              st, n, EXPECT_N[[st]], mb, EXPECT_MB[[st]], 100 * (mb - EXPECT_MB[[st]]) / EXPECT_MB[[st]]))
  if (n != EXPECT_N[[st]])
    stop("Segment count for '", st, "' is ", n, ", expected ", EXPECT_N[[st]],
         ". The file or the state column is not what we think it is.")
  if (dev > QC_TOL)
    warning("Span for '", st, "' is off by ", round(100 * dev, 1), "%.")
}

# =============================================================================
# STEP 3 - Define callable space
#
# callable = covered chromosomes minus exclusion. The exclusion set contains
# both assembly gaps (chr1 1-10000 is the telomeric N-run) and uncertainty
# buffers between calls, so subtracting it is what makes "not called" safely
# mean "neutral" rather than "unknown".
# =============================================================================
cat("\n== STEP 3: callable space ==\n")

genome_gr <- set_lengths(GRanges(chroms, IRanges(1, CHROM_LEN[chroms])))
callable  <- GenomicRanges::setdiff(genome_gr, exclusion, ignore.strand = TRUE)

cat(sprintf("  genome (covered chroms): %8.1f Mb\n", sum(as.numeric(width(genome_gr))) / 1e6))
cat(sprintf("  after exclusion:         %8.1f Mb callable (%.1f%% retained)\n",
            sum(as.numeric(width(callable))) / 1e6,
            100 * sum(as.numeric(width(callable))) / sum(as.numeric(width(genome_gr)))))

# =============================================================================
# STEP 4 - The three state tracks
#
# reduce() merges touching/overlapping intervals within a state; intersect()
# clips each to callable space. gain and loss are mutually exclusive; LOH is a
# SEPARATE AXIS - a region can be gained and LOH at once - so it never competes
# with gain/loss for the arm call.
# =============================================================================
cat("\n== STEP 4: state tracks ==\n")

track <- function(st) GenomicRanges::intersect(
  reduce(benchmark[benchmark$state == st]), callable, ignore.strand = TRUE)

gain <- track("gain"); loss <- track("loss"); loh <- track("loh")

# gain and loss should be mutually exclusive. A trace of overlap is a rounding
# artefact at segment boundaries; a lot of it means the parse is wrong.
ov_gl <- sum(as.numeric(width(GenomicRanges::intersect(gain, loss, ignore.strand = TRUE))))
if (ov_gl > 0) {
  frac <- ov_gl / sum(as.numeric(width(gain)))
  msg <- sprintf("gain and loss overlap over %.2f Mb (%.2f%% of gain)", ov_gl / 1e6, 100 * frac)
  if (frac > 0.01) stop(msg, " - a region cannot be both. Check the parse.")
  warning(msg, " - small enough to be boundary rounding; loss takes precedence.")
  gain <- GenomicRanges::setdiff(gain, loss, ignore.strand = TRUE)
}

for (nm in c("gain", "loss", "loh"))
  cat(sprintf("  %-5s %8.1f Mb callable\n", nm,
              sum(as.numeric(width(get(nm)))) / 1e6))

# =============================================================================
# STEP 5 - Copy-number magnitude
#
# The median-CN files give a continuous value per region, which lets you
# correlate against inferCNV/CopyKAT continuous scores rather than only checking
# categorical agreement. The track MUST be disjoint: overlapping intervals
# double-count base pairs and silently inflate mean_cn.
# =============================================================================
cat("\n== STEP 5: CN track ==\n")

read_cn <- function(path) {
  gr <- import(path, format = "bed")       # column 4 -> $name
  gr$cn <- suppressWarnings(as.numeric(gr$name))
  if (any(is.na(gr$cn)))
    stop(basename(path), ": column 4 is not numeric - wrong column or wrong file.")
  gr$name <- NULL
  gr
}

cn_track <- c(read_cn(f("Additional file 3_cnv_gain_cn_median.txt")),
              read_cn(f("Additional file 4_cnv_loss_cn_median.txt")))
cn_track <- cn_track[as.character(seqnames(cn_track)) %in% chroms]
seqlevels(cn_track, pruning.mode = "coarse") <- chroms
cn_track <- set_lengths(cn_track)

if (!isDisjoint(cn_track)) {
  d    <- disjoin(cn_track)
  hits <- findOverlaps(d, cn_track)
  d$cn <- as.numeric(tapply(cn_track$cn[subjectHits(hits)], queryHits(hits), mean))
  cat(sprintf("  flattened %d overlapping intervals -> %d disjoint\n",
              length(cn_track), length(d)))
  cn_track <- d
}
cat(sprintf("  %d intervals, CN range %.1f - %.1f\n",
            length(cn_track), min(cn_track$cn), max(cn_track$cn)))

# =============================================================================
# STEP 6 - Build arm ranges
# =============================================================================
build_arms <- function(chroms) {
  d <- CEN[CEN$chrom %in% chroms, , drop = FALSE]
  id <- sub("^chr", "", d$chrom)

  p_ok <- !(d$chrom %in% ACRO)          # acrocentric p-arms are satellite repeat
  chrom_v <- c(d$chrom[p_ok], d$chrom)
  start_v <- c(rep(1, sum(p_ok)),        d$end)
  end_v   <- c(d$start[p_ok],            CHROM_LEN[d$chrom])
  arm_v   <- c(paste0(id[p_ok], "p"),    paste0(id, "q"))

  gr <- GRanges(chrom_v, IRanges(start_v, end_v), arm_id = arm_v)
  sort(gr)
}
arms <- set_lengths(build_arms(chroms))
cat(sprintf("\n== STEP 6: %d arms in scope ==\n", length(arms)))

# =============================================================================
# STEP 7 - Collapse to arms
#
# Every number is a width ratio over the CALLABLE part of the arm, so excluded
# base pairs neither count as altered nor dilute the denominator.
# =============================================================================
collapse_arms <- function(arms, gain, loss, loh, callable, cn_track, label,
                          all_calls = NULL) {
  rows <- lapply(seq_along(arms), function(i) {
    a  <- arms[i]
    ca <- GenomicRanges::intersect(a, callable, ignore.strand = TRUE)
    W  <- sum(as.numeric(width(ca)))

    if (W == 0) {
      return(data.frame(arm_id = a$arm_id, chrom = as.character(seqnames(a)),
                        callable_mb = 0, n_segments = 0L,
                        frac_gain = NA_real_, frac_loss = NA_real_,
                        frac_loh = NA_real_, frac_neutral = NA_real_,
                        mean_cn = NA_real_, state_call = NA_character_,
                        masked = a$arm_id %in% MASKED_ARMS, source = label,
                        stringsAsFactors = FALSE))
    }

    ov_bp <- function(x) sum(as.numeric(width(
      GenomicRanges::intersect(ca, x, ignore.strand = TRUE))))
    fg <- ov_bp(gain) / W
    fl <- ov_bp(loss) / W
    fh <- ov_bp(loh)  / W

    # length-weighted mean CN; callable bp the CN track does not cover sit at
    # the baseline, because "not called gain or loss" means diploid here
    hits <- findOverlaps(cn_track, ca)
    if (length(hits)) {
      pin <- pintersect(cn_track[queryHits(hits)], ca[subjectHits(hits)])
      den <- sum(as.numeric(width(pin)))
      num <- sum(as.numeric(width(pin)) * cn_track$cn[queryHits(hits)])
    } else { den <- 0; num <- 0 }
    mean_cn <- (num + BASELINE_CN * max(0, W - den)) / W

    # How many source segments touch this arm? Zero calls across a callable arm
    # is NO CALL, not a diploid assertion - the distinction that separates 6p and
    # 16q (silent in the consensus) from genuinely copy-neutral arms like 8p,
    # which carries 97.5% LOH. Reporting these as "neutral, mean_cn 2" would
    # assert the truth set says they are diploid. It does not.
    n_seg <- if (is.null(all_calls)) NA_integer_ else
      length(subsetByOverlaps(all_calls, ca, ignore.strand = TRUE))

    altered <- fg + fl
    state <- if (!is.na(n_seg) && n_seg == 0L) {
               NA_character_
             } else if (fg >= DOMINANT_FRAC && fg >= fl) {
               "gain"
             } else if (fl >= DOMINANT_FRAC) {
               "loss"
             } else if (altered < NEUTRAL_MAX) {
               "neutral"
             } else {
               "mixed"
             }
    if (!is.na(n_seg) && n_seg == 0L) mean_cn <- NA_real_

    data.frame(arm_id = a$arm_id, chrom = as.character(seqnames(a)),
               callable_mb = round(W / 1e6, 2), n_segments = n_seg,
               frac_gain = round(fg, 4), frac_loss = round(fl, 4),
               frac_loh = round(fh, 4), frac_neutral = round(max(0, 1 - altered), 4),
               mean_cn = round(mean_cn, 3), state_call = state,
               masked = a$arm_id %in% MASKED_ARMS, source = label,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

truth_arms <- collapse_arms(arms, gain, loss, loh, callable, cn_track,
                            "masood_consensus", all_calls = benchmark)

cat("\n== STEP 7: consensus arm states ==\n")
print(table(truth_arms$state_call,
            ifelse(truth_arms$masked, "MASKED", "scored"), useNA = "ifany"))

# =============================================================================
# STEP 8 - ascatNgs cross-check
#
# The caveman CSV is headerless and 1-based inclusive (unlike BED), so its
# coordinates go straight into IRanges with no shift. Columns:
#   1 index | 2 chr | 3 start | 4 end | 5 normal CN | 6 normal minor
#   7 tumour total CN | 8 tumour minor CN
# =============================================================================
ascat_arms <- NULL
if (file.exists(f("WGS_IL_T_1_copynumber_caveman.csv"))) {
  cat("\n== STEP 8: ascatNgs cross-check ==\n")
  a <- read.csv(f("WGS_IL_T_1_copynumber_caveman.csv"), header = FALSE,
                stringsAsFactors = FALSE)
  names(a) <- c("idx", "chrom", "start", "end", "n_cn", "n_minor", "t_cn", "t_minor")

  ag <- GRanges(a$chrom, IRanges(a$start, a$end), cn = a$t_cn, minor = a$t_minor)
  ag <- ag[as.character(seqnames(ag)) %in% chroms]      # autosomes, to match scope
  seqlevels(ag, pruning.mode = "coarse") <- chroms

  ploidy <- sum(as.numeric(width(ag)) * ag$cn) / sum(as.numeric(width(ag)))
  cat(sprintf("  %d segments, length-weighted ploidy %.3f\n", length(ag), ploidy))

  a_gain <- GenomicRanges::intersect(reduce(ag[ag$cn > BASELINE_CN]), callable, ignore.strand = TRUE)
  a_loss <- GenomicRanges::intersect(reduce(ag[ag$cn < BASELINE_CN]), callable, ignore.strand = TRUE)
  a_loh  <- GenomicRanges::intersect(reduce(ag[ag$minor == 0]),       callable, ignore.strand = TRUE)

  # CN track built from ag itself - it already tiles the genome and is disjoint,
  # so LOH must NOT be appended to it (same intervals, would double-count).
  a_cn <- granges(ag); a_cn$cn <- ag$cn
  if (!isDisjoint(a_cn)) {
    d <- disjoin(a_cn); h <- findOverlaps(d, a_cn)
    d$cn <- as.numeric(tapply(a_cn$cn[subjectHits(h)], queryHits(h), mean))
    a_cn <- d
  }

  # ascat tiles the genome, so every arm has calls and none go NA - correct,
  # since ascat does make claims where the consensus is silent.
  ascat_arms <- collapse_arms(arms, a_gain, a_loss, a_loh, callable, a_cn,
                              "ascat_IL_T_1", all_calls = ag)
}

# =============================================================================
# STEP 9 - Write and compare
# =============================================================================
write.csv(truth_arms, f("truth_arms.csv"), row.names = FALSE)

if (!is.null(ascat_arms)) {
  write.csv(ascat_arms, f("ascat_arms.csv"), row.names = FALSE)

  cmp <- merge(truth_arms[, c("arm_id", "masked", "state_call", "mean_cn")],
               ascat_arms[, c("arm_id", "state_call", "mean_cn")],
               by = "arm_id", suffixes = c("_consensus", "_ascat"))
  cmp$agree <- cmp$state_call_consensus == cmp$state_call_ascat
  write.csv(cmp, f("arm_comparison.csv"), row.names = FALSE)

  sc <- cmp[!cmp$masked & !is.na(cmp$agree), ]
  cat(sprintf("\n== STEP 9: cross-check on %d scored arms ==\n", nrow(sc)))
  cat(sprintf("  state agreement : %d/%d (%.0f%%)\n",
              sum(sc$agree), nrow(sc), 100 * mean(sc$agree)))
  cat(sprintf("  mean_cn Spearman: %.3f\n",
              suppressWarnings(cor(sc$mean_cn_consensus, sc$mean_cn_ascat,
                                   method = "spearman", use = "complete.obs"))))

  mk <- cmp[cmp$masked, ]
  if (nrow(mk)) {
    cat("\n  masked arms (diagnostic only - not FP ground truth, not required to be silent):\n")
    print(mk[, c("arm_id", "state_call_consensus", "state_call_ascat",
                 "mean_cn_consensus", "mean_cn_ascat")], row.names = FALSE)
  }
}

cat("\nWrote truth_arms.csv",
    if (!is.null(ascat_arms)) ", ascat_arms.csv, arm_comparison.csv" else "",
    "\n", sep = "")
