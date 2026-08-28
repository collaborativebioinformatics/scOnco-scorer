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

segs_dir <- get_arg("--segs_dir", "/w/c1_loh")
segments_file <- get_arg("--segments", NA_character_)
truth_file <- get_arg("--truth", "/w/truth_arms.csv")
outdir <- get_arg("--outdir", segs_dir)
truth_threshold <- as.numeric(get_arg("--truth_threshold", "0.5"))
call_threshold <- as.numeric(get_arg("--call_threshold", "0.5"))
posterior_threshold <- as.numeric(get_arg("--posterior_threshold", "0.5"))

if (!requireNamespace("numbat", quietly = TRUE)) {
  stop(
    "The numbat package is required. Run in pkharchenkolab/numbat-rbase:latest.",
    call. = FALSE
  )
}

for (x in c(truth_threshold, call_threshold, posterior_threshold)) {
  if (!is.finite(x) || x < 0 || x > 1) stop("thresholds must be in [0,1]", call. = FALSE)
}

if (!file.exists(truth_file)) stop("truth file does not exist: ", truth_file, call. = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

normalize_chr <- function(x) {
  z <- sub("^chr", "", as.character(x), ignore.case = TRUE)
  z[z %in% c("M", "MT")] <- "M"
  z
}

normalize_arm_id <- function(x) {
  z <- tolower(trimws(as.character(x)))
  sub("^chr", "", z)
}

as_logical_strict <- function(x, label) {
  if (is.logical(x)) {
    if (anyNA(x)) stop(label, " contains NA", call. = FALSE)
    return(x)
  }
  z <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(z))
  out[z %in% c("true", "t", "1", "yes", "y")] <- TRUE
  out[z %in% c("false", "f", "0", "no", "n")] <- FALSE
  if (anyNA(out)) {
    stop(
      label, " contains unrecognized values: ",
      paste(unique(z[is.na(out)]), collapse = ", "),
      call. = FALSE
    )
  }
  out
}

pick_col <- function(dt, candidates, object_name) {
  low <- tolower(names(dt))
  hit <- match(tolower(candidates), low)
  hit <- hit[!is.na(hit)]
  if (!length(hit)) {
    stop(
      object_name, " missing expected column; tried: ",
      paste(candidates, collapse = ", "),
      "; found: ", paste(names(dt), collapse = ", "),
      call. = FALSE
    )
  }
  names(dt)[hit[1L]]
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(x)
}

union_length <- function(start, end) {
  ok <- is.finite(start) & is.finite(end) & start <= end
  start <- as.numeric(start[ok])
  end <- as.numeric(end[ok])
  if (!length(start)) return(0)

  ord <- order(start, end)
  start <- start[ord]
  end <- end[ord]

  cur_start <- start[1L]
  cur_end <- end[1L]
  total <- 0

  if (length(start) >= 2L) {
    for (i in 2:length(start)) {
      if (start[i] <= cur_end + 1) {
        cur_end <- max(cur_end, end[i])
      } else {
        total <- total + (cur_end - cur_start + 1)
        cur_start <- start[i]
        cur_end <- end[i]
      }
    }
  }

  total + (cur_end - cur_start + 1)
}

find_final_consensus <- function(path) {
  if (!dir.exists(path)) stop("segments directory does not exist: ", path, call. = FALSE)

  ff <- list.files(
    path,
    pattern = "^segs_consensus_[0-9]+\\.tsv(\\.gz)?$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (!length(ff)) {
    stop("no segs_consensus_<iteration>.tsv[.gz] files found under ", path, call. = FALSE)
  }

  iter <- as.integer(
    sub(
      "^segs_consensus_([0-9]+)\\.tsv(?:\\.gz)?$",
      "\\1",
      basename(ff)
    )
  )

  if (anyNA(iter)) stop("failed to parse Numbat iteration numbers", call. = FALSE)

  max_iter <- max(iter)
  final <- ff[iter == max_iter]

  if (length(final) != 1L) {
    stop(
      "expected one consensus file at final iteration ", max_iter,
      "; found: ", paste(final, collapse = ", "),
      call. = FALSE
    )
  }

  list(file = final, iteration = max_iter)
}

if (!is.na(segments_file) && nzchar(segments_file)) {
  if (!file.exists(segments_file)) {
    stop("explicit --segments file does not exist: ", segments_file, call. = FALSE)
  }
  selected_segments <- segments_file
  bn <- basename(segments_file)
  selected_iteration <- if (
    grepl("^segs_consensus_[0-9]+\\.tsv(\\.gz)?$", bn)
  ) {
    as.integer(
      sub(
        "^segs_consensus_([0-9]+)\\.tsv(?:\\.gz)?$",
        "\\1",
        bn
      )
    )
  } else {
    NA_integer_
  }
} else {
  sel <- find_final_consensus(segs_dir)
  selected_segments <- sel$file
  selected_iteration <- sel$iteration
}

message("Numbat version: ", as.character(utils::packageVersion("numbat")))
message("selected consensus file: ", selected_segments)
message("NOTE: segs_consensus_<i> uses iteration number i; it is not a clone identifier.")

segs <- fread(selected_segments, showProgress = FALSE)

required_seg <- c("CHROM", "seg_start", "seg_end")
missing_seg <- setdiff(required_seg, names(segs))
if (length(missing_seg)) {
  stop(
    "consensus segments missing columns: ",
    paste(missing_seg, collapse = ", "),
    call. = FALSE
  )
}

if (!"cnv_state" %in% names(segs) && !"cnv_state_post" %in% names(segs)) {
  stop("consensus segments contain neither cnv_state nor cnv_state_post", call. = FALSE)
}
if (!"p_loh" %in% names(segs)) segs[, p_loh := NA_real_]
if (!"p_del" %in% names(segs)) segs[, p_del := NA_real_]

segs[, CHROM := normalize_chr(CHROM)]
segs[, seg_start := suppressWarnings(as.numeric(seg_start))]
segs[, seg_end := suppressWarnings(as.numeric(seg_end))]
segs[, p_loh := suppressWarnings(as.numeric(p_loh))]
segs[, p_del := suppressWarnings(as.numeric(p_del))]

segs <- segs[
  CHROM %in% as.character(1:22) &
    is.finite(seg_start) &
    is.finite(seg_end)
]

if (!nrow(segs)) stop("no valid autosomal consensus segments remain", call. = FALSE)
if (any(segs$seg_start < 1) || any(segs$seg_end < segs$seg_start)) {
  stop("invalid segment coordinates detected", call. = FALSE)
}

for (pp in c("p_loh", "p_del")) {
  finite <- is.finite(segs[[pp]])
  if (any(segs[[pp]][finite] < -1e-8 | segs[[pp]][finite] > 1 + 1e-8)) {
    stop(pp, " contains values outside [0,1]", call. = FALSE)
  }
}

state_original <- if ("cnv_state" %in% names(segs)) {
  tolower(trimws(as.character(segs$cnv_state)))
} else {
  rep(NA_character_, nrow(segs))
}

state_post <- if ("cnv_state_post" %in% names(segs)) {
  tolower(trimws(as.character(segs$cnv_state_post)))
} else {
  rep(NA_character_, nrow(segs))
}

state_original[is.na(state_original) | !nzchar(state_original)] <- NA_character_
state_post[is.na(state_post) | !nzchar(state_post)] <- NA_character_

segs[, state_final := fifelse(!is.na(state_post), state_post, state_original)]

valid_states <- c("neu", "del", "amp", "loh", "bamp", "bdel")
bad_states <- setdiff(unique(segs[!is.na(state_final), state_final]), valid_states)
if (length(bad_states)) {
  stop("unexpected Numbat states: ", paste(bad_states, collapse = ", "), call. = FALSE)
}

segs[
  ,
  is_cnloh := (
    state_final == "loh"
  ) | (
    is.finite(p_loh) &
      p_loh > posterior_threshold
  )
]

segs[
  ,
  p_loh_or_del := fifelse(
    is.finite(p_loh) | is.finite(p_del),
    fifelse(is.finite(p_loh), p_loh, 0) +
      fifelse(is.finite(p_del), p_del, 0),
    NA_real_
  )
]

if (any(is.finite(segs$p_loh_or_del) & segs$p_loh_or_del > 1 + 1e-6)) {
  stop("p_loh + p_del exceeds 1 for at least one segment", call. = FALSE)
}

segs[
  ,
  is_loh_or_del := (
    state_final %in% c("loh", "del")
  ) | (
    is.finite(p_loh_or_del) &
      p_loh_or_del > posterior_threshold
  )
]

state_p_disagreement <- segs[
  !is.na(state_final) &
    is.finite(p_loh) &
    (
      (state_final == "loh" & p_loh <= posterior_threshold) |
        (state_final != "loh" & p_loh > posterior_threshold)
    )
]

if (nrow(state_p_disagreement)) {
  warning(
    nrow(state_p_disagreement),
    " segments have MAP-state/p_loh-threshold disagreement.",
    call. = FALSE
  )
}

env <- new.env(parent = emptyenv())
utils::data("acen_hg38", package = "numbat", envir = env)
utils::data("chrom_sizes_hg38", package = "numbat", envir = env)

if (
  !exists("acen_hg38", envir = env, inherits = FALSE) ||
    !exists("chrom_sizes_hg38", envir = env, inherits = FALSE)
) {
  stop("could not load Numbat hg38 genome annotations", call. = FALSE)
}

acen <- as.data.table(get("acen_hg38", envir = env))
cs <- as.data.table(get("chrom_sizes_hg38", envir = env))

acen_chr_col <- pick_col(acen, c("CHROM", "chrom", "chr", "seqnames", "chromosome"), "acen_hg38")
acen_start_col <- pick_col(acen, c("START", "start", "chromStart", "acen_start"), "acen_hg38")
acen_end_col <- pick_col(acen, c("END", "end", "stop", "chromEnd", "acen_end"), "acen_hg38")

cs_chr_col <- pick_col(cs, c("CHROM", "chrom", "chr", "seqnames", "chromosome"), "chrom_sizes_hg38")
cs_size_col <- pick_col(cs, c("SIZE", "size", "length", "chrom_size", "chromLength", "end"), "chrom_sizes_hg38")

acen2 <- data.table(
  CHROM = normalize_chr(acen[[acen_chr_col]]),
  start = suppressWarnings(as.numeric(acen[[acen_start_col]])),
  end = suppressWarnings(as.numeric(acen[[acen_end_col]]))
)

acen2 <- acen2[
  CHROM %in% as.character(1:22) &
    is.finite(start) &
    is.finite(end)
]

acen2 <- acen2[
  ,
  .(
    centromere_start = min(start),
    centromere_end = max(end)
  ),
  by = CHROM
]

cs2 <- data.table(
  CHROM = normalize_chr(cs[[cs_chr_col]]),
  chrom_size = suppressWarnings(as.numeric(cs[[cs_size_col]]))
)

cs2 <- unique(
  cs2[
    CHROM %in% as.character(1:22) &
      is.finite(chrom_size)
  ],
  by = "CHROM"
)

hg38 <- merge(cs2, acen2, by = "CHROM", all = FALSE)

if (nrow(hg38) != 22L) {
  stop(
    "expected hg38 annotations for 22 autosomes; found ",
    nrow(hg38),
    call. = FALSE
  )
}

if (
  any(
    hg38$centromere_start < 1 |
      hg38$centromere_end < hg38$centromere_start |
      hg38$centromere_end >= hg38$chrom_size
  )
) {
  stop("invalid Numbat hg38 centromere annotation", call. = FALSE)
}

arms <- rbindlist(
  lapply(
    seq_len(nrow(hg38)),
    function(i) {
      chr <- hg38$CHROM[i]
      data.table(
        CHROM = chr,
        arm = c(paste0(chr, "p"), paste0(chr, "q")),
        arm_start = c(1, hg38$centromere_end[i] + 1),
        arm_end = c(hg38$centromere_start[i] - 1, hg38$chrom_size[i])
      )
    }
  )
)

arms[, arm_bp := arm_end - arm_start + 1]

if (any(!is.finite(arms$arm_bp) | arms$arm_bp <= 0)) {
  stop("invalid hg38 chromosome-arm intervals", call. = FALSE)
}

overlap_parts <- vector("list", nrow(arms))

for (i in seq_len(nrow(arms))) {
  a <- arms[i]

  x <- segs[
    CHROM == a$CHROM &
      seg_end >= a$arm_start &
      seg_start <= a$arm_end
  ]

  if (!nrow(x)) next

  x[, overlap_start := pmax(seg_start, a$arm_start)]
  x[, overlap_end := pmin(seg_end, a$arm_end)]
  x <- x[overlap_start <= overlap_end]

  if (!nrow(x)) next

  x[
    ,
    `:=`(
      arm = a$arm,
      arm_start = a$arm_start,
      arm_end = a$arm_end,
      arm_bp = a$arm_bp,
      overlap_bp = overlap_end - overlap_start + 1
    )
  ]

  overlap_parts[[i]] <- x
}

ov <- rbindlist(overlap_parts, use.names = TRUE, fill = TRUE)

if (!nrow(ov)) {
  stop("no final Numbat consensus segment overlaps an autosomal hg38 arm", call. = FALSE)
}

arm_calls <- rbindlist(
  lapply(
    seq_len(nrow(arms)),
    function(i) {
      a <- arms[i]
      x <- ov[arm == a$arm]

      if (!nrow(x)) {
        return(
          data.table(
            arm_id = a$arm,
            n_consensus_segments = 0L,
            n_cnloh_segments = 0L,
            n_loh_or_del_segments = 0L,
            any_cnloh_segment = FALSE,
            any_loh_or_del_segment = FALSE,
            max_p_loh = NA_real_,
            max_p_loh_or_del = NA_real_,
            cnloh_bp = 0,
            loh_or_del_bp = 0,
            arm_bp = a$arm_bp,
            cnloh_frac_arm = 0,
            loh_or_del_frac_arm = 0
          )
        )
      }

      cnloh_bp <- union_length(
        x[is_cnloh == TRUE, overlap_start],
        x[is_cnloh == TRUE, overlap_end]
      )

      loh_or_del_bp <- union_length(
        x[is_loh_or_del == TRUE, overlap_start],
        x[is_loh_or_del == TRUE, overlap_end]
      )

      data.table(
        arm_id = a$arm,
        n_consensus_segments = nrow(x),
        n_cnloh_segments = sum(x$is_cnloh, na.rm = TRUE),
        n_loh_or_del_segments = sum(x$is_loh_or_del, na.rm = TRUE),
        any_cnloh_segment = any(x$is_cnloh, na.rm = TRUE),
        any_loh_or_del_segment = any(x$is_loh_or_del, na.rm = TRUE),
        max_p_loh = safe_max(x$p_loh),
        max_p_loh_or_del = safe_max(x$p_loh_or_del),
        cnloh_bp = cnloh_bp,
        loh_or_del_bp = loh_or_del_bp,
        arm_bp = a$arm_bp,
        cnloh_frac_arm = cnloh_bp / a$arm_bp,
        loh_or_del_frac_arm = loh_or_del_bp / a$arm_bp
      )
    }
  )
)

arm_calls[, numbat_cnloh_call := cnloh_frac_arm > call_threshold]
arm_calls[, numbat_loh_or_del_call := loh_or_del_frac_arm > call_threshold]

truth <- fread(truth_file, showProgress = FALSE)

required_truth <- c("arm_id", "frac_loh", "masked")
missing_truth <- setdiff(required_truth, names(truth))
if (length(missing_truth)) {
  stop(
    "truth file missing columns: ",
    paste(missing_truth, collapse = ", "),
    call. = FALSE
  )
}

truth[, arm_id := normalize_arm_id(arm_id)]
truth[, frac_loh := suppressWarnings(as.numeric(frac_loh))]
truth[, masked_bool := as_logical_strict(masked, "truth$masked")]
truth[, chrom_from_arm := sub("[pq]$", "", arm_id)]

truth_eval <- truth[
  masked_bool == FALSE &
    chrom_from_arm %in% as.character(1:22) &
    is.finite(frac_loh)
]

if (!nrow(truth_eval)) {
  stop("no unmasked autosomal truth arms with finite frac_loh", call. = FALSE)
}

if (any(truth_eval$frac_loh < 0 | truth_eval$frac_loh > 1)) {
  stop("truth frac_loh contains values outside [0,1]", call. = FALSE)
}

if (anyDuplicated(truth_eval$arm_id)) {
  stop("truth contains duplicated evaluable arm_id rows", call. = FALSE)
}

truth_eval[, truth_loh := frac_loh > truth_threshold]

m <- merge(
  truth_eval,
  arm_calls,
  by = "arm_id",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(m$arm_bp)) {
  stop(
    "truth arms cannot be matched to hg38 arm definitions: ",
    paste(m[is.na(arm_bp), arm_id], collapse = ", "),
    call. = FALSE
  )
}

metric_block <- function(call, truth) {
  if (anyNA(call) || anyNA(truth)) stop("NA binary call/truth", call. = FALSE)

  tp <- sum(call & truth)
  fn <- sum(!call & truth)
  tn <- sum(!call & !truth)
  fp <- sum(call & !truth)

  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  ppv <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  npv <- if ((tn + fn) > 0) tn / (tn + fn) else NA_real_
  accuracy <- (tp + tn) / length(truth)

  data.table(
    TP = tp,
    FN = fn,
    TN = tn,
    FP = fp,
    sensitivity_pct = 100 * sensitivity,
    specificity_pct = 100 * specificity,
    PPV_pct = 100 * ppv,
    NPV_pct = 100 * npv,
    accuracy_pct = 100 * accuracy,
    balanced_accuracy_pct = if (
      is.finite(sensitivity) && is.finite(specificity)
    ) {
      100 * mean(c(sensitivity, specificity))
    } else {
      NA_real_
    }
  )
}

cnloh_metrics <- metric_block(m$numbat_cnloh_call, m$truth_loh)
cnloh_metrics[, endpoint := "CNLOH_ONLY"]

lohdel_metrics <- metric_block(m$numbat_loh_or_del_call, m$truth_loh)
lohdel_metrics[, endpoint := "CNLOH_OR_HEMIZYGOUS_DEL"]

metrics <- rbindlist(
  list(cnloh_metrics, lohdel_metrics),
  use.names = TRUE
)

setcolorder(metrics, c("endpoint", setdiff(names(metrics), "endpoint")))

cat("\n================ NUMBAT C1 LOH CONCORDANCE vs SEQC2 ARM TRUTH ================\n")
cat("Numbat version: ", as.character(utils::packageVersion("numbat")), "\n", sep = "")
cat("final consensus file: ", selected_segments, "\n", sep = "")
cat(
  "selected Numbat iteration: ",
  ifelse(is.na(selected_iteration), "NA", selected_iteration),
  "\n",
  sep = ""
)
cat(sprintf("truth definition: frac_loh > %.3f\n", truth_threshold))
cat(
  sprintf(
    "Numbat arm call: qualifying LOH coverage / non-centromeric hg38 arm bp > %.3f\n",
    call_threshold
  )
)
cat(sprintf("segment posterior threshold: > %.3f\n", posterior_threshold))
cat(sprintf("unmasked autosomal arms compared: %d\n", nrow(m)))
cat(sprintf("truth LOH arms: %d\n", sum(m$truth_loh)))
cat(sprintf("truth non-LOH arms: %d\n\n", sum(!m$truth_loh)))

cat("=== PRIMARY: Numbat copy-neutral LOH state only ===\n")
print(metrics[endpoint == "CNLOH_ONLY"])

cat("\n=== SECONDARY: CNLoH OR hemizygous deletion ===\n")
print(metrics[endpoint == "CNLOH_OR_HEMIZYGOUS_DEL"])

if ("17q" %in% m$arm_id) {
  z <- m[arm_id == "17q"]

  cat("\n=== 17q ===\n")
  cat(
    sprintf(
      paste0(
        "truth frac_loh=%.4f; truth_LOH=%s; ",
        "CNLoH fraction=%.4f; CNLoH call=%s; ",
        "LOH-or-del fraction=%.4f; LOH-or-del call=%s; ",
        "max p_loh=%s\n"
      ),
      z$frac_loh,
      z$truth_loh,
      z$cnloh_frac_arm,
      z$numbat_cnloh_call,
      z$loh_or_del_frac_arm,
      z$numbat_loh_or_del_call,
      ifelse(
        is.finite(z$max_p_loh),
        sprintf("%.4f", z$max_p_loh),
        "NA"
      )
    )
  )
}

cat("\n=== per-arm detail: truth LOH arms ===\n")
print(
  m[
    truth_loh == TRUE
  ][
    order(-frac_loh),
    .(
      arm_id,
      frac_loh,
      numbat_cnloh_call,
      cnloh_frac_arm,
      max_p_loh,
      numbat_loh_or_del_call,
      loh_or_del_frac_arm,
      max_p_loh_or_del,
      n_consensus_segments
    )
  ]
)

per_arm_path <- file.path(outdir, "numbat_loh_concordance_per_arm.tsv")
metrics_path <- file.path(outdir, "numbat_loh_concordance_metrics.tsv")
segments_path <- file.path(outdir, "numbat_loh_final_segments_scored.tsv")
overlap_path <- file.path(outdir, "numbat_loh_segment_arm_overlaps.tsv")
metadata_path <- file.path(outdir, "numbat_loh_concordance_metadata.tsv")

fwrite(m, per_arm_path, sep = "\t")
fwrite(metrics, metrics_path, sep = "\t")
fwrite(segs, segments_path, sep = "\t")
fwrite(ov, overlap_path, sep = "\t")

metadata <- data.table(
  numbat_version = as.character(utils::packageVersion("numbat")),
  selected_consensus_file = normalizePath(selected_segments, mustWork = TRUE),
  selected_iteration = selected_iteration,
  truth_file = normalizePath(truth_file, mustWork = TRUE),
  truth_threshold = truth_threshold,
  call_threshold = call_threshold,
  posterior_threshold = posterior_threshold,
  arm_annotation = "numbat::acen_hg38 + numbat::chrom_sizes_hg38",
  final_state_rule = "cnv_state_post if present/nonmissing, otherwise cnv_state",
  cnloh_segment_rule = "state_final==loh OR p_loh>posterior_threshold",
  loh_or_del_segment_rule = "state_final in {loh,del} OR (p_loh+p_del)>posterior_threshold",
  arm_call_rule = "union qualifying bp / non-centromeric hg38 arm bp > call_threshold",
  truth_filter = "masked==FALSE; autosomes; finite frac_loh; state_call not used",
  consensus_scope = "final Numbat iteration only"
)

fwrite(metadata, metadata_path, sep = "\t")

cat(
  "\nwrote:\n",
  "  ", per_arm_path, "\n",
  "  ", metrics_path, "\n",
  "  ", segments_path, "\n",
  "  ", overlap_path, "\n",
  "  ", metadata_path, "\n",
  sep = ""
)
