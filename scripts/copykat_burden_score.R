#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (!length(hit)) return(default)
  if (hit[length(hit)] == length(args)) stop(flag, " requires a value", call. = FALSE)
  args[hit[length(hit)] + 1L]
}

runs_dir <- get_arg("--runs_dir", "/w/runs")
mask_file <- get_arg("--mask", "/w/mask_regions.tsv")
outdir <- get_arg("--outdir", "/w")

if (!dir.exists(runs_dir)) stop("runs directory does not exist: ", runs_dir, call. = FALSE)
if (!file.exists(mask_file)) stop("mask file does not exist: ", mask_file, call. = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

rds <- sort(list.files(runs_dir, pattern = "_full\\.rds$", full.names = TRUE))
if (!length(rds)) stop("no *_full.rds files found in ", runs_dir, call. = FALSE)

mask <- fread(mask_file, showProgress = FALSE)
required_mask <- c("CHROM", "START", "END")
missing_mask <- setdiff(required_mask, names(mask))
if (length(missing_mask)) {
  stop("mask file missing columns: ", paste(missing_mask, collapse = ", "), call. = FALSE)
}

mask[, CHROM := sub("^chr", "", as.character(CHROM), ignore.case = TRUE)]
mask[, START := suppressWarnings(as.numeric(START))]
mask[, END := suppressWarnings(as.numeric(END))]

if (
  anyNA(mask$CHROM) ||
  anyNA(mask$START) ||
  anyNA(mask$END) ||
  any(mask$START < 1) ||
  any(mask$END < mask$START)
) {
  stop("invalid mask coordinates", call. = FALSE)
}

auroc <- function(score, y) {
  ok <- is.finite(score) & !is.na(y)
  score <- score[ok]
  y <- y[ok]
  if (length(unique(y)) != 2L) return(NA_real_)
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

average_precision <- function(score, y) {
  ok <- is.finite(score) & !is.na(y)
  score <- score[ok]
  y <- y[ok]
  if (length(unique(y)) != 2L) return(NA_real_)
  n_positive <- sum(y == 1L)
  if (n_positive == 0L) return(NA_real_)

  z <- data.table(score = score, y = y)
  z <- z[
    ,
    .(
      positives = sum(y == 1L),
      negatives = sum(y == 0L)
    ),
    by = score
  ][order(-score)]

  z[, tp := cumsum(positives)]
  z[, fp := cumsum(negatives)]
  z[, recall := tp / n_positive]
  z[, precision := tp / (tp + fp)]

  delta_recall <- diff(c(0, z$recall))
  sum(delta_recall * z$precision)
}

normalize_truth <- function(x) {
  z <- tolower(trimws(as.character(x)))
  y <- rep(NA_integer_, length(z))
  y[z %in% c("tumor", "tumour", "a", "hcc1395")] <- 1L
  y[z %in% c("normal", "b", "hcc1395bl")] <- 0L

  if (anyNA(y)) {
    stop(
      "unrecognized truth labels: ",
      paste(unique(z[is.na(y)]), collapse = ", "),
      call. = FALSE
    )
  }
  y
}

mask_rows <- function(cna) {
  required_meta <- c("chrom", "chrompos", "abspos")
  missing_meta <- setdiff(required_meta, colnames(cna))

  if (length(missing_meta)) {
    stop(
      "CNAmat missing metadata columns: ",
      paste(missing_meta, collapse = ", "),
      call. = FALSE
    )
  }

  chr <- sub("^chr", "", as.character(cna$chrom), ignore.case = TRUE)
  pos <- suppressWarnings(as.numeric(cna$chrompos))

  if (anyNA(chr) || anyNA(pos)) {
    stop("invalid chrom/chrompos values in CNAmat", call. = FALSE)
  }

  masked <- rep(FALSE, nrow(cna))

  for (i in seq_len(nrow(mask))) {
    masked <- masked | (
      chr == mask$CHROM[i] &
      pos >= mask$START[i] &
      pos <= mask$END[i]
    )
  }

  masked
}

one_run <- function(f) {
  run <- sub("_full\\.rds$", "", basename(f))
  pred_file <- file.path(runs_dir, paste0(run, "_pred.tsv"))

  if (!file.exists(pred_file)) {
    stop("missing prediction file for ", run, ": ", pred_file, call. = FALSE)
  }

  x <- readRDS(f)

  if (!is.list(x) || is.null(x$CNAmat)) {
    stop(run, ": RDS does not contain CNAmat", call. = FALSE)
  }

  cna <- x$CNAmat

  if (!(is.matrix(cna) || is.data.frame(cna))) {
    stop(run, ": CNAmat is not a matrix/data.frame", call. = FALSE)
  }

  required_meta <- c("chrom", "chrompos", "abspos")
  missing_meta <- setdiff(required_meta, colnames(cna))

  if (length(missing_meta)) {
    stop(
      run,
      ": missing CNAmat columns: ",
      paste(missing_meta, collapse = ", "),
      call. = FALSE
    )
  }

  cells <- setdiff(colnames(cna), required_meta)

  if (!length(cells)) stop(run, ": CNAmat contains no cell columns", call. = FALSE)
  if (anyDuplicated(cells)) stop(run, ": duplicated CNAmat cell names", call. = FALSE)

  masked <- mask_rows(cna)
  if (all(masked)) stop(run, ": mask removes every CNAmat row", call. = FALSE)

  cna_use <- cna[!masked, cells, drop = FALSE]

  M <- suppressWarnings(
    matrix(
      as.numeric(as.matrix(cna_use)),
      nrow = nrow(cna_use),
      ncol = ncol(cna_use),
      dimnames = list(NULL, cells)
    )
  )

  if (any(!is.finite(M))) {
    stop(run, ": non-finite CNA values detected", call. = FALSE)
  }

  burden <- colMeans(abs(M))

  burden_dt <- data.table(
    cell = cells,
    burden = as.numeric(burden)
  )

  pred <- fread(pred_file, showProgress = FALSE)

  required_pred <- c("cell_id", "truth", "tumor_fraction", "replicate")
  missing_pred <- setdiff(required_pred, names(pred))

  if (length(missing_pred)) {
    stop(
      run,
      ": prediction file missing columns: ",
      paste(missing_pred, collapse = ", "),
      call. = FALSE
    )
  }

  pred[, cell_id := as.character(cell_id)]

  if (anyNA(pred$cell_id) || any(!nzchar(pred$cell_id))) {
    stop(run, ": NA/empty cell IDs in prediction file", call. = FALSE)
  }

  if (anyDuplicated(pred$cell_id)) {
    stop(run, ": duplicated cell_id entries in prediction file", call. = FALSE)
  }

  pred[, y := normalize_truth(truth)]

  dt <- merge(
    pred[
      ,
      .(
        cell = cell_id,
        truth,
        y,
        tumor_fraction,
        replicate
      )
    ],
    burden_dt,
    by = "cell",
    all.x = TRUE,
    sort = FALSE
  )

  if (nrow(dt) != nrow(pred)) {
    stop(run, ": one-to-one manifest merge failed", call. = FALSE)
  }

  dt[, copykat_scored := is.finite(burden)]
  dt[, run_id := run]
  dt[, n_bins_total := nrow(cna)]
  dt[, n_bins_masked := sum(masked)]
  dt[, n_bins_used := sum(!masked)]

  message(sprintf(
    "%s: manifest=%d scored=%d filtered=%d bins=%d masked=%d used=%d",
    run,
    nrow(dt),
    sum(dt$copykat_scored),
    sum(!dt$copykat_scored),
    nrow(cna),
    sum(masked),
    sum(!masked)
  ))

  dt[]
}

all <- rbindlist(
  lapply(rds, one_run),
  use.names = TRUE,
  fill = FALSE
)

if (!nrow(all)) stop("no per-cell results produced", call. = FALSE)

per_rep <- all[
  ,
  {
    q <- .SD[copykat_scored == TRUE]

    list(
      n_total = .N,
      n_scored = nrow(q),
      n_filtered = .N - nrow(q),
      filtered_fraction = (.N - nrow(q)) / .N,
      n_tumor = sum(q$y == 1L),
      n_normal = sum(q$y == 0L),
      AUROC = auroc(q$burden, q$y),
      AUPRC_AP = average_precision(q$burden, q$y),
      median_burden_tumor = if (sum(q$y == 1L)) median(q$burden[q$y == 1L]) else NA_real_,
      median_burden_normal = if (sum(q$y == 0L)) median(q$burden[q$y == 0L]) else NA_real_
    )
  },
  by = .(tumor_fraction, replicate)
][order(tumor_fraction, replicate)]

by_fraction <- per_rep[
  ,
  .(
    n_replicates = .N,
    mean_AUROC = if (all(is.na(AUROC))) NA_real_ else mean(AUROC, na.rm = TRUE),
    median_AUROC = if (all(is.na(AUROC))) NA_real_ else median(AUROC, na.rm = TRUE),
    mean_AUPRC_AP = if (all(is.na(AUPRC_AP))) NA_real_ else mean(AUPRC_AP, na.rm = TRUE),
    median_AUPRC_AP = if (all(is.na(AUPRC_AP))) NA_real_ else median(AUPRC_AP, na.rm = TRUE),
    mean_filtered_fraction = mean(filtered_fraction)
  ),
  by = tumor_fraction
][order(tumor_fraction)]

print(per_rep)
cat("\n")
print(by_fraction)

percell_path <- file.path(outdir, "copykat_burden_percell.tsv")
perrep_path <- file.path(outdir, "copykat_burden_metrics_per_replicate.tsv")
summary_path <- file.path(outdir, "copykat_burden_summary_by_fraction.tsv")

fwrite(all, percell_path, sep = "\t")
fwrite(per_rep, perrep_path, sep = "\t")
fwrite(by_fraction, summary_path, sep = "\t")

cat(
  "\nwrote:\n",
  "  ", percell_path, "\n",
  "  ", perrep_path, "\n",
  "  ", summary_path, "\n",
  sep = ""
)
