#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (!length(hit)) return(default)
  i <- hit[length(hit)]
  if (i == length(args)) stop(flag, " requires a value", call. = FALSE)
  args[i + 1L]
}

pred_root <- get_arg("--pred_root", "/w/infercnv_out")
manifest_file <- get_arg("--manifest", "/w/inputs/LLU_all_mixtures_with_p00.tsv")
outdir <- get_arg("--outdir", pred_root)
burden_col <- get_arg("--burden_col", "infercnv_burden_masked")

if (!dir.exists(pred_root)) stop("prediction root does not exist: ", pred_root, call. = FALSE)
if (!file.exists(manifest_file)) stop("manifest does not exist: ", manifest_file, call. = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

files <- sort(list.files(
  pred_root,
  pattern = "_infercnv_burden\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
))

if (!length(files)) {
  stop("no *_infercnv_burden.tsv files found under ", pred_root, call. = FALSE)
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

  z <- data.table(score = score, y = y)[
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

  sum(diff(c(0, z$recall)) * z$precision)
}

read_burden <- function(path) {
  x <- fread(path, showProgress = FALSE)

  required <- c("cell_id", "dataset_id", "in_mixture", burden_col)
  missing <- setdiff(required, names(x))

  if (length(missing)) {
    stop(
      basename(path),
      " missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  x[, source_burden_file := normalizePath(path, mustWork = TRUE)]
  x
}

pred <- rbindlist(
  lapply(files, read_burden),
  use.names = TRUE,
  fill = TRUE
)

pred[, cell_id := as.character(cell_id)]
pred[, dataset_id := as.character(dataset_id)]
pred[, score := suppressWarnings(as.numeric(get(burden_col)))]

if (is.logical(pred$in_mixture)) {
  pred[, in_mixture_flag := in_mixture]
} else {
  z <- tolower(trimws(as.character(pred$in_mixture)))
  pred[, in_mixture_flag := z %in% c("true", "t", "1", "yes", "y")]
}

if (anyNA(pred$cell_id) || any(!nzchar(pred$cell_id))) {
  stop("inferCNV burden files contain NA/empty cell_id", call. = FALSE)
}
if (anyNA(pred$dataset_id) || any(!nzchar(pred$dataset_id))) {
  stop("inferCNV burden files contain NA/empty dataset_id", call. = FALSE)
}

pred_key <- paste(pred$dataset_id, pred$cell_id, sep = "||")

if (anyDuplicated(pred_key)) {
  dup <- pred[
    duplicated(pred_key) | duplicated(pred_key, fromLast = TRUE),
    .(dataset_id, cell_id)
  ]
  stop(
    "duplicated dataset_id/cell_id burden rows detected; first: ",
    dup$dataset_id[1L],
    "/",
    dup$cell_id[1L],
    call. = FALSE
  )
}

manifest <- fread(manifest_file, showProgress = FALSE)

required_manifest <- c(
  "dataset_id",
  "cell_id",
  "truth",
  "tumor_fraction",
  "replicate"
)

missing_manifest <- setdiff(required_manifest, names(manifest))

if (length(missing_manifest)) {
  stop(
    "manifest missing required columns: ",
    paste(missing_manifest, collapse = ", "),
    call. = FALSE
  )
}

manifest[, dataset_id := as.character(dataset_id)]
manifest[, cell_id := as.character(cell_id)]
manifest[, tumor_fraction := suppressWarnings(as.numeric(tumor_fraction))]
manifest[, replicate := as.character(replicate)]
manifest[, y := normalize_truth(truth)]

if (
  anyNA(manifest$dataset_id) ||
  any(!nzchar(manifest$dataset_id)) ||
  anyNA(manifest$cell_id) ||
  any(!nzchar(manifest$cell_id))
) {
  stop("manifest contains NA/empty dataset_id or cell_id", call. = FALSE)
}

if (
  anyNA(manifest$tumor_fraction) ||
  any(manifest$tumor_fraction < 0 | manifest$tumor_fraction > 1)
) {
  stop("manifest tumor_fraction must be numeric in [0,1]", call. = FALSE)
}

manifest_key <- paste(manifest$dataset_id, manifest$cell_id, sep = "||")

if (anyDuplicated(manifest_key)) {
  stop("manifest contains duplicated dataset_id/cell_id rows", call. = FALSE)
}

available_ds <- unique(pred$dataset_id)
manifest_sub <- manifest[dataset_id %in% available_ds]

if (!nrow(manifest_sub)) {
  stop("none of the inferCNV dataset IDs occur in the manifest", call. = FALSE)
}

mix_key <- paste(manifest_sub$dataset_id, manifest_sub$cell_id, sep = "||")

pred_mixture <- pred[in_mixture_flag == TRUE]
pred_mixture_key <- paste(
  pred_mixture$dataset_id,
  pred_mixture$cell_id,
  sep = "||"
)

missing_predictions <- manifest_sub[!mix_key %in% pred_mixture_key]

if (nrow(missing_predictions)) {
  stop(
    nrow(missing_predictions),
    " manifest mixture cells have no inferCNV mixture burden row; first: ",
    missing_predictions$dataset_id[1L],
    "/",
    missing_predictions$cell_id[1L],
    call. = FALSE
  )
}

unexpected_mixture <- pred_mixture[!pred_mixture_key %in% mix_key]

if (nrow(unexpected_mixture)) {
  stop(
    nrow(unexpected_mixture),
    " inferCNV rows are flagged in_mixture=TRUE but absent from manifest; first: ",
    unexpected_mixture$dataset_id[1L],
    "/",
    unexpected_mixture$cell_id[1L],
    call. = FALSE
  )
}

merged <- merge(
  manifest_sub,
  pred_mixture[
    ,
    .(
      dataset_id,
      cell_id,
      score,
      source_burden_file
    )
  ],
  by = c("dataset_id", "cell_id"),
  all.x = TRUE,
  sort = FALSE
)

if (nrow(merged) != nrow(manifest_sub)) {
  stop("manifest/inferCNV burden merge is not one-to-one", call. = FALSE)
}

merged[, infercnv_scored := is.finite(score)]

dataset_qc <- merged[
  ,
  .(
    n_manifest = .N,
    n_scored = sum(infercnv_scored),
    n_unscored = sum(!infercnv_scored),
    unscored_fraction = mean(!infercnv_scored),
    n_tumor_total = sum(y == 1L),
    n_normal_total = sum(y == 0L),
    n_tumor_scored = sum(y == 1L & infercnv_scored),
    n_normal_scored = sum(y == 0L & infercnv_scored)
  ),
  by = .(dataset_id, tumor_fraction, replicate)
][order(tumor_fraction, replicate, dataset_id)]

if (any(dataset_qc$n_unscored > 0L)) {
  bad <- dataset_qc[n_unscored > 0L]
  stop(
    "non-finite/missing ",
    burden_col,
    " values detected in inferCNV mixture cells. ",
    "The corrected inferCNV runner should produce a finite primary burden for every mixture cell. ",
    "Affected dataset(s): ",
    paste(
      paste0(
        bad$dataset_id,
        " [",
        bad$n_unscored,
        "/",
        bad$n_manifest,
        " unscored]"
      ),
      collapse = "; "
    ),
    call. = FALSE
  )
}

per_replicate <- merged[
  ,
  {
    q <- .SD[infercnv_scored == TRUE]

    list(
      n_total = .N,
      n_scored = nrow(q),
      n_unscored = .N - nrow(q),
      unscored_fraction = (.N - nrow(q)) / .N,
      n_tumor = sum(q$y == 1L),
      n_normal = sum(q$y == 0L),
      prevalence = if (nrow(q)) mean(q$y) else NA_real_,
      AUROC = auroc(q$score, q$y),
      AUPRC_AP = average_precision(q$score, q$y),
      median_burden_tumor = if (sum(q$y == 1L) > 0L) {
        median(q$score[q$y == 1L])
      } else {
        NA_real_
      },
      median_burden_normal = if (sum(q$y == 0L) > 0L) {
        median(q$score[q$y == 0L])
      } else {
        NA_real_
      },
      median_burden_delta = if (
        sum(q$y == 1L) > 0L &&
        sum(q$y == 0L) > 0L
      ) {
        median(q$score[q$y == 1L]) -
          median(q$score[q$y == 0L])
      } else {
        NA_real_
      }
    )
  },
  by = .(dataset_id, tumor_fraction, replicate)
][order(tumor_fraction, replicate, dataset_id)]

summary_by_fraction <- per_replicate[
  ,
  .(
    n_replicates = .N,
    mean_AUROC = if (all(is.na(AUROC))) NA_real_ else mean(AUROC, na.rm = TRUE),
    median_AUROC = if (all(is.na(AUROC))) NA_real_ else median(AUROC, na.rm = TRUE),
    min_AUROC = if (all(is.na(AUROC))) NA_real_ else min(AUROC, na.rm = TRUE),
    max_AUROC = if (all(is.na(AUROC))) NA_real_ else max(AUROC, na.rm = TRUE),
    mean_AUPRC_AP = if (all(is.na(AUPRC_AP))) NA_real_ else mean(AUPRC_AP, na.rm = TRUE),
    median_AUPRC_AP = if (all(is.na(AUPRC_AP))) NA_real_ else median(AUPRC_AP, na.rm = TRUE),
    mean_prevalence = if (all(is.na(prevalence))) NA_real_ else mean(prevalence, na.rm = TRUE),
    mean_unscored_fraction = mean(unscored_fraction),
    median_burden_delta = if (all(is.na(median_burden_delta))) {
      NA_real_
    } else {
      median(median_burden_delta, na.rm = TRUE)
    }
  ),
  by = tumor_fraction
][order(tumor_fraction)]

merged[, burden_column_used := burden_col]

if (all(c("source", "original_barcode") %in% names(merged))) {
  merged[
    ,
    source_cell_id := paste(
      as.character(source),
      as.character(original_barcode),
      sep = "||"
    )
  ]

  source_truth_check <- merged[
    ,
    .(n_truth = uniqueN(y)),
    by = source_cell_id
  ]

  if (any(source_truth_check$n_truth != 1L)) {
    stop(
      "source/original_barcode grouping produces inconsistent truth labels",
      call. = FALSE
    )
  }

  message(
    "source-cell identifier retained as source_cell_id = source || original_barcode"
  )
} else {
  warning(
    paste0(
      "manifest lacks source and/or original_barcode; ",
      "per-cell scoring is valid, but the final source-cell cluster bootstrap ",
      "cannot be reconstructed from this output alone"
    ),
    call. = FALSE
  )
}

message("burden column used: ", burden_col)
message("burden files: ", length(files))
message("mixture rows: ", nrow(merged))
message("datasets: ", uniqueN(merged$dataset_id))

print(per_replicate)
cat("\n")
print(summary_by_fraction)

percell_path <- file.path(outdir, "infercnv_burden_percell_all.tsv")
perrep_path <- file.path(outdir, "infercnv_burden_metrics_per_replicate.tsv")
summary_path <- file.path(outdir, "infercnv_burden_summary_by_fraction.tsv")
qc_path <- file.path(outdir, "infercnv_burden_dataset_qc.tsv")

fwrite(merged, percell_path, sep = "\t")
fwrite(per_replicate, perrep_path, sep = "\t")
fwrite(summary_by_fraction, summary_path, sep = "\t")
fwrite(dataset_qc, qc_path, sep = "\t")

cat(
  "\nwrote:\n",
  "  ", percell_path, "\n",
  "  ", perrep_path, "\n",
  "  ", summary_path, "\n",
  "  ", qc_path, "\n",
  sep = ""
)
