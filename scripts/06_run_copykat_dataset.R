library(Matrix)
library(copykat)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript 06_run_copykat_dataset.R <dataset_id>")
}

dataset_id <- args[1]

cat("========================================\n")
cat("CopyKAT dataset:", dataset_id, "\n")
cat("Start:", format(Sys.time()), "\n")
cat("========================================\n")

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

out_root <- file.path("/work/results/copykat", dataset_id)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

full_rds <- file.path(out_root, paste0(dataset_id, "_copykat_full.rds"))
pred_tsv <- file.path(out_root, paste0(dataset_id, "_predictions.tsv"))
input_tsv <- file.path(out_root, paste0(dataset_id, "_input_summary.tsv"))
done_file <- file.path(out_root, "DONE")

# Resumable behavior
if (file.exists(done_file) &&
    file.exists(full_rds) &&
    file.exists(pred_tsv)) {
  cat("SKIP:", dataset_id, "already completed.\n")
  quit(save = "no", status = 0)
}

# ------------------------------------------------------------
# Load prepared expression matrices
# ------------------------------------------------------------

A <- readRDS("/work/data/LLU_A_prefixed.rds")
B <- readRDS("/work/data/LLU_B_prefixed.rds")

# ------------------------------------------------------------
# Pull this dataset from master manifest
# ------------------------------------------------------------

master <- read.delim(
  "/work/manifests/LLU_all_mixtures.tsv",
  stringsAsFactors = FALSE
)

manifest <- master[master$dataset_id == dataset_id, , drop = FALSE]

if (nrow(manifest) == 0) {
  stop("Dataset not found in master manifest: ", dataset_id)
}

if (nrow(manifest) != 1000) {
  stop("Expected 1000 mixture cells, found ", nrow(manifest))
}

ref_cells <- readLines(
  "/work/metadata/LLU_B_reference_cells.txt"
)

tumor_cells <- manifest$cell_id[manifest$source == "LLU_A"]
normal_mix_cells <- manifest$cell_id[manifest$source == "LLU_B"]

# ------------------------------------------------------------
# Input validation
# ------------------------------------------------------------

stopifnot(
  length(ref_cells) == 432,
  length(intersect(ref_cells, manifest$cell_id)) == 0,
  !anyDuplicated(manifest$cell_id),
  all(tumor_cells %in% colnames(A)),
  all(normal_mix_cells %in% colnames(B)),
  all(ref_cells %in% colnames(B))
)

expected_tumor <- unique(round(1000 * manifest$tumor_fraction))
stopifnot(length(expected_tumor) == 1)
stopifnot(length(tumor_cells) == expected_tumor)
stopifnot(length(normal_mix_cells) == 1000 - expected_tumor)

# ------------------------------------------------------------
# Build CopyKAT input
# ------------------------------------------------------------

rawmat <- cbind(
  A[, tumor_cells, drop = FALSE],
  B[, normal_mix_cells, drop = FALSE],
  B[, ref_cells, drop = FALSE]
)

stopifnot(
  ncol(rawmat) == 1432,
  !anyDuplicated(colnames(rawmat))
)

# CopyKAT 1.2.5 annotation functions require a base matrix.
rawmat <- as.matrix(rawmat)

input_summary <- data.frame(
  dataset_id = dataset_id,
  tumor_fraction = unique(manifest$tumor_fraction),
  replicate = unique(manifest$replicate),
  mixture_total = nrow(manifest),
  tumor_cells = length(tumor_cells),
  mixture_normal_cells = length(normal_mix_cells),
  heldout_reference_cells = length(ref_cells),
  total_copykat_input = ncol(rawmat),
  genes = nrow(rawmat),
  reference_mixture_overlap =
    length(intersect(ref_cells, manifest$cell_id))
)

write.table(
  input_summary,
  input_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Tumor fraction:", unique(manifest$tumor_fraction), "\n")
cat("Replicate:", unique(manifest$replicate), "\n")
cat("Tumor cells:", length(tumor_cells), "\n")
cat("Mixture normal cells:", length(normal_mix_cells), "\n")
cat("Held-out reference:", length(ref_cells), "\n")
cat("Total CopyKAT input:", ncol(rawmat), "\n")
cat("Genes:", nrow(rawmat), "\n")

# ------------------------------------------------------------
# Run CopyKAT
# ------------------------------------------------------------

start_time <- Sys.time()

ck <- copykat(
  rawmat = rawmat,
  id.type = "E",
  cell.line = "no",
  norm.cell.names = ref_cells,
  sam.name = dataset_id,
  genome = "hg20",
  n.cores = 16
)

end_time <- Sys.time()

saveRDS(
  ck,
  full_rds
)

# ------------------------------------------------------------
# Build compact per-cell output
#
# Keep ALL 1000 mixture cells in the table.
# Cells removed by CopyKAT before prediction remain explicit.
# ------------------------------------------------------------

pred <- ck$prediction

compact <- merge(
  manifest,
  pred[, c("cell.names", "copykat.pred")],
  by.x = "cell_id",
  by.y = "cell.names",
  all.x = TRUE,
  sort = FALSE
)

# Preserve manifest ordering
compact <- compact[
  match(manifest$cell_id, compact$cell_id),
  ,
  drop = FALSE
]

compact$copykat_status <- ifelse(
  is.na(compact$copykat.pred),
  "filtered_before_prediction",
  ifelse(
    compact$copykat.pred == "not.defined",
    "not_defined",
    "called"
  )
)

# Binary label only for definitive CopyKAT calls.
# This is NOT a calibrated probability.
compact$copykat_binary <- NA_integer_

compact$copykat_binary[
  compact$copykat.pred == "aneuploid"
] <- 1L

compact$copykat_binary[
  compact$copykat.pred == "diploid"
] <- 0L

write.table(
  compact,
  pred_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# Reference-cell diagnostic
# ------------------------------------------------------------

ref_pred <- pred[pred$cell.names %in% ref_cells, , drop = FALSE]

ref_summary <- data.frame(
  dataset_id = dataset_id,
  reference_supplied = length(ref_cells),
  reference_predicted = nrow(ref_pred),
  reference_aneuploid =
    sum(ref_pred$copykat.pred == "aneuploid"),
  reference_diploid =
    sum(ref_pred$copykat.pred == "diploid"),
  reference_not_defined =
    sum(ref_pred$copykat.pred == "not.defined"),
  reference_absent =
    sum(!ref_cells %in% pred$cell.names)
)

write.table(
  ref_summary,
  file.path(out_root, paste0(dataset_id, "_reference_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Run summary
# ------------------------------------------------------------

run_summary <- data.frame(
  dataset_id = dataset_id,
  tumor_fraction = unique(manifest$tumor_fraction),
  replicate = unique(manifest$replicate),
  runtime_minutes =
    as.numeric(difftime(end_time, start_time, units = "mins")),
  mixture_cells = nrow(manifest),
  mixture_predicted =
    sum(manifest$cell_id %in% pred$cell.names),
  mixture_filtered =
    sum(!manifest$cell_id %in% pred$cell.names),
  aneuploid =
    sum(compact$copykat.pred == "aneuploid", na.rm = TRUE),
  diploid =
    sum(compact$copykat.pred == "diploid", na.rm = TRUE),
  not_defined =
    sum(compact$copykat.pred == "not.defined", na.rm = TRUE)
)

write.table(
  run_summary,
  file.path(out_root, paste0(dataset_id, "_run_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

writeLines(
  paste(
    "completed",
    dataset_id,
    format(Sys.time()),
    sep = "\t"
  ),
  done_file
)

cat("\nCopyKAT completed:", dataset_id, "\n")
print(run_summary)
cat("Output:", out_root, "\n")
