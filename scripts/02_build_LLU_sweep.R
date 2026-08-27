library(Matrix)

# ============================================================
# LLU in-silico tumour purity sweep
# Lead-approved design:
#   - Cell Ranger 3.1
#   - LLU_B split 30% reference / 70% mixture pool
#   - N = 1000 cells per mixture
#   - Tumour fractions = 1, 5, 10, 20, 40, 80%
#   - 10 replicates per fraction
#   - no replacement within any mixture
#   - reference cells never appear in mixtures
#   - preserve source barcode + replicate
# ============================================================

library(Matrix)

BASE_SEED <- 20260826
fractions <- c(0.01, 0.05, 0.10, 0.20, 0.40, 0.80)
n_total <- 1000
n_reps <- 10

dir.create("/work/metadata", showWarnings = FALSE, recursive = TRUE)
dir.create("/work/manifests", showWarnings = FALSE, recursive = TRUE)
dir.create("/work/logs", showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Load prepared matrices
# ----------------------------

A <- readRDS("/work/data/LLU_A_cellranger3.1_aligned.rds")
B <- readRDS("/work/data/LLU_B_cellranger3.1_aligned.rds")

A_original <- colnames(A)
B_original <- colnames(B)

# Prefix barcodes to guarantee uniqueness between A and B
colnames(A) <- paste0("LLU_A__", A_original)
colnames(B) <- paste0("LLU_B__", B_original)

n_B <- ncol(B)

# ----------------------------
# Depth summary
# ----------------------------

A_umi <- Matrix::colSums(A)
B_umi <- Matrix::colSums(B)

depth_summary <- data.frame(
  sample = c("LLU_A", "LLU_B"),
  total_cells = c(length(A_umi), length(B_umi)),
  cells_lt_2000_UMI = c(sum(A_umi < 2000), sum(B_umi < 2000)),
  cells_ge_2000_UMI = c(sum(A_umi >= 2000), sum(B_umi >= 2000)),
  median_UMI = c(median(A_umi), median(B_umi))
)

write.table(
  depth_summary,
  "/work/logs/LLU_depth_2000UMI.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# NOTE:
# 2000 UMI is logged as a CopyKAT stability diagnostic,
# not imposed as a hard filtering threshold.

# ----------------------------
# Barcode map
# ----------------------------

barcode_map <- rbind(
  data.frame(
    cell_id = colnames(A),
    original_barcode = A_original,
    source = "LLU_A",
    truth = "tumor"
  ),
  data.frame(
    cell_id = colnames(B),
    original_barcode = B_original,
    source = "LLU_B",
    truth = "normal"
  )
)

write.table(
  barcode_map,
  "/work/metadata/LLU_barcode_map.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------
# Fixed 30:70 LLU_B split
# ----------------------------

set.seed(BASE_SEED)

n_ref <- round(0.30 * n_B)

ref_idx <- sample(
  seq_len(n_B),
  size = n_ref,
  replace = FALSE
)

B_ref_cells <- colnames(B)[ref_idx]
B_mix_cells <- setdiff(colnames(B), B_ref_cells)

stopifnot(
  length(intersect(B_ref_cells, B_mix_cells)) == 0,
  length(B_ref_cells) + length(B_mix_cells) == n_B
)

split_meta <- data.frame(
  cell_id = colnames(B),
  original_barcode = sub("^LLU_B__", "", colnames(B)),
  source = "LLU_B",
  role = ifelse(
    colnames(B) %in% B_ref_cells,
    "reference",
    "mixture_pool"
  )
)

write.table(
  split_meta,
  "/work/metadata/LLU_B_split.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

writeLines(
  B_ref_cells,
  "/work/metadata/LLU_B_reference_cells.txt"
)

writeLines(
  B_mix_cells,
  "/work/metadata/LLU_B_mixture_pool_cells.txt"
)

# Save prefixed matrices for caller scripts
saveRDS(A, "/work/data/LLU_A_prefixed.rds")
saveRDS(B, "/work/data/LLU_B_prefixed.rds")

# ----------------------------
# Sweep design table
# ----------------------------

design <- data.frame(
  tumor_fraction = fractions,
  tumor_cells = round(n_total * fractions),
  normal_cells = n_total - round(n_total * fractions)
)

stopifnot(
  all(design$tumor_cells <= ncol(A)),
  all(design$normal_cells <= length(B_mix_cells))
)

write.table(
  design,
  "/work/metadata/LLU_sweep_design.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------
# Generate manifests
# ----------------------------

all_manifests <- list()
manifest_counter <- 1

for (f_idx in seq_along(fractions)) {

  f <- fractions[f_idx]
  n_tumor <- round(n_total * f)
  n_normal <- n_total - n_tumor

  # Where possible, tumour cells are kept disjoint across the
  # 10 replicates for this purity level.
  tumor_disjoint_possible <- (n_tumor * n_reps <= ncol(A))

  if (tumor_disjoint_possible) {
    set.seed(BASE_SEED + f_idx * 10000)
    tumor_order <- sample(colnames(A), ncol(A), replace = FALSE)
  }

  for (rep_id in seq_len(n_reps)) {

    dataset_id <- sprintf(
      "LLU_p%02d_r%02d",
      round(f * 100),
      rep_id
    )

    if (tumor_disjoint_possible) {

      start_i <- (rep_id - 1) * n_tumor + 1
      end_i <- rep_id * n_tumor

      tumor_sel <- tumor_order[start_i:end_i]

    } else {

      # No replacement within a replicate.
      # Across-replicate overlap is permitted when the requested
      # tumour-cell total exceeds the available A-cell pool.
      set.seed(BASE_SEED + f_idx * 10000 + rep_id)

      tumor_sel <- sample(
        colnames(A),
        size = n_tumor,
        replace = FALSE
      )
    }

    # Normal cells may recur across replicates, but never duplicate
    # within a replicate and never overlap with held-out reference.
    set.seed(BASE_SEED + 500000 + f_idx * 10000 + rep_id)

    normal_sel <- sample(
      B_mix_cells,
      size = n_normal,
      replace = FALSE
    )

    manifest <- rbind(
      data.frame(
        dataset_id = dataset_id,
        cell_id = tumor_sel,
        original_barcode = sub("^LLU_A__", "", tumor_sel),
        source = "LLU_A",
        truth = "tumor",
        tumor_fraction = f,
        replicate = rep_id
      ),
      data.frame(
        dataset_id = dataset_id,
        cell_id = normal_sel,
        original_barcode = sub("^LLU_B__", "", normal_sel),
        source = "LLU_B",
        truth = "normal",
        tumor_fraction = f,
        replicate = rep_id
      )
    )

    stopifnot(
      nrow(manifest) == n_total,
      sum(duplicated(manifest$cell_id)) == 0,
      sum(manifest$truth == "tumor") == n_tumor,
      sum(manifest$truth == "normal") == n_normal,
      length(intersect(manifest$cell_id, B_ref_cells)) == 0
    )

    write.table(
      manifest,
      paste0("/work/manifests/", dataset_id, ".tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    all_manifests[[manifest_counter]] <- manifest
    manifest_counter <- manifest_counter + 1
  }
}

master_manifest <- do.call(rbind, all_manifests)

write.table(
  master_manifest,
  "/work/manifests/LLU_all_mixtures.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------
# Reproducibility / validation summary
# ----------------------------

sweep_summary <- data.frame(
  tumor_fraction = fractions,
  n_replicates = n_reps,
  cells_per_replicate = n_total,
  tumor_cells_per_replicate = round(n_total * fractions),
  normal_cells_per_replicate = n_total - round(n_total * fractions),
  tumor_disjoint_across_replicates =
    round(n_total * fractions) * n_reps <= ncol(A)
)

write.table(
  sweep_summary,
  "/work/logs/LLU_sweep_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("LLU sweep generation completed successfully.\n")
cat("LLU_B total:", n_B, "\n")
cat("LLU_B reference:", length(B_ref_cells), "\n")
cat("LLU_B mixture pool:", length(B_mix_cells), "\n")
cat("Mixtures generated:", length(all_manifests), "\n")
cat("Rows in master manifest:", nrow(master_manifest), "\n")
cat("Reference/mixture overlap:",
    length(intersect(B_ref_cells, B_mix_cells)), "\n")
cat("A <2000 UMI:", sum(A_umi < 2000), "\n")
cat("B <2000 UMI:", sum(B_umi < 2000), "\n")
