library(Matrix)

set.seed(20260826)

A <- readRDS("/work/data/LLU_A_cellranger3.1_aligned.rds")
B <- readRDS("/work/data/LLU_B_cellranger3.1_aligned.rds")

# Preserve original barcodes before prefixing
A_original <- colnames(A)
B_original <- colnames(B)

# Prefix cell IDs so A and B can never collide
colnames(A) <- paste0("LLU_A__", A_original)
colnames(B) <- paste0("LLU_B__", B_original)

# Save barcode map
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

# Fixed 50:50 split of LLU_B
n_B <- ncol(B)

n_ref <- floor(n_B / 2)

ref_idx <- sample(
  seq_len(n_B),
  size = n_ref,
  replace = FALSE
)

B_ref_cells <- colnames(B)[ref_idx]
B_mix_cells <- setdiff(colnames(B), B_ref_cells)

split_meta <- data.frame(
  cell_id = colnames(B),
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

# Save matrices with unique prefixed cell IDs
saveRDS(A, "/work/data/LLU_A_prefixed.rds")
saveRDS(B, "/work/data/LLU_B_prefixed.rds")

# ------------------------------------------------------------
# 80% tumor pilot: 3000 cells total
# 2400 tumor + 600 normal
# This is feasible without replacement from the held-out B pool.
# ------------------------------------------------------------

set.seed(20260826 + 800 + 1)

tumor_sel <- sample(
  colnames(A),
  size = 2400,
  replace = FALSE
)

normal_sel <- sample(
  B_mix_cells,
  size = 600,
  replace = FALSE
)

pilot_manifest <- rbind(
  data.frame(
    dataset_id = "LLU_p80_r01",
    cell_id = tumor_sel,
    source = "LLU_A",
    truth = "tumor",
    tumor_fraction = 0.80,
    replicate = 1
  ),
  data.frame(
    dataset_id = "LLU_p80_r01",
    cell_id = normal_sel,
    source = "LLU_B",
    truth = "normal",
    tumor_fraction = 0.80,
    replicate = 1
  )
)

write.table(
  pilot_manifest,
  "/work/manifests/LLU_p80_r01.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("LLU_B total:", n_B, "\n")
cat("LLU_B reference:", length(B_ref_cells), "\n")
cat("LLU_B mixture pool:", length(B_mix_cells), "\n")
cat("80% pilot tumor cells:", sum(pilot_manifest$truth == "tumor"), "\n")
cat("80% pilot normal cells:", sum(pilot_manifest$truth == "normal"), "\n")
cat("80% pilot total:", nrow(pilot_manifest), "\n")
cat("Duplicate pilot cell IDs:", sum(duplicated(pilot_manifest$cell_id)), "\n")
