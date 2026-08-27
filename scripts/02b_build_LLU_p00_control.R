# LLU B-vs-B negative control for CopyKAT
# 1000 scored LLU_B mixture-pool cells + 432 disjoint held-out LLU_B reference cells

BASE_SEED <- 20260827
dataset_id <- "LLU_p00_r01"

B <- readRDS("/work/data/LLU_B_prefixed.rds")

ref_cells <- readLines(
  "/work/metadata/LLU_B_reference_cells.txt"
)

mix_pool <- readLines(
  "/work/metadata/LLU_B_mixture_pool_cells.txt"
)

stopifnot(
  length(ref_cells) == 432,
  length(mix_pool) == 1007,
  length(intersect(ref_cells, mix_pool)) == 0,
  all(ref_cells %in% colnames(B)),
  all(mix_pool %in% colnames(B))
)

set.seed(BASE_SEED)

scored_B <- sample(
  mix_pool,
  size = 1000,
  replace = FALSE
)

manifest <- data.frame(
  dataset_id = dataset_id,
  cell_id = scored_B,
  original_barcode = sub("^LLU_B__", "", scored_B),
  source = "LLU_B",
  truth = "normal",
  tumor_fraction = 0,
  replicate = 1,
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(manifest) == 1000,
  sum(duplicated(manifest$cell_id)) == 0,
  length(intersect(manifest$cell_id, ref_cells)) == 0
)

write.table(
  manifest,
  "/work/manifests/LLU_p00_r01.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Append p00 locally to the master manifest used by the generic runner.
master_file <- "/work/manifests/LLU_all_mixtures.tsv"
master <- read.delim(master_file, stringsAsFactors = FALSE)

master <- master[master$dataset_id != dataset_id, , drop = FALSE]
master <- rbind(master, manifest)

write.table(
  master,
  master_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Built", dataset_id, "\n")
cat("Scored B cells:", nrow(manifest), "\n")
cat("Reference B cells:", length(ref_cells), "\n")
cat("Scored/reference overlap:",
    length(intersect(manifest$cell_id, ref_cells)), "\n")
cat("Master manifest rows:", nrow(master), "\n")
