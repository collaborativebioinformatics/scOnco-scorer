library(Matrix)
library(copykat)

dataset_id <- "LLU_p80_r01"

dir.create("/work/results/copykat", recursive = TRUE, showWarnings = FALSE)

# Load prefixed matrices
A <- readRDS("/work/data/LLU_A_prefixed.rds")
B <- readRDS("/work/data/LLU_B_prefixed.rds")

# Load mixture manifest
manifest <- read.delim(
  paste0("/work/manifests/", dataset_id, ".tsv"),
  stringsAsFactors = FALSE
)

# Load permanently held-out LLU_B reference cells
ref_cells <- readLines("/work/metadata/LLU_B_reference_cells.txt")

tumor_cells <- manifest$cell_id[manifest$source == "LLU_A"]
normal_mix_cells <- manifest$cell_id[manifest$source == "LLU_B"]
mixture_cells <- manifest$cell_id

# Validate IDs
stopifnot(
  length(mixture_cells) == 1000,
  length(tumor_cells) == 800,
  length(normal_mix_cells) == 200,
  length(ref_cells) == 432,
  length(intersect(ref_cells, mixture_cells)) == 0,
  all(tumor_cells %in% colnames(A)),
  all(normal_mix_cells %in% colnames(B)),
  all(ref_cells %in% colnames(B))
)

# CopyKAT input:
# 1000 mixture cells + 432 held-out normal reference cells
rawmat <- cbind(
  A[, tumor_cells, drop = FALSE],
  B[, normal_mix_cells, drop = FALSE],
  B[, ref_cells, drop = FALSE]
)

rawmat <- as.matrix(rawmat)


stopifnot(
  ncol(rawmat) == 1432,
  !anyDuplicated(colnames(rawmat))
)

cat("Dataset:", dataset_id, "\n")
cat("Tumor mixture cells:", length(tumor_cells), "\n")
cat("Normal mixture cells:", length(normal_mix_cells), "\n")
cat("Held-out reference cells:", length(ref_cells), "\n")
cat("CopyKAT input cells:", ncol(rawmat), "\n")
cat("Input genes:", nrow(rawmat), "\n")
cat("Reference/mixture overlap:",
    length(intersect(ref_cells, mixture_cells)), "\n")

set.seed(20260826)

ck <- copykat(
  rawmat = rawmat,
  id.type = "E",
  cell.line = "no",
  norm.cell.names = ref_cells,
  sam.name = dataset_id,
  genome = "hg20",
  n.cores = 16
)

saveRDS(
  ck,
  paste0(
    "/work/results/copykat/",
    dataset_id,
    "_copykat_full.rds"
  )
)

cat("CopyKAT pilot completed successfully.\n")
cat("Returned object names:\n")
print(names(ck))
