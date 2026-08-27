library(Matrix)
library(infercnv)

dataset_id <- "LLU_p80_r01"

counts_A <- readRDS("/work/data/LLU_A_prefixed.rds")
counts_B <- readRDS("/work/data/LLU_B_prefixed.rds")

manifest <- read.delim(
  "/work/manifests/LLU_p80_r01.tsv",
  stringsAsFactors = FALSE
)

B_ref <- readLines("/work/metadata/LLU_B_reference_cells.txt")

# Assemble exact caller input:
# held-out B reference + 3000 mixture cells
test_cells <- manifest$cell_id
input_cells <- c(B_ref, test_cells)

counts <- cbind(
  counts_B[, B_ref, drop = FALSE],
  counts_A[, intersect(test_cells, colnames(counts_A)), drop = FALSE],
  counts_B[, intersect(test_cells, colnames(counts_B)), drop = FALSE]
)

stopifnot(
  ncol(counts) == length(input_cells),
  all(B_ref %in% colnames(counts)),
  all(test_cells %in% colnames(counts))
)

annotations <- data.frame(
  cell = colnames(counts),
  group = ifelse(
    colnames(counts) %in% B_ref,
    "normal_ref",
    "observed"
  )
)

write.table(
  annotations,
  file = paste0("/work/results/infercnv/", dataset_id, "_annotations.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

gene_order_file <- "/work/reference/GRCh38_gene_order.tsv"

if (!file.exists(gene_order_file)) {
  stop(
    "Missing GRCh38 gene-order file: ",
    gene_order_file,
    ". Obtain the team-approved pinned resource before running."
  )
}

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts,
  gene_order_file = gene_order_file,
  annotations_file = paste0(
    "/work/results/infercnv/",
    dataset_id,
    "_annotations.tsv"
  ),
  ref_group_names = "normal_ref"
)

saveRDS(
  infercnv_obj,
  paste0(
    "/work/results/infercnv/",
    dataset_id,
    "_input_object.rds"
  )
)

cat("inferCNV pilot input prepared successfully\n")
