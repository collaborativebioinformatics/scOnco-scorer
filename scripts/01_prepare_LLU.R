library(Matrix)

input_file <- "/work/data/gene_counts_all.rdata"

load(input_file)

LLU_A_raw <- gene_counts[["10X_LLU_A_cellranger3.1"]]
LLU_B_raw <- gene_counts[["10X_LLU_B_cellranger3.1"]]

# Force gene IDs from dimnames onto rownames
rownames(LLU_A_raw) <- dimnames(LLU_A_raw)[[1]]
rownames(LLU_B_raw) <- dimnames(LLU_B_raw)[[1]]

# Strip Ensembl version suffixes if present
rownames(LLU_A_raw) <- sub("\\.[0-9]+$", "", rownames(LLU_A_raw))
rownames(LLU_B_raw) <- sub("\\.[0-9]+$", "", rownames(LLU_B_raw))

# Log ID issues before duplicate collapse
id_log <- data.frame(
  sample = c("LLU_A", "LLU_B"),
  original_genes = c(nrow(LLU_A_raw), nrow(LLU_B_raw)),
  missing_or_blank_ids = c(
    sum(is.na(rownames(LLU_A_raw)) | rownames(LLU_A_raw) == ""),
    sum(is.na(rownames(LLU_B_raw)) | rownames(LLU_B_raw) == "")
  ),
  duplicated_gene_ids = c(
    sum(duplicated(rownames(LLU_A_raw))),
    sum(duplicated(rownames(LLU_B_raw)))
  )
)

write.table(
  id_log,
  "/work/logs/LLU_gene_id_log.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Drop missing/blank IDs if any
keep_A <- !(is.na(rownames(LLU_A_raw)) | rownames(LLU_A_raw) == "")
keep_B <- !(is.na(rownames(LLU_B_raw)) | rownames(LLU_B_raw) == "")

LLU_A_raw <- LLU_A_raw[keep_A, , drop = FALSE]
LLU_B_raw <- LLU_B_raw[keep_B, , drop = FALSE]

# Sparse-safe duplicate gene collapse
collapse_duplicate_rows <- function(mat) {
  ids <- rownames(mat)

  if (!anyDuplicated(ids)) {
    return(mat)
  }

  idx <- match(ids, unique(ids))

  G <- sparseMatrix(
    i = idx,
    j = seq_along(idx),
    x = 1,
    dims = c(length(unique(ids)), length(ids))
  )

  out <- G %*% mat
  rownames(out) <- unique(ids)
  colnames(out) <- colnames(mat)

  out
}

LLU_A_raw <- collapse_duplicate_rows(LLU_A_raw)
LLU_B_raw <- collapse_duplicate_rows(LLU_B_raw)

# QC metrics
A_umi   <- Matrix::colSums(LLU_A_raw)
B_umi   <- Matrix::colSums(LLU_B_raw)

A_genes <- Matrix::colSums(LLU_A_raw > 0)
B_genes <- Matrix::colSums(LLU_B_raw > 0)

qc_summary <- data.frame(
  sample = c("LLU_A", "LLU_B"),
  genes = c(nrow(LLU_A_raw), nrow(LLU_B_raw)),
  cells = c(ncol(LLU_A_raw), ncol(LLU_B_raw)),
  median_UMI = c(median(A_umi), median(B_umi)),
  median_genes = c(median(A_genes), median(B_genes)),
  UMI_lt_1000 = c(sum(A_umi < 1000), sum(B_umi < 1000)),
  genes_lt_200 = c(sum(A_genes < 200), sum(B_genes < 200)),
  versioned_gene_IDs = c(
    sum(grepl("\\.[0-9]+$", rownames(LLU_A_raw))),
    sum(grepl("\\.[0-9]+$", rownames(LLU_B_raw)))
  ),
  duplicate_gene_IDs = c(
    sum(duplicated(rownames(LLU_A_raw))),
    sum(duplicated(rownames(LLU_B_raw)))
  ),
  duplicate_cell_IDs = c(
    sum(duplicated(colnames(LLU_A_raw))),
    sum(duplicated(colnames(LLU_B_raw)))
  )
)

write.table(
  qc_summary,
  "/work/logs/LLU_cellranger3.1_QC.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Align A and B to union gene universe with sparse zero-fill
all_genes <- union(rownames(LLU_A_raw), rownames(LLU_B_raw))

align_to_genes <- function(mat, genes) {

  missing_genes <- setdiff(genes, rownames(mat))

  if (length(missing_genes) > 0) {

    zeros <- Matrix(
      0,
      nrow = length(missing_genes),
      ncol = ncol(mat),
      sparse = TRUE,
      dimnames = list(missing_genes, colnames(mat))
    )

    mat <- rbind(mat, zeros)
  }

  mat[genes, , drop = FALSE]
}

LLU_A <- align_to_genes(LLU_A_raw, all_genes)
LLU_B <- align_to_genes(LLU_B_raw, all_genes)

# Alignment summary
A_only <- setdiff(rownames(LLU_A_raw), rownames(LLU_B_raw))
B_only <- setdiff(rownames(LLU_B_raw), rownames(LLU_A_raw))

gene_alignment_summary <- data.frame(
  metric = c(
    "A_original_genes",
    "B_original_genes",
    "shared_genes",
    "A_only_genes",
    "B_only_genes",
    "union_genes"
  ),
  value = c(
    nrow(LLU_A_raw),
    nrow(LLU_B_raw),
    length(intersect(rownames(LLU_A_raw), rownames(LLU_B_raw))),
    length(A_only),
    length(B_only),
    length(all_genes)
  )
)

write.table(
  gene_alignment_summary,
  "/work/logs/LLU_gene_alignment.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Validation
stopifnot(
  identical(rownames(LLU_A), rownames(LLU_B)),
  sum(Matrix::colSums(LLU_A)) == sum(Matrix::colSums(LLU_A_raw)),
  sum(Matrix::colSums(LLU_B)) == sum(Matrix::colSums(LLU_B_raw))
)

saveRDS(
  LLU_A,
  "/work/data/LLU_A_cellranger3.1_aligned.rds"
)

saveRDS(
  LLU_B,
  "/work/data/LLU_B_cellranger3.1_aligned.rds"
)

cat("LLU preparation completed successfully.\n")
cat("LLU_A:", nrow(LLU_A), "genes x", ncol(LLU_A), "cells\n")
cat("LLU_B:", nrow(LLU_B), "genes x", ncol(LLU_B), "cells\n")
