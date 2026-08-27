library(Matrix)
library(SCEVAN)

dataset_id <- "LLU_p80_r01"

counts_A <- readRDS("/work/data/LLU_A_prefixed.rds")
counts_B <- readRDS("/work/data/LLU_B_prefixed.rds")

manifest <- read.delim(
  "/work/manifests/LLU_p80_r01.tsv",
  stringsAsFactors = FALSE
)

B_ref <- readLines("/work/metadata/LLU_B_reference_cells.txt")
test_cells <- manifest$cell_id

counts <- cbind(
  counts_B[, B_ref, drop = FALSE],
  counts_A[, intersect(test_cells, colnames(counts_A)), drop = FALSE],
  counts_B[, intersect(test_cells, colnames(counts_B)), drop = FALSE]
)

stopifnot(
  all(B_ref %in% colnames(counts)),
  all(test_cells %in% colnames(counts))
)

wrapper_file <- "/work/code/scevan_wrapper_pinned.R"

if (!file.exists(wrapper_file)) {
  stop(
    "Missing pinned SCEVAN wrapper: ",
    wrapper_file,
    ". Do not run default pipelineCNA because it discards CNAmat."
  )
}

source(wrapper_file)

cat("SCEVAN pilot input prepared; pinned wrapper loaded\n")

# Actual wrapper function call will be added once the team provides
# the pinned wrapper and confirms its function name/signature.
