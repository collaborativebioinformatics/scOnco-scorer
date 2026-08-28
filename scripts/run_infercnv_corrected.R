#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(infercnv)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (!length(hit)) return(default)
  i <- hit[length(hit)]
  if (i == length(args)) stop(flag, " requires a value", call. = FALSE)
  args[i + 1L]
}

if (!length(args)) {
  stop(
    paste0(
      "usage: Rscript run_infercnv.R <dataset_id> ",
      "[--inputs /w/inputs] [--outroot /w/infercnv_out] ",
      "[--gene_order /w/inputs/gene_order.txt] ",
      "[--mask /w/mask_regions.tsv] [--cores 8]"
    ),
    call. = FALSE
  )
}

ds <- args[1L]

if (!grepl("^[A-Za-z0-9_.-]+$", ds)) {
  stop("unsafe dataset_id: ", ds, call. = FALSE)
}

inputs <- get_arg("--inputs", "/w/inputs")
outroot <- get_arg("--outroot", "/w/infercnv_out")
gene_order_file <- get_arg("--gene_order", file.path(inputs, "gene_order.txt"))
mask_file <- get_arg("--mask", NA_character_)
cores <- as.integer(get_arg("--cores", "8"))
min_gene_overlap <- as.numeric(get_arg("--min_gene_overlap", "0.90"))

if (is.na(cores) || cores < 1L) stop("--cores must be >=1", call. = FALSE)
if (
  !is.finite(min_gene_overlap) ||
  min_gene_overlap <= 0 ||
  min_gene_overlap > 1
) {
  stop("--min_gene_overlap must be in (0,1]", call. = FALSE)
}

A_file <- file.path(inputs, "LLU_A_prefixed.rds")
B_file <- file.path(inputs, "LLU_B_prefixed.rds")
manifest_file <- file.path(inputs, "LLU_all_mixtures_with_p00.tsv")
reference_file <- file.path(inputs, "LLU_B_reference_cells.txt")

for (p in c(
  A_file,
  B_file,
  manifest_file,
  reference_file,
  gene_order_file
)) {
  if (!file.exists(p)) stop("missing required input: ", p, call. = FALSE)
  if (file.info(p)$size <= 0) stop("empty required input: ", p, call. = FALSE)
}

outdir <- file.path(outroot, ds)

if (dir.exists(outdir) && file.exists(file.path(outdir, "DONE"))) {
  stop("completed inferCNV output already exists: ", outdir, call. = FALSE)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("=== inferCNV ", ds, " ===")
message("infercnv version: ", as.character(packageVersion("infercnv")))
message("cores: ", cores)
message("output: ", outdir)

collapse_duplicate_rows <- function(x, ids) {
  if (length(ids) != nrow(x)) {
    stop("row-ID length mismatch", call. = FALSE)
  }

  u <- unique(ids)
  map <- match(ids, u)

  S <- Matrix::sparseMatrix(
    i = map,
    j = seq_along(ids),
    x = 1,
    dims = c(length(u), length(ids)),
    dimnames = list(u, NULL)
  )

  y <- S %*% x
  colnames(y) <- colnames(x)
  as(y, "dgCMatrix")
}

prepare_counts <- function(x, label) {
  if (!(is.matrix(x) || inherits(x, "Matrix"))) {
    stop(label, " must be a matrix/Matrix object", call. = FALSE)
  }

  x <- as(Matrix::Matrix(x, sparse = TRUE), "dgCMatrix")

  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop(label, " must have rownames and colnames", call. = FALSE)
  }

  if (
    anyNA(rownames(x)) ||
    any(!nzchar(rownames(x))) ||
    anyNA(colnames(x)) ||
    any(!nzchar(colnames(x)))
  ) {
    stop(label, " has NA/empty identifiers", call. = FALSE)
  }

  if (anyDuplicated(colnames(x))) {
    stop(label, " has duplicated cell IDs", call. = FALSE)
  }

  if (
    any(!is.finite(x@x)) ||
    any(x@x < 0) ||
    any(abs(x@x - round(x@x)) > 1e-8)
  ) {
    stop(
      label,
      " must contain finite, nonnegative integer raw counts",
      call. = FALSE
    )
  }

  ids <- sub("\\.\\d+$", "", rownames(x))
  ndup <- length(ids) - length(unique(ids))

  if (ndup > 0L) {
    message(label, ": collapsing ", ndup, " duplicate rows after Ensembl version stripping")
    x <- collapse_duplicate_rows(x, ids)
  } else {
    rownames(x) <- ids
  }

  if (anyDuplicated(rownames(x))) {
    stop(label, ": duplicate genes remain after collapse", call. = FALSE)
  }

  if (any(Matrix::colSums(x) <= 0)) {
    bad <- colnames(x)[Matrix::colSums(x) <= 0]
    stop(
      label,
      ": zero-count cells present; first: ",
      paste(head(bad, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  x
}

prepare_gene_order <- function(path) {
  go <- fread(
    path,
    header = FALSE,
    sep = "\t",
    fill = TRUE,
    quote = "",
    showProgress = FALSE
  )

  if (ncol(go) < 4L || nrow(go) == 0L) {
    stop(
      "gene_order must contain at least four columns: gene, chr, start, end",
      call. = FALSE
    )
  }

  go <- go[, 1:4]
  setnames(go, c("gene", "chr", "start", "stop"))

  go[, gene := sub("\\.\\d+$", "", as.character(gene))]
  go[, chr := as.character(chr)]
  go[, chr := sub("^chr", "", chr, ignore.case = TRUE)]
  go[, chr := ifelse(chr %in% c("M", "MT"), "M", chr)]
  go[, chr := paste0("chr", chr)]
  go[, start := suppressWarnings(as.numeric(start))]
  go[, stop := suppressWarnings(as.numeric(stop))]

  go <- go[
    !is.na(gene) &
      nzchar(gene) &
      !is.na(chr) &
      nzchar(chr) &
      !is.na(start) &
      !is.na(stop)
  ]

  if (!nrow(go)) stop("no valid rows remain in gene_order", call. = FALSE)

  if (any(go$start < 1) || any(go$stop < go$start)) {
    stop("gene_order contains invalid coordinates", call. = FALSE)
  }

  cross_chr <- go[
    ,
    .(n_chr = uniqueN(chr)),
    by = gene
  ][n_chr > 1L]

  if (nrow(cross_chr)) {
    stop(
      "version-stripped gene_order has genes mapping to multiple chromosomes; first: ",
      paste(head(cross_chr$gene, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  go <- go[
    ,
    .(
      chr = first(chr),
      start = min(start),
      stop = max(stop)
    ),
    by = gene
  ]

  canonical <- c(
    paste0("chr", 1:22),
    "chrX",
    "chrY",
    "chrM"
  )

  go <- go[chr %in% canonical]

  if (!nrow(go)) {
    stop("gene_order has no canonical GRCh38 chromosomes after normalization", call. = FALSE)
  }

  chr_levels <- canonical[canonical %in% unique(go$chr)]
  go[, chr_rank := match(chr, chr_levels)]
  setorder(go, chr_rank, start, stop, gene)
  go[, chr_rank := NULL]

  out <- as.data.frame(go[, .(chr, start, stop)])
  rownames(out) <- go$gene

  if (anyDuplicated(rownames(out))) {
    stop("gene_order still has duplicated genes", call. = FALSE)
  }

  out
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

read_mask <- function(path) {
  if (is.na(path) || !nzchar(path)) return(NULL)

  if (!file.exists(path)) {
    stop("--mask file does not exist: ", path, call. = FALSE)
  }

  m <- fread(path, showProgress = FALSE)

  req <- c("CHROM", "START", "END")
  miss <- setdiff(req, names(m))

  if (length(miss)) {
    stop(
      "--mask missing columns: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }

  m[, CHROM := sub("^chr", "", as.character(CHROM), ignore.case = TRUE)]
  m[, CHROM := ifelse(CHROM %in% c("M", "MT"), "M", CHROM)]
  m[, CHROM := paste0("chr", CHROM)]
  m[, START := suppressWarnings(as.numeric(START))]
  m[, END := suppressWarnings(as.numeric(END))]

  if (
    anyNA(m$CHROM) ||
    anyNA(m$START) ||
    anyNA(m$END) ||
    any(m$START < 1) ||
    any(m$END < m$START)
  ) {
    stop("invalid mask coordinates", call. = FALSE)
  }

  m
}

auroc <- function(score, y) {
  ok <- is.finite(score) & !is.na(y)
  score <- score[ok]
  y <- y[ok]

  if (length(unique(y)) != 2L) return(NA_real_)

  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)

  if (!n1 || !n0) return(NA_real_)

  r <- rank(score, ties.method = "average")

  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

A <- prepare_counts(
  readRDS(A_file),
  "LLU_A_prefixed.rds"
)

B <- prepare_counts(
  readRDS(B_file),
  "LLU_B_prefixed.rds"
)

if (length(intersect(colnames(A), colnames(B))) > 0L) {
  stop("A and B matrices contain overlapping cell IDs", call. = FALSE)
}

common_AB <- intersect(rownames(A), rownames(B))

if (!length(common_AB)) {
  stop("A and B matrices have no common genes", call. = FALSE)
}

frac_A <- length(common_AB) / nrow(A)
frac_B <- length(common_AB) / nrow(B)

if (frac_A < 0.90 || frac_B < 0.90) {
  stop(
    sprintf(
      "A/B common-gene overlap too low: A %.2f%%, B %.2f%%",
      100 * frac_A,
      100 * frac_B
    ),
    call. = FALSE
  )
}

A <- A[common_AB, , drop = FALSE]
B <- B[common_AB, , drop = FALSE]

gene_order <- prepare_gene_order(gene_order_file)

shared_GO <- intersect(common_AB, rownames(gene_order))
gene_overlap <- length(shared_GO) / length(common_AB)

message(
  sprintf(
    "count/gene-order overlap: %d/%d genes (%.2f%%)",
    length(shared_GO),
    length(common_AB),
    100 * gene_overlap
  )
)

if (gene_overlap < min_gene_overlap) {
  stop(
    sprintf(
      "count/gene-order overlap %.3f is below threshold %.3f",
      gene_overlap,
      min_gene_overlap
    ),
    call. = FALSE
  )
}

A <- A[shared_GO, , drop = FALSE]
B <- B[shared_GO, , drop = FALSE]
gene_order <- gene_order[shared_GO, , drop = FALSE]

man <- fread(manifest_file, showProgress = FALSE)

required_manifest <- c(
  "dataset_id",
  "cell_id",
  "truth",
  "tumor_fraction",
  "replicate"
)

missing_manifest <- setdiff(required_manifest, names(man))

if (length(missing_manifest)) {
  stop(
    "manifest missing columns: ",
    paste(missing_manifest, collapse = ", "),
    call. = FALSE
  )
}

sub <- man[dataset_id == ds]

if (!nrow(sub)) {
  stop("dataset_id not found in manifest: ", ds, call. = FALSE)
}

sub[, cell_id := as.character(cell_id)]
sub[, tumor_fraction := suppressWarnings(as.numeric(tumor_fraction))]
sub[, replicate := as.character(replicate)]
sub[, y := normalize_truth(truth)]

if (
  anyNA(sub$cell_id) ||
  any(!nzchar(sub$cell_id)) ||
  anyDuplicated(sub$cell_id)
) {
  stop(ds, ": invalid or duplicated mixture cell IDs", call. = FALSE)
}

if (
  anyNA(sub$tumor_fraction) ||
  length(unique(sub$tumor_fraction)) != 1L
) {
  stop(ds, ": manifest must have one valid tumor_fraction", call. = FALSE)
}

if (length(unique(sub$replicate)) != 1L) {
  stop(ds, ": manifest must have one replicate value", call. = FALSE)
}

mix <- sub$cell_id

refcells <- trimws(
  readLines(reference_file, warn = FALSE)
)
refcells <- refcells[nzchar(refcells)]

if (!length(refcells)) {
  stop("reference-B list is empty", call. = FALSE)
}

if (anyDuplicated(refcells)) {
  stop("reference-B list contains duplicates", call. = FALSE)
}

if (length(intersect(mix, refcells)) > 0L) {
  stop(
    ds,
    ": mixture/reference overlap detected",
    call. = FALSE
  )
}

all_cells <- c(colnames(A), colnames(B))

missing_mix <- setdiff(mix, all_cells)

if (length(missing_mix)) {
  stop(
    ds,
    ": mixture cells absent from A/B matrices; first: ",
    paste(head(missing_mix, 10L), collapse = ", "),
    call. = FALSE
  )
}

missing_ref <- setdiff(refcells, colnames(B))

if (length(missing_ref)) {
  stop(
    ds,
    ": reference-B cells absent from B matrix; first: ",
    paste(head(missing_ref, 10L), collapse = ", "),
    call. = FALSE
  )
}

if (length(intersect(refcells, colnames(A))) > 0L) {
  stop("reference list contains A cells", call. = FALSE)
}

mix_in_A <- mix[mix %in% colnames(A)]
mix_in_B <- mix[mix %in% colnames(B)]

if (length(mix_in_A) + length(mix_in_B) != length(mix)) {
  stop(ds, ": inconsistent source assignment for mixture cells", call. = FALSE)
}

comb <- cbind(A, B)
use_cells <- c(mix, refcells)

if (anyDuplicated(use_cells)) {
  stop("duplicate cells after mixture + reference construction", call. = FALSE)
}

mat <- comb[, use_cells, drop = FALSE]

if (!identical(colnames(mat), use_cells)) {
  stop("matrix column-order mismatch", call. = FALSE)
}

annotation <- data.frame(
  group = ifelse(
    use_cells %in% refcells,
    "reference",
    "observation"
  ),
  row.names = use_cells,
  stringsAsFactors = FALSE
)

write.table(
  data.frame(
    cell_id = rownames(annotation),
    group = annotation$group
  ),
  file = file.path(outdir, "annotations.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

message(
  sprintf(
    "matrix: %d genes x %d cells; mixture=%d (A=%d, B=%d), reference-B=%d",
    nrow(mat),
    ncol(mat),
    length(mix),
    length(mix_in_A),
    length(mix_in_B),
    length(refcells)
  )
)

obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = mat,
  annotations_file = annotation,
  gene_order_file = gene_order,
  ref_group_names = "reference",
  min_max_counts_per_cell = c(1, Inf),
  chr_exclude = c("chrX", "chrY", "chrM")
)

created_cells <- colnames(obj@expr.data)

missing_after_create <- setdiff(use_cells, created_cells)
extra_after_create <- setdiff(created_cells, use_cells)

if (length(missing_after_create) || length(extra_after_create)) {
  stop(
    ds,
    ": unexpected cell filtering/relabeling during CreateInfercnvObject; missing=",
    length(missing_after_create),
    ", extra=",
    length(extra_after_create),
    call. = FALSE
  )
}

message(
  "CreateInfercnvObject retained all ",
  length(created_cells),
  " requested cells"
)

obj <- infercnv::run(
  obj,
  cutoff = 0.1,
  out_dir = outdir,
  cluster_by_groups = FALSE,
  cluster_references = FALSE,
  denoise = TRUE,
  HMM = FALSE,
  analysis_mode = "samples",
  num_threads = cores,
  plot_steps = FALSE,
  no_prelim_plot = TRUE,
  no_plot = TRUE,
  resume_mode = FALSE
)

expr <- obj@expr.data
go_final <- obj@gene_order

if (
  is.null(expr) ||
  !(is.matrix(expr) || inherits(expr, "Matrix")) ||
  nrow(expr) == 0L ||
  ncol(expr) == 0L
) {
  stop("inferCNV returned an empty final expr.data matrix", call. = FALSE)
}

if (!identical(rownames(expr), rownames(go_final))) {
  stop("inferCNV expr.data and gene_order rows are not aligned", call. = FALSE)
}

if (any(!is.finite(expr))) {
  stop("inferCNV expr.data contains non-finite values", call. = FALSE)
}

missing_final_cells <- setdiff(use_cells, colnames(expr))
extra_final_cells <- setdiff(colnames(expr), use_cells)

if (length(missing_final_cells) || length(extra_final_cells)) {
  stop(
    ds,
    ": final inferCNV matrix cell mismatch; missing=",
    length(missing_final_cells),
    ", extra=",
    length(extra_final_cells),
    call. = FALSE
  )
}

expr <- expr[, use_cells, drop = FALSE]

neutral_center <- 1
burden_unmasked <- colMeans(abs(expr - neutral_center))

mask <- read_mask(mask_file)
masked_gene <- rep(FALSE, nrow(expr))
mask_applied <- FALSE

if (!is.null(mask)) {
  go_chr <- as.character(go_final$chr)
  go_start <- suppressWarnings(as.numeric(go_final$start))
  go_stop <- suppressWarnings(as.numeric(go_final$stop))

  if (
    anyNA(go_chr) ||
    anyNA(go_start) ||
    anyNA(go_stop)
  ) {
    stop("invalid final inferCNV gene-order coordinates", call. = FALSE)
  }

  for (i in seq_len(nrow(mask))) {
    masked_gene <- masked_gene | (
      go_chr == mask$CHROM[i] &
      go_stop >= mask$START[i] &
      go_start <= mask$END[i]
    )
  }

  if (all(masked_gene)) {
    stop("mask removes every final inferCNV gene", call. = FALSE)
  }

  burden_masked <- colMeans(
    abs(expr[!masked_gene, , drop = FALSE] - neutral_center)
  )

  mask_applied <- TRUE

  message(
    "mask applied: ",
    sum(masked_gene),
    "/",
    nrow(expr),
    " inferCNV genes excluded from burden"
  )
} else {
  burden_masked <- rep(
    NA_real_,
    ncol(expr)
  )
  names(burden_masked) <- colnames(expr)

  message(
    "no --mask supplied: infercnv_burden_masked will be NA"
  )
}

bd <- data.table(
  cell_id = colnames(expr),
  infercnv_burden_unmasked = as.numeric(
    burden_unmasked[colnames(expr)]
  ),
  infercnv_burden_masked = as.numeric(
    burden_masked[colnames(expr)]
  )
)

bd[, dataset_id := ds]
bd[, in_mixture := cell_id %in% mix]
bd[, in_reference := cell_id %in% refcells]

mix_meta <- copy(sub)

keep_meta <- setdiff(
  names(mix_meta),
  c("dataset_id")
)

setnames(
  mix_meta,
  old = "cell_id",
  new = "cell_id"
)

bd <- merge(
  bd,
  mix_meta[
    ,
    c(
      list(cell_id = cell_id),
      .SD
    ),
    .SDcols = setdiff(
      keep_meta,
      "cell_id"
    )
  ],
  by = "cell_id",
  all.x = TRUE,
  sort = FALSE
)

if (nrow(bd) != length(use_cells)) {
  stop("output row count mismatch after manifest merge", call. = FALSE)
}

if (anyNA(bd$dataset_id)) {
  bd[
    in_reference == TRUE,
    dataset_id := ds
  ]
}

primary_score_col <- if (mask_applied) {
  "infercnv_burden_masked"
} else {
  "infercnv_burden_unmasked"
}

m <- bd[in_mixture == TRUE]

if (nrow(m) != length(mix)) {
  stop("mixture output row count does not match manifest", call. = FALSE)
}

if (anyNA(m$truth) || anyNA(m$y)) {
  stop("mixture truth labels missing after output merge", call. = FALSE)
}

primary_score <- m[[primary_score_col]]

if (any(!is.finite(primary_score))) {
  stop("primary inferCNV burden contains non-finite mixture values", call. = FALSE)
}

auc_value <- auroc(
  primary_score,
  m$y
)

summary <- data.table(
  dataset_id = ds,
  tumor_fraction = unique(m$tumor_fraction),
  replicate = unique(m$replicate),
  n_mixture = nrow(m),
  n_tumor = sum(m$y == 1L),
  n_normal = sum(m$y == 0L),
  n_reference = length(refcells),
  n_final_genes = nrow(expr),
  n_masked_genes = sum(masked_gene),
  mask_applied = mask_applied,
  burden_column = primary_score_col,
  AUROC = auc_value,
  median_burden_tumor = if (sum(m$y == 1L)) {
    median(primary_score[m$y == 1L])
  } else {
    NA_real_
  },
  median_burden_normal = if (sum(m$y == 0L)) {
    median(primary_score[m$y == 0L])
  } else {
    NA_real_
  },
  infercnv_version = as.character(
    packageVersion("infercnv")
  )
)

burden_path <- file.path(
  outdir,
  paste0(ds, "_infercnv_burden.tsv")
)

expr_path <- file.path(
  outdir,
  paste0(ds, "_infercnv_exprdata.rds")
)

gene_order_path <- file.path(
  outdir,
  paste0(ds, "_infercnv_gene_order_final.tsv")
)

object_path <- file.path(
  outdir,
  paste0(ds, "_infercnv_object_final.rds")
)

summary_path <- file.path(
  outdir,
  paste0(ds, "_infercnv_summary.tsv")
)

session_path <- file.path(
  outdir,
  "sessionInfo.txt"
)

fwrite(
  bd,
  burden_path,
  sep = "\t"
)

saveRDS(
  expr,
  expr_path
)

fwrite(
  data.table(
    gene = rownames(go_final),
    chr = as.character(go_final$chr),
    start = as.numeric(go_final$start),
    stop = as.numeric(go_final$stop),
    masked_from_primary_burden = masked_gene
  ),
  gene_order_path,
  sep = "\t"
)

saveRDS(
  obj,
  object_path
)

fwrite(
  summary,
  summary_path,
  sep = "\t"
)

capture.output(
  sessionInfo(),
  file = session_path
)

writeLines(
  c(
    paste0("dataset_id\t", ds),
    paste0(
      "infercnv_version\t",
      as.character(packageVersion("infercnv"))
    ),
    "cutoff\t0.1",
    "denoise\tTRUE",
    "HMM\tFALSE",
    "analysis_mode\tsamples",
    "cluster_by_groups\tFALSE",
    paste0("mixture_cells\t", length(mix)),
    paste0("reference_B_cells\t", length(refcells)),
    paste0("final_genes\t", nrow(expr)),
    paste0("mask_applied\t", mask_applied),
    paste0("n_masked_genes\t", sum(masked_gene)),
    paste0("primary_burden\t", primary_score_col)
  ),
  file.path(
    outdir,
    "RUN_METADATA.tsv"
  )
)

writeLines(
  "DONE",
  file.path(
    outdir,
    "DONE"
  )
)

message(
  sprintf(
    "inferCNV burden AUROC (%s, mixture): %s",
    ds,
    ifelse(
      is.na(auc_value),
      "NA",
      sprintf("%.4f", auc_value)
    )
  )
)

print(
  m[
    ,
    .(
      n = .N,
      median_burden = median(
        get(primary_score_col)
      )
    ),
    by = truth
  ]
)

message("wrote: ", burden_path)
message("wrote: ", expr_path)
message("wrote: ", gene_order_path)
message("wrote: ", object_path)
message("wrote: ", summary_path)
