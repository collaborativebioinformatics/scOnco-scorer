#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SCEVAN)
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
      "usage: Rscript run_scevan.R <dataset_id> ",
      "[--cores 8] [--mask /w/mask_regions.tsv] ",
      "[--inputs /w/inputs] [--outroot /w/scevan_out]"
    ),
    call. = FALSE
  )
}

ds <- args[1L]
if (!grepl("^[A-Za-z0-9_.-]+$", ds)) {
  stop("dataset_id contains unsafe characters: ", ds, call. = FALSE)
}

cores <- as.integer(get_arg("--cores", "8"))
inputs <- get_arg("--inputs", "/w/inputs")
outroot <- get_arg("--outroot", "/w/scevan_out")
mask_file <- get_arg("--mask", NA_character_)
min_ref_after <- as.integer(get_arg("--min_ref_after", "20"))

if (is.na(cores) || cores < 1L) stop("--cores must be >=1", call. = FALSE)
if (is.na(min_ref_after) || min_ref_after < 1L) {
  stop("--min_ref_after must be >=1", call. = FALSE)
}

A_file <- file.path(inputs, "LLU_A_prefixed.rds")
B_file <- file.path(inputs, "LLU_B_prefixed.rds")
manifest_file <- file.path(inputs, "LLU_all_mixtures_with_p00.tsv")
reference_file <- file.path(inputs, "LLU_B_reference_cells.txt")

for (p in c(A_file, B_file, manifest_file, reference_file)) {
  if (!file.exists(p)) stop("missing required input: ", p, call. = FALSE)
  if (file.info(p)$size <= 0) stop("empty required input: ", p, call. = FALSE)
}

outdir <- file.path(outroot, ds)
if (dir.exists(outdir) && file.exists(file.path(outdir, "DONE"))) {
  stop("completed SCEVAN output already exists for ", ds, ": ", outdir, call. = FALSE)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("=== SCEVAN ", ds, " ===")
message("SCEVAN version: ", as.character(packageVersion("SCEVAN")))
message("cores: ", cores)
message("output: ", outdir)

as_count_matrix <- function(x, label) {
  if (!(is.matrix(x) || inherits(x, "Matrix"))) {
    stop(label, " must be a matrix/Matrix object", call. = FALSE)
  }

  x <- as(Matrix::Matrix(x, sparse = TRUE), "dgCMatrix")

  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop(label, " must have gene rownames and cell colnames", call. = FALSE)
  }
  if (anyNA(rownames(x)) || any(!nzchar(rownames(x)))) {
    stop(label, " has NA/empty gene identifiers", call. = FALSE)
  }
  if (anyNA(colnames(x)) || any(!nzchar(colnames(x)))) {
    stop(label, " has NA/empty cell identifiers", call. = FALSE)
  }
  if (anyDuplicated(colnames(x))) {
    stop(label, " has duplicated cell identifiers", call. = FALSE)
  }
  if (any(!is.finite(x@x)) || any(x@x < 0)) {
    stop(label, " contains non-finite or negative counts", call. = FALSE)
  }
  if (any(abs(x@x - round(x@x)) > 1e-8)) {
    stop(label, " contains non-integer values; SCEVAN requires raw counts", call. = FALSE)
  }

  old_cols <- colnames(x)
  gene_id <- sub("\\.\\d+$", "", rownames(x))

  if (anyDuplicated(gene_id)) {
    u <- unique(gene_id)
    map <- match(gene_id, u)

    S <- Matrix::sparseMatrix(
      i = map,
      j = seq_along(gene_id),
      x = 1,
      dims = c(length(u), length(gene_id)),
      dimnames = list(u, NULL)
    )

    x <- S %*% x
    colnames(x) <- old_cols
  } else {
    rownames(x) <- gene_id
  }

  as(x, "dgCMatrix")
}

A <- as_count_matrix(readRDS(A_file), "LLU_A_prefixed.rds")
B <- as_count_matrix(readRDS(B_file), "LLU_B_prefixed.rds")

if (length(intersect(colnames(A), colnames(B))) > 0L) {
  stop("A and B matrices contain overlapping cell identifiers", call. = FALSE)
}

common <- intersect(rownames(A), rownames(B))
if (!length(common)) stop("A and B matrices have no common genes", call. = FALSE)

frac_A <- length(common) / nrow(A)
frac_B <- length(common) / nrow(B)

if (frac_A < 0.90 || frac_B < 0.90) {
  stop(
    sprintf(
      "A/B common-gene overlap is too low: A %.1f%%, B %.1f%%",
      100 * frac_A,
      100 * frac_B
    ),
    call. = FALSE
  )
}

A <- A[common, , drop = FALSE]
B <- B[common, , drop = FALSE]

man <- fread(manifest_file, showProgress = FALSE)

required_manifest <- c("dataset_id", "cell_id", "truth")
missing_manifest <- setdiff(required_manifest, names(man))

if (length(missing_manifest)) {
  stop(
    "manifest missing columns: ",
    paste(missing_manifest, collapse = ", "),
    call. = FALSE
  )
}

sub <- man[dataset_id == ds]

if (!nrow(sub)) stop("dataset_id not found in manifest: ", ds, call. = FALSE)

sub[, cell_id := as.character(cell_id)]
sub[, truth := as.character(truth)]

if (anyNA(sub$cell_id) || any(!nzchar(sub$cell_id))) {
  stop(ds, ": manifest has NA/empty cell_id", call. = FALSE)
}

if (anyDuplicated(sub$cell_id)) {
  dup <- unique(sub$cell_id[duplicated(sub$cell_id)])
  stop(
    ds,
    ": duplicated mixture cell IDs: ",
    paste(head(dup, 10L), collapse = ", "),
    call. = FALSE
  )
}

mix_cells <- sub$cell_id
truth <- setNames(sub$truth, sub$cell_id)

refcells <- trimws(readLines(reference_file, warn = FALSE))
refcells <- refcells[nzchar(refcells)]

if (!length(refcells)) stop("reference-B list is empty", call. = FALSE)
if (anyDuplicated(refcells)) stop("reference-B list contains duplicates", call. = FALSE)

if (length(intersect(mix_cells, refcells)) > 0L) {
  stop(
    ds,
    ": mixture cells overlap the held-out B reference; design leakage detected",
    call. = FALSE
  )
}

all_matrix_cells <- c(colnames(A), colnames(B))

missing_mix <- setdiff(mix_cells, all_matrix_cells)
if (length(missing_mix)) {
  stop(
    ds,
    ": ",
    length(missing_mix),
    " mixture cells are absent from A/B matrices; first: ",
    paste(head(missing_mix, 10L), collapse = ", "),
    call. = FALSE
  )
}

missing_ref <- setdiff(refcells, colnames(B))
if (length(missing_ref)) {
  stop(
    ds,
    ": ",
    length(missing_ref),
    " reference cells are absent from the B matrix; first: ",
    paste(head(missing_ref, 10L), collapse = ", "),
    call. = FALSE
  )
}

if (length(intersect(refcells, colnames(A)))) {
  stop("reference list contains A cells", call. = FALSE)
}

mix_in_A <- mix_cells[mix_cells %in% colnames(A)]
mix_in_B <- mix_cells[mix_cells %in% colnames(B)]

if (length(mix_in_A) + length(mix_in_B) != length(mix_cells)) {
  stop(ds, ": mixture cell source assignment is inconsistent", call. = FALSE)
}

comb <- cbind(A, B)
use_cells <- c(mix_cells, refcells)

if (anyDuplicated(use_cells)) {
  stop(ds, ": duplicate cells after mixture + reference construction", call. = FALSE)
}

mat_sparse <- comb[, use_cells, drop = FALSE]

if (!identical(colnames(mat_sparse), use_cells)) {
  stop(ds, ": matrix column order mismatch", call. = FALSE)
}

dense_gb <- as.numeric(nrow(mat_sparse)) * as.numeric(ncol(mat_sparse)) * 8 / 1024^3
message(
  sprintf(
    "input matrix: %d genes x %d cells; mixture=%d (A=%d, B=%d), reference-B=%d; dense size ~%.2f GiB",
    nrow(mat_sparse),
    ncol(mat_sparse),
    length(mix_cells),
    length(mix_in_A),
    length(mix_in_B),
    length(refcells),
    dense_gb
  )
)

mat <- as.matrix(mat_sparse)
rm(mat_sparse, comb, A, B)
gc()

if (any(!is.finite(mat)) || any(mat < 0)) {
  stop(ds, ": dense count matrix contains invalid values", call. = FALSE)
}

preprocess_fun <- getFromNamespace("preprocessingMtx", "SCEVAN")

if (!all(
  c("count_mtx", "sample", "par_cores", "findConfident", "organism") %in%
    names(formals(preprocess_fun))
)) {
  stop("installed SCEVAN preprocessingMtx interface is incompatible with this wrapper", call. = FALSE)
}

if (!all(
  c(
    "count_mtx", "annot_mtx", "sample", "par_cores", "norm_cell_names",
    "SEGMENTATION_CLASS", "SMOOTH", "beta_vega", "FIXED_NORMAL_CELLS"
  ) %in% names(formals(SCEVAN::classifyTumorCells))
)) {
  stop("installed SCEVAN classifyTumorCells interface is incompatible with this wrapper", call. = FALSE)
}

message("running SCEVAN preprocessing")

res_proc <- preprocess_fun(
  count_mtx = mat,
  sample = ds,
  par_cores = cores,
  findConfident = FALSE,
  organism = "human"
)

if (is.null(res_proc$count_mtx_norm) || is.null(res_proc$count_mtx_annot)) {
  stop(ds, ": preprocessingMtx returned incomplete output", call. = FALSE)
}

processed_cells <- colnames(res_proc$count_mtx_norm)
norm_cells_after <- intersect(refcells, processed_cells)

message(
  "reference-B cells surviving SCEVAN preprocessing: ",
  length(norm_cells_after),
  "/",
  length(refcells)
)

if (length(norm_cells_after) < min_ref_after) {
  stop(
    ds,
    ": only ",
    length(norm_cells_after),
    " reference-B cells survived preprocessing; minimum required=",
    min_ref_after,
    call. = FALSE
  )
}

if (length(norm_cells_after) < 0.8 * length(refcells)) {
  warning(
    ds,
    ": fewer than 80% of reference-B cells survived SCEVAN preprocessing",
    call. = FALSE
  )
}

message("running SCEVAN classification with FIXED_NORMAL_CELLS=FALSE")

res_class <- SCEVAN::classifyTumorCells(
  count_mtx = res_proc$count_mtx_norm,
  annot_mtx = res_proc$count_mtx_annot,
  sample = ds,
  par_cores = cores,
  norm_cell_names = norm_cells_after,
  SEGMENTATION_CLASS = TRUE,
  SMOOTH = TRUE,
  beta_vega = 0.5,
  FIXED_NORMAL_CELLS = FALSE
)

if (!is.list(res_class) || is.null(res_class$tum_cells) || is.null(res_class$CNAmat)) {
  stop(ds, ": classifyTumorCells returned incomplete output", call. = FALSE)
}

class_cells <- colnames(res_class$CNAmat)
if (length(class_cells) < 4L) {
  stop(ds, ": CNAmat has no cell columns", call. = FALSE)
}
class_cells <- class_cells[-c(1L, 2L, 3L)]

unknown_class_cells <- setdiff(class_cells, use_cells)
if (length(unknown_class_cells)) {
  stop(
    ds,
    ": SCEVAN returned unknown cells: ",
    paste(head(unknown_class_cells, 10L), collapse = ", "),
    call. = FALSE
  )
}

class_df <- data.table(
  cell_id = use_cells,
  scevan_class = "filtered"
)

class_df[cell_id %in% class_cells, scevan_class := "normal"]
class_df[cell_id %in% res_class$tum_cells, scevan_class := "tumor"]

cna_file <- file.path(outdir, "output", paste0(ds, "_CNAmtx.RData"))

if (!file.exists(cna_file) || file.info(cna_file)$size <= 0) {
  stop(
    ds,
    ": classifyTumorCells did not produce the expected CNA matrix file: ",
    cna_file,
    call. = FALSE
  )
}

env <- new.env(parent = emptyenv())
loaded_names <- load(cna_file, envir = env)

if ("CNA_mtx_relat" %in% loaded_names) {
  cna_name <- "CNA_mtx_relat"
} else {
  numeric_matrix_names <- loaded_names[
    vapply(
      loaded_names,
      function(nm) {
        obj <- get(nm, envir = env)
        is.matrix(obj) && is.numeric(obj)
      },
      logical(1)
    )
  ]

  if (length(numeric_matrix_names) != 1L) {
    stop(
      ds,
      ": could not uniquely identify the saved SCEVAN CNA matrix; objects: ",
      paste(loaded_names, collapse = ", "),
      call. = FALSE
    )
  }

  cna_name <- numeric_matrix_names[1L]
}

cna <- get(cna_name, envir = env)

if (!is.matrix(cna) || !is.numeric(cna)) {
  stop(ds, ": saved CNA object is not a numeric matrix", call. = FALSE)
}

if (is.null(colnames(cna))) {
  stop(ds, ": saved CNA matrix has no cell column names", call. = FALSE)
}

if (any(!is.finite(cna))) {
  stop(ds, ": saved CNA matrix contains non-finite values", call. = FALSE)
}

if (!setequal(colnames(cna), class_cells)) {
  stop(ds, ": saved CNA matrix cells do not match classification output", call. = FALSE)
}

cna <- cna[, class_cells, drop = FALSE]

score_unmasked <- colMeans(abs(cna))

score_dt <- data.table(
  cell_id = names(score_unmasked),
  scevan_burden_unmasked = as.numeric(score_unmasked),
  scevan_burden_masked = NA_real_
)

mask_applied <- FALSE
n_masked_genes <- 0L

if (!is.na(mask_file) && nzchar(mask_file)) {
  if (!file.exists(mask_file)) stop("--mask file does not exist: ", mask_file, call. = FALSE)

  mask <- fread(mask_file, showProgress = FALSE)
  required_mask <- c("CHROM", "START", "END")
  missing_mask <- setdiff(required_mask, names(mask))

  if (length(missing_mask)) {
    stop("--mask missing columns: ", paste(missing_mask, collapse = ", "), call. = FALSE)
  }

  mask[, CHROM := sub("^chr", "", as.character(CHROM), ignore.case = TRUE)]
  mask[, START := suppressWarnings(as.numeric(START))]
  mask[, END := suppressWarnings(as.numeric(END))]

  if (
    anyNA(mask$CHROM) ||
    anyNA(mask$START) ||
    anyNA(mask$END) ||
    any(mask$START < 1) ||
    any(mask$END < mask$START)
  ) {
    stop("invalid mask coordinates", call. = FALSE)
  }

  annot <- as.data.table(res_proc$count_mtx_annot)

  required_annot <- c("gene_name", "seqnames", "start", "end")
  missing_annot <- setdiff(required_annot, names(annot))

  if (length(missing_annot)) {
    stop(
      "SCEVAN annotation missing columns required for masking: ",
      paste(missing_annot, collapse = ", "),
      call. = FALSE
    )
  }

  gene_names <- rownames(cna)

  if (is.null(gene_names)) {
    if (nrow(cna) != nrow(annot)) {
      stop(ds, ": CNA matrix lacks rownames and cannot be aligned to annotation", call. = FALSE)
    }
    annot_use <- annot
  } else {
    idx <- match(gene_names, annot$gene_name)

    if (anyNA(idx)) {
      stop(
        ds,
        ": ",
        sum(is.na(idx)),
        " CNA rows cannot be matched to SCEVAN gene annotation",
        call. = FALSE
      )
    }

    annot_use <- annot[idx]
  }

  chr <- sub("^chr", "", as.character(annot_use$seqnames), ignore.case = TRUE)
  start <- suppressWarnings(as.numeric(annot_use$start))
  end <- suppressWarnings(as.numeric(annot_use$end))

  if (anyNA(chr) || anyNA(start) || anyNA(end)) {
    stop(ds, ": invalid SCEVAN gene coordinates", call. = FALSE)
  }

  masked_gene <- rep(FALSE, nrow(cna))

  for (i in seq_len(nrow(mask))) {
    masked_gene <- masked_gene | (
      chr == mask$CHROM[i] &
      end >= mask$START[i] &
      start <= mask$END[i]
    )
  }

  n_masked_genes <- sum(masked_gene)

  if (all(masked_gene)) {
    stop(ds, ": mask removes every SCEVAN CNA gene", call. = FALSE)
  }

  score_masked <- colMeans(abs(cna[!masked_gene, , drop = FALSE]))

  score_dt[
    match(names(score_masked), cell_id),
    scevan_burden_masked := as.numeric(score_masked)
  ]

  mask_applied <- TRUE

  message(
    "mask applied to continuous SCEVAN burden: ",
    n_masked_genes,
    "/",
    nrow(cna),
    " genes excluded"
  )
} else {
  message("no --mask supplied: masked continuous burden is NA")
}

out_all <- merge(
  class_df,
  score_dt,
  by = "cell_id",
  all.x = TRUE,
  sort = FALSE
)

out_all[, dataset_id := ds]
out_all[, in_mixture := cell_id %in% mix_cells]
out_all[, in_reference := cell_id %in% refcells]
out_all[, truth := truth[cell_id]]

if ("tumor_fraction" %in% names(sub)) {
  tf <- unique(sub$tumor_fraction)
  if (length(tf) == 1L) out_all[, tumor_fraction := tf]
}

if ("replicate" %in% names(sub)) {
  rr <- unique(sub$replicate)
  if (length(rr) == 1L) out_all[, replicate := rr]
}

for (cc in c("source_barcode", "source_replicate")) {
  if (cc %in% names(sub)) {
    map <- setNames(sub[[cc]], sub$cell_id)
    out_all[, (cc) := map[cell_id]]
  }
}

out_mix <- out_all[in_mixture == TRUE]

if (nrow(out_mix) != length(mix_cells)) {
  stop(
    ds,
    ": mixture prediction output has ",
    nrow(out_mix),
    " rows but manifest has ",
    length(mix_cells),
    call. = FALSE
  )
}

if (anyNA(out_mix$truth)) {
  stop(ds, ": truth labels are missing for mixture cells", call. = FALSE)
}

pred_path <- file.path(outdir, paste0(ds, "_scevan_predictions.tsv"))
mix_path <- file.path(outdir, paste0(ds, "_scevan_mixture_predictions.tsv"))
full_path <- file.path(outdir, paste0(ds, "_scevan_full.rds"))
cna_rds_path <- file.path(outdir, paste0(ds, "_scevan_CNAmtx.rds"))
summary_path <- file.path(outdir, paste0(ds, "_scevan_summary.tsv"))

fwrite(out_all, pred_path, sep = "\t")
fwrite(out_mix, mix_path, sep = "\t")

saveRDS(
  list(
    preprocessing = res_proc,
    classification = res_class,
    class_table = class_df,
    mask_applied = mask_applied,
    mask_file = if (mask_applied) normalizePath(mask_file) else NA_character_
  ),
  full_path
)

saveRDS(cna, cna_rds_path)

summary <- data.table(
  dataset_id = ds,
  n_manifest_mixture = length(mix_cells),
  n_input_reference_B = length(refcells),
  n_processed_reference_B = length(norm_cells_after),
  n_mixture_filtered = sum(out_mix$scevan_class == "filtered"),
  n_mixture_normal = sum(out_mix$scevan_class == "normal"),
  n_mixture_tumor = sum(out_mix$scevan_class == "tumor"),
  n_reference_filtered = sum(out_all$in_reference & out_all$scevan_class == "filtered"),
  n_reference_normal = sum(out_all$in_reference & out_all$scevan_class == "normal"),
  n_reference_tumor = sum(out_all$in_reference & out_all$scevan_class == "tumor"),
  n_CNA_genes = nrow(cna),
  n_masked_CNA_genes = n_masked_genes,
  mask_applied = mask_applied,
  SCEVAN_version = as.character(packageVersion("SCEVAN"))
)

fwrite(summary, summary_path, sep = "\t")
capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo.txt"))

writeLines(
  c(
    paste0("dataset_id\t", ds),
    paste0("SCEVAN_version\t", as.character(packageVersion("SCEVAN"))),
    paste0("mixture_cells\t", length(mix_cells)),
    paste0("reference_B_input\t", length(refcells)),
    paste0("reference_B_after_preprocessing\t", length(norm_cells_after)),
    paste0("mask_applied\t", mask_applied),
    paste0("n_masked_CNA_genes\t", n_masked_genes)
  ),
  file.path(outdir, "RUN_METADATA.tsv")
)

writeLines("DONE", file.path(outdir, "DONE"))

message("wrote: ", pred_path)
message("wrote: ", mix_path)
message("wrote: ", full_path)
message("wrote: ", cna_rds_path)
message("wrote: ", summary_path)

message("mixture classification:")
print(table(out_mix$scevan_class, out_mix$truth, useNA = "ifany"))

message("reference-B classification:")
print(table(out_all[in_reference == TRUE]$scevan_class, useNA = "ifany"))
