#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2L) {
  stop("usage: Rscript arm_all.R <scevan|infercnv|copykat> <dataset_id>", call. = FALSE)
}

caller <- tolower(args[1L])
ds <- args[2L]

if (!caller %in% c("scevan", "infercnv", "copykat")) {
  stop("caller must be one of: scevan, infercnv, copykat", call. = FALSE)
}

if (!grepl("^[A-Za-z0-9_.-]+$", ds)) {
  stop("unsafe dataset_id: ", ds, call. = FALSE)
}

INPUT <- "/w/inputs"
OUTROOT <- "/w/arm_out"
CORES <- 8L
REF_ALPHA <- 0.05
MIN_REF_CELLS <- 20L
MIN_FEATURES_ARM <- 5L

A_FILE <- file.path(INPUT, "LLU_A_prefixed.rds")
B_FILE <- file.path(INPUT, "LLU_B_prefixed.rds")
MANIFEST_FILE <- file.path(INPUT, "manifest.tsv")
REF_FILE <- file.path(INPUT, "refcells.txt")
GTF_FILE <- file.path(INPUT, "genes.gtf")
GENE_ORDER_FILE <- file.path(INPUT, "gene_order_ensembl.txt")
TRUTH_FILE <- file.path(INPUT, "truth_arms.csv")

for (p in c(A_FILE, B_FILE, MANIFEST_FILE, REF_FILE, GTF_FILE, TRUTH_FILE)) {
  if (!file.exists(p) || file.info(p)$size <= 0) {
    stop("missing/empty required input: ", p, call. = FALSE)
  }
}

if (caller == "infercnv" && (!file.exists(GENE_ORDER_FILE) || file.info(GENE_ORDER_FILE)$size <= 0)) {
  stop("missing/empty inferCNV gene-order file: ", GENE_ORDER_FILE, call. = FALSE)
}

wd <- file.path(OUTROOT, paste0(caller, "_", ds))

if (dir.exists(wd)) {
  unlink(wd, recursive = TRUE, force = TRUE)
}

dir.create(wd, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(wd, "output"), recursive = TRUE, showWarnings = FALSE)

message("=== ARM FIDELITY: ", caller, " ", ds, " ===")

collapse_duplicate_rows <- function(x, ids) {
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
    stop(label, " must be a matrix/Matrix", call. = FALSE)
  }

  x <- as(Matrix::Matrix(x, sparse = TRUE), "dgCMatrix")

  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop(label, " requires gene rownames and cell colnames", call. = FALSE)
  }

  if (
    anyNA(rownames(x)) ||
    any(!nzchar(rownames(x))) ||
    anyNA(colnames(x)) ||
    any(!nzchar(colnames(x)))
  ) {
    stop(label, " contains NA/empty identifiers", call. = FALSE)
  }

  if (anyDuplicated(colnames(x))) {
    stop(label, " contains duplicated cell IDs", call. = FALSE)
  }

  if (
    any(!is.finite(x@x)) ||
    any(x@x < 0) ||
    any(abs(x@x - round(x@x)) > 1e-8)
  ) {
    stop(label, " must contain raw nonnegative integer counts", call. = FALSE)
  }

  ids <- sub("\\.\\d+$", "", rownames(x))

  if (anyDuplicated(ids)) {
    message(
      label,
      ": collapsing ",
      length(ids) - length(unique(ids)),
      " duplicate rows after Ensembl-version stripping"
    )
    x <- collapse_duplicate_rows(x, ids)
  } else {
    rownames(x) <- ids
  }

  if (anyDuplicated(rownames(x))) {
    stop(label, ": duplicate genes remain after collapse", call. = FALSE)
  }

  x
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

as_logical_strict <- function(x, label) {
  if (is.logical(x)) {
    if (anyNA(x)) stop(label, " contains NA", call. = FALSE)
    return(x)
  }

  z <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(z))
  out[z %in% c("true", "t", "1", "yes", "y")] <- TRUE
  out[z %in% c("false", "f", "0", "no", "n")] <- FALSE

  if (anyNA(out)) {
    stop(
      label,
      " contains unrecognized values: ",
      paste(unique(z[is.na(out)]), collapse = ", "),
      call. = FALSE
    )
  }

  out
}

q_start_bp <- c(
  "1" = 123400000,
  "2" = 93900000,
  "3" = 90900000,
  "4" = 50000000,
  "5" = 48800000,
  "6" = 59800000,
  "7" = 60100000,
  "8" = 45200000,
  "9" = 43000000,
  "10" = 39800000,
  "11" = 53400000,
  "12" = 35500000,
  "13" = 17700000,
  "14" = 17200000,
  "15" = 19000000,
  "16" = 36800000,
  "17" = 25100000,
  "18" = 18500000,
  "19" = 26200000,
  "20" = 28100000,
  "21" = 12000000,
  "22" = 15000000,
  "X" = 61000000
)

assign_arm <- function(chr, pos) {
  chr <- sub("^chr", "", as.character(chr), ignore.case = TRUE)
  pos <- suppressWarnings(as.numeric(pos))

  out <- rep(NA_character_, length(chr))

  ok <- chr %in% names(q_start_bp) & is.finite(pos)

  out[ok] <- paste0(
    chr[ok],
    ifelse(
      pos[ok] < unname(q_start_bp[chr[ok]]),
      "p",
      "q"
    )
  )

  out
}

read_gtf_gene_map <- function(path) {
  cmd <- sprintf(
    "awk -F '\\t' '$3==\"gene\"' %s",
    shQuote(path)
  )

  gl <- fread(
    cmd = cmd,
    sep = "\t",
    header = FALSE,
    quote = "",
    fill = TRUE,
    showProgress = FALSE
  )

  if (nrow(gl) == 0L || ncol(gl) < 9L) {
    stop("failed to parse gene features from GTF", call. = FALSE)
  }

  attrs <- as.character(gl[[9L]])
  has_gid <- grepl('gene_id "[^"]+"', attrs)

  gl <- gl[has_gid]
  attrs <- attrs[has_gid]

  gid <- sub(
    '.*gene_id "([^"]+)".*',
    '\\1',
    attrs
  )
  gid <- sub("\\.\\d+$", "", gid)

  chr <- sub("^chr", "", as.character(gl[[1L]]), ignore.case = TRUE)
  start <- suppressWarnings(as.numeric(gl[[4L]]))
  stop <- suppressWarnings(as.numeric(gl[[5L]]))

  gm <- data.table(
    gene = gid,
    chr = chr,
    start = start,
    stop = stop
  )

  gm <- gm[
    !is.na(gene) &
      nzchar(gene) &
      chr %in% names(q_start_bp) &
      is.finite(start) &
      is.finite(stop) &
      start >= 1 &
      stop >= start
  ]

  cross_chr <- gm[
    ,
    .(n_chr = uniqueN(chr)),
    by = gene
  ][n_chr > 1L]

  if (nrow(cross_chr) > 0L) {
    stop(
      "version-stripped GTF genes map to multiple chromosomes; first: ",
      paste(head(cross_chr$gene, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  gm <- gm[
    ,
    .(
      chr = first(chr),
      start = min(start),
      stop = max(stop)
    ),
    by = gene
  ]

  gm[, midpoint := (start + stop) / 2]
  gm[, arm := assign_arm(chr, midpoint)]

  gm[!is.na(arm)]
}

make_arm_cell_table <- function(signal_matrix, feature_arm, cells) {
  if (!length(cells)) {
    stop("zero requested cells for arm aggregation", call. = FALSE)
  }

  cells <- intersect(cells, colnames(signal_matrix))

  if (!length(cells)) {
    stop("none of the requested cells are present in signal matrix", call. = FALSE)
  }

  if (length(feature_arm) != nrow(signal_matrix)) {
    stop("feature-arm vector length does not match signal matrix rows", call. = FALSE)
  }

  valid <- !is.na(feature_arm) & nzchar(feature_arm)

  signal_matrix <- signal_matrix[valid, cells, drop = FALSE]
  feature_arm <- feature_arm[valid]

  arm_levels <- sort(unique(feature_arm))

  out <- rbindlist(
    lapply(
      arm_levels,
      function(a) {
        idx <- which(feature_arm == a)

        if (length(idx) < MIN_FEATURES_ARM) {
          return(NULL)
        }

        vals <- colMeans(
          signal_matrix[idx, , drop = FALSE]
        )

        data.table(
          arm = a,
          cell_id = names(vals),
          arm_signal = as.numeric(vals),
          n_features = length(idx)
        )
      }
    )
  )

  if (!nrow(out)) {
    stop("no chromosome arms passed MIN_FEATURES_ARM", call. = FALSE)
  }

  out
}

summarize_arm_signal <- function(tum_table, ref_table) {
  tum <- tum_table[
    ,
    .(
      n_tumor_cells = .N,
      n_features = first(n_features),
      tumor_mean_signal = mean(arm_signal),
      tumor_median_signal = median(arm_signal),
      tumor_sd_signal = if (.N > 1L) sd(arm_signal) else NA_real_
    ),
    by = arm
  ]

  ref <- ref_table[
    ,
    .(
      n_ref_cells = .N,
      ref_mean_signal = mean(arm_signal),
      ref_median_signal = median(arm_signal),
      ref_sd_signal = if (.N > 1L) sd(arm_signal) else NA_real_,
      ref_lower = as.numeric(
        quantile(
          arm_signal,
          probs = REF_ALPHA / 2,
          names = FALSE,
          type = 8
        )
      ),
      ref_upper = as.numeric(
        quantile(
          arm_signal,
          probs = 1 - REF_ALPHA / 2,
          names = FALSE,
          type = 8
        )
      )
    ),
    by = arm
  ]

  z <- merge(
    tum,
    ref,
    by = "arm",
    all = FALSE
  )

  if (!nrow(z)) {
    stop("no arms have both tumor and reference signal", call. = FALSE)
  }

  z[, delta_vs_ref := tumor_median_signal - ref_median_signal]

  z[
    ,
    caller_direction := fifelse(
      delta_vs_ref > 0,
      "gain",
      fifelse(
        delta_vs_ref < 0,
        "loss",
        "zero"
      )
    )
  ]

  z[
    ,
    caller_state_ref95 := fifelse(
      tumor_median_signal > ref_upper,
      "gain",
      fifelse(
        tumor_median_signal < ref_lower,
        "loss",
        "neutral"
      )
    )
  ]

  z
}

A <- prepare_counts(
  readRDS(A_FILE),
  "LLU_A_prefixed.rds"
)

B <- prepare_counts(
  readRDS(B_FILE),
  "LLU_B_prefixed.rds"
)

if (length(intersect(colnames(A), colnames(B))) > 0L) {
  stop("A and B matrices share cell IDs", call. = FALSE)
}

common <- intersect(rownames(A), rownames(B))

if (!length(common)) {
  stop("A and B matrices have zero common genes", call. = FALSE)
}

overlap_A <- length(common) / nrow(A)
overlap_B <- length(common) / nrow(B)

if (overlap_A < 0.90 || overlap_B < 0.90) {
  stop(
    sprintf(
      "A/B common-gene overlap too low: A %.2f%% B %.2f%%",
      100 * overlap_A,
      100 * overlap_B
    ),
    call. = FALSE
  )
}

A <- A[common, , drop = FALSE]
B <- B[common, , drop = FALSE]

man <- fread(MANIFEST_FILE, showProgress = FALSE)

required_man <- c(
  "dataset_id",
  "cell_id",
  "truth",
  "tumor_fraction",
  "replicate"
)

missing_man <- setdiff(required_man, names(man))

if (length(missing_man)) {
  stop(
    "manifest missing columns: ",
    paste(missing_man, collapse = ", "),
    call. = FALSE
  )
}

sub <- man[dataset_id == ds]

if (!nrow(sub)) {
  stop("dataset not found in manifest: ", ds, call. = FALSE)
}

sub[, cell_id := as.character(cell_id)]
sub[, y := normalize_truth(truth)]
sub[, tumor_fraction := suppressWarnings(as.numeric(tumor_fraction))]
sub[, replicate := as.character(replicate)]

if (
  anyNA(sub$cell_id) ||
  any(!nzchar(sub$cell_id)) ||
  anyDuplicated(sub$cell_id)
) {
  stop(ds, ": invalid/duplicated mixture cell IDs", call. = FALSE)
}

if (
  anyNA(sub$tumor_fraction) ||
  length(unique(sub$tumor_fraction)) != 1L
) {
  stop(ds, ": tumor_fraction must be one valid value", call. = FALSE)
}

if (length(unique(sub$replicate)) != 1L) {
  stop(ds, ": replicate must be one value", call. = FALSE)
}

mix <- sub$cell_id
tum <- sub[y == 1L, cell_id]
normal_mix <- sub[y == 0L, cell_id]

ref <- trimws(
  readLines(
    REF_FILE,
    warn = FALSE
  )
)
ref <- ref[nzchar(ref)]

if (!length(ref)) {
  stop("reference-B list is empty", call. = FALSE)
}

if (anyDuplicated(ref)) {
  stop("reference-B list contains duplicates", call. = FALSE)
}

if (length(intersect(mix, ref)) > 0L) {
  stop(ds, ": mixture/reference overlap detected", call. = FALSE)
}

all_cells <- c(
  colnames(A),
  colnames(B)
)

missing_mix <- setdiff(
  mix,
  all_cells
)

if (length(missing_mix)) {
  stop(
    ds,
    ": mixture cells absent from A/B matrices; first: ",
    paste(head(missing_mix, 10L), collapse = ", "),
    call. = FALSE
  )
}

missing_ref <- setdiff(
  ref,
  colnames(B)
)

if (length(missing_ref)) {
  stop(
    ds,
    ": reference-B cells absent from B matrix; first: ",
    paste(head(missing_ref, 10L), collapse = ", "),
    call. = FALSE
  )
}

if (length(intersect(ref, colnames(A))) > 0L) {
  stop("reference list contains A cells", call. = FALSE)
}

if (!all(tum %in% colnames(A))) {
  stop(ds, ": at least one truth=tumor cell is not in A matrix", call. = FALSE)
}

if (!all(normal_mix %in% colnames(B))) {
  stop(ds, ": at least one truth=normal cell is not in B matrix", call. = FALSE)
}

comb <- cbind(
  A,
  B
)

use <- c(
  mix,
  ref
)

if (anyDuplicated(use)) {
  stop("duplicate cells after mixture + reference construction", call. = FALSE)
}

mat_sparse <- comb[
  ,
  use,
  drop = FALSE
]

if (!identical(colnames(mat_sparse), use)) {
  stop("input matrix column order mismatch", call. = FALSE)
}

dense_gb <- (
  as.numeric(nrow(mat_sparse)) *
    as.numeric(ncol(mat_sparse)) *
    8
) / 1024^3

message(
  sprintf(
    "matrix=%d genes x %d cells; mixture=%d tumor=%d mixture-normal=%d reference-B=%d dense~%.2f GiB",
    nrow(mat_sparse),
    ncol(mat_sparse),
    length(mix),
    length(tum),
    length(normal_mix),
    length(ref),
    dense_gb
  )
)

gtf_map <- read_gtf_gene_map(GTF_FILE)

signal_matrix <- NULL
feature_arm <- NULL
signal_center <- NA_real_
signal_source <- NA_character_
caller_version <- NA_character_
native_target_cells <- character()
native_ref_cells <- character()

if (caller == "scevan") {
  suppressPackageStartupMessages(
    library(SCEVAN)
  )

  caller_version <- as.character(
    packageVersion("SCEVAN")
  )

  mat <- as.matrix(mat_sparse)

  preprocess_fun <- getFromNamespace(
    "preprocessingMtx",
    "SCEVAN"
  )

  res_proc <- preprocess_fun(
    count_mtx = mat,
    sample = ds,
    par_cores = CORES,
    findConfident = FALSE,
    organism = "human",
    output_dir = file.path(wd, "output")
  )

  if (
    is.null(res_proc$count_mtx_norm) ||
    is.null(res_proc$count_mtx_annot)
  ) {
    stop("SCEVAN preprocessing returned incomplete output", call. = FALSE)
  }

  norm_after <- intersect(
    ref,
    colnames(
      res_proc$count_mtx_norm
    )
  )

  if (length(norm_after) < MIN_REF_CELLS) {
    stop(
      "only ",
      length(norm_after),
      " reference-B cells survived SCEVAN preprocessing; minimum=",
      MIN_REF_CELLS,
      call. = FALSE
    )
  }

  res_class <- SCEVAN::classifyTumorCells(
    count_mtx = res_proc$count_mtx_norm,
    annot_mtx = res_proc$count_mtx_annot,
    sample = ds,
    par_cores = CORES,
    norm_cell_names = norm_after,
    SEGMENTATION_CLASS = TRUE,
    SMOOTH = TRUE,
    beta_vega = 0.5,
    FIXED_NORMAL_CELLS = FALSE,
    output_dir = file.path(wd, "output")
  )

  cna_file <- file.path(
    wd,
    "output",
    paste0(
      ds,
      "_CNAmtx.RData"
    )
  )

  if (!file.exists(cna_file) || file.info(cna_file)$size <= 0) {
    stop("SCEVAN CNA matrix file missing: ", cna_file, call. = FALSE)
  }

  e <- new.env(
    parent = emptyenv()
  )

  loaded <- load(
    cna_file,
    envir = e
  )

  if (!"CNA_mtx_relat" %in% loaded) {
    stop(
      "SCEVAN CNA file does not contain expected object CNA_mtx_relat; found: ",
      paste(loaded, collapse = ", "),
      call. = FALSE
    )
  }

  cna <- get(
    "CNA_mtx_relat",
    envir = e
  )

  if (
    !(is.matrix(cna) || inherits(cna, "Matrix")) ||
    !is.numeric(cna)
  ) {
    stop("SCEVAN CNA_mtx_relat is not a numeric matrix", call. = FALSE)
  }

  if (any(!is.finite(cna))) {
    stop("SCEVAN CNA_mtx_relat contains non-finite values", call. = FALSE)
  }

  annot <- as.data.table(
    res_proc$count_mtx_annot
  )

  required_annot <- c(
    "seqnames",
    "start",
    "end"
  )

  missing_annot <- setdiff(
    required_annot,
    names(annot)
  )

  if (length(missing_annot)) {
    stop(
      "SCEVAN annotation missing columns: ",
      paste(missing_annot, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(annot) != nrow(cna)) {
    stop(
      "SCEVAN annotation/CNA row mismatch: ",
      nrow(annot),
      " vs ",
      nrow(cna),
      call. = FALSE
    )
  }

  feature_arm <- assign_arm(
    annot$seqnames,
    (
      as.numeric(annot$start) +
        as.numeric(annot$end)
    ) / 2
  )

  signal_matrix <- cna
  signal_center <- 0
  signal_source <- "SCEVAN_CNA_mtx_relat"

  native_target_cells <- intersect(
    mix,
    colnames(signal_matrix)
  )

  native_ref_cells <- intersect(
    ref,
    colnames(signal_matrix)
  )

  rm(
    mat,
    res_proc,
    res_class,
    cna
  )
  gc()
}

if (caller == "infercnv") {
  suppressPackageStartupMessages(
    library(infercnv)
  )

  caller_version <- as.character(
    packageVersion("infercnv")
  )

  ann <- data.frame(
    cell = use,
    group = ifelse(
      use %in% ref,
      "reference",
      "observation"
    ),
    stringsAsFactors = FALSE
  )

  ann_file <- file.path(
    wd,
    "annotations.txt"
  )

  write.table(
    ann,
    ann_file,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = mat_sparse,
    annotations_file = ann_file,
    gene_order_file = GENE_ORDER_FILE,
    ref_group_names = "reference",
    min_max_counts_per_cell = c(1, Inf),
    chr_exclude = c(
      "chrX",
      "chrY",
      "chrM"
    )
  )

  missing_created <- setdiff(
    use,
    colnames(obj@expr.data)
  )

  if (length(missing_created)) {
    stop(
      "inferCNV object creation removed requested cells; first: ",
      paste(head(missing_created, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  obj <- infercnv::run(
    obj,
    cutoff = 0.1,
    out_dir = file.path(
      wd,
      "infercnv"
    ),
    cluster_by_groups = FALSE,
    cluster_references = FALSE,
    denoise = TRUE,
    HMM = FALSE,
    analysis_mode = "samples",
    num_threads = CORES,
    plot_steps = FALSE,
    no_prelim_plot = TRUE,
    no_plot = TRUE,
    resume_mode = FALSE
  )

  ex <- obj@expr.data
  go <- obj@gene_order

  if (
    is.null(ex) ||
    nrow(ex) == 0L ||
    ncol(ex) == 0L
  ) {
    stop("inferCNV returned empty expr.data", call. = FALSE)
  }

  if (!identical(rownames(ex), rownames(go))) {
    stop("inferCNV expr.data/gene_order row mismatch", call. = FALSE)
  }

  if (any(!is.finite(ex))) {
    stop("inferCNV expr.data contains non-finite values", call. = FALSE)
  }

  feature_arm <- assign_arm(
    go$chr,
    (
      as.numeric(go$start) +
        as.numeric(go$stop)
    ) / 2
  )

  signal_matrix <- ex
  signal_center <- 1
  signal_source <- "inferCNV_final_expr_data"

  native_target_cells <- intersect(
    mix,
    colnames(signal_matrix)
  )

  native_ref_cells <- intersect(
    ref,
    colnames(signal_matrix)
  )

  saveRDS(
    obj,
    file.path(
      wd,
      paste0(
        ds,
        "_infercnv_object_final.rds"
      )
    )
  )

  rm(
    obj,
    ex,
    go
  )
  gc()
}

if (caller == "copykat") {
  suppressPackageStartupMessages(
    library(copykat)
  )

  caller_version <- as.character(
    packageVersion("copykat")
  )

  mat <- as.matrix(
    mat_sparse
  )

  ck <- copykat::copykat(
    rawmat = mat,
    id.type = "E",
    cell.line = "no",
    ngene.chr = 5,
    win.size = 25,
    KS.cut = 0.1,
    sam.name = ds,
    distance = "euclidean",
    norm.cell.names = ref,
    n.cores = CORES,
    output.seg = "FALSE",
    plot.genes = FALSE,
    genome = "hg20"
  )

  if (
    !is.list(ck) ||
    is.null(ck$CNAmat)
  ) {
    stop("CopyKAT returned no CNAmat", call. = FALSE)
  }

  cm <- ck$CNAmat

  required_meta <- c(
    "chrom",
    "chrompos",
    "abspos"
  )

  missing_meta <- setdiff(
    required_meta,
    colnames(cm)
  )

  if (length(missing_meta)) {
    stop(
      "CopyKAT CNAmat missing metadata columns: ",
      paste(missing_meta, collapse = ", "),
      call. = FALSE
    )
  }

  cell_cols <- setdiff(
    colnames(cm),
    required_meta
  )

  if (!length(cell_cols)) {
    stop("CopyKAT CNAmat has no cell columns", call. = FALSE)
  }

  values <- as.matrix(
    cm[
      ,
      cell_cols,
      drop = FALSE
    ]
  )

  storage.mode(values) <- "double"

  if (any(!is.finite(values))) {
    stop("CopyKAT CNAmat contains non-finite values", call. = FALSE)
  }

  feature_arm <- assign_arm(
    cm$chrom,
    cm$chrompos
  )

  signal_matrix <- values
  signal_center <- 0
  signal_source <- "CopyKAT_CNAmat"

  native_target_cells <- intersect(
    mix,
    colnames(signal_matrix)
  )

  native_ref_cells <- intersect(
    ref,
    colnames(signal_matrix)
  )

  saveRDS(
    ck,
    file.path(
      wd,
      paste0(
        ds,
        "_copykat_full.rds"
      )
    )
  )

  rm(
    mat,
    ck,
    cm,
    values
  )
  gc()
}

if (is.null(signal_matrix) || is.null(feature_arm)) {
  stop("caller did not produce an arm-scoring signal matrix", call. = FALSE)
}

if (length(native_ref_cells) < MIN_REF_CELLS) {
  stop(
    caller,
    ": only ",
    length(native_ref_cells),
    " reference-B cells are present in caller signal matrix; minimum=",
    MIN_REF_CELLS,
    call. = FALSE
  )
}

missing_tumor_signal <- setdiff(
  tum,
  native_target_cells
)

missing_normal_signal <- setdiff(
  normal_mix,
  native_target_cells
)

if (length(missing_tumor_signal) > 0L) {
  message(
    caller,
    ": ",
    length(missing_tumor_signal),
    "/",
    length(tum),
    " true tumor cells were filtered before arm scoring"
  )
}

if (length(missing_normal_signal) > 0L) {
  message(
    caller,
    ": ",
    length(missing_normal_signal),
    "/",
    length(normal_mix),
    " mixture-normal cells were filtered before arm scoring"
  )
}

tum_signal_cells <- intersect(
  tum,
  colnames(signal_matrix)
)

ref_signal_cells <- intersect(
  ref,
  colnames(signal_matrix)
)

if (!length(tum_signal_cells)) {
  stop("zero true tumor cells remain for arm scoring", call. = FALSE)
}

signal_dev <- signal_matrix - signal_center

tum_arm_cells <- make_arm_cell_table(
  signal_dev,
  feature_arm,
  tum_signal_cells
)

ref_arm_cells <- make_arm_cell_table(
  signal_dev,
  feature_arm,
  ref_signal_cells
)

arm_signal <- summarize_arm_signal(
  tum_arm_cells,
  ref_arm_cells
)

truth <- fread(
  TRUTH_FILE,
  showProgress = FALSE
)

required_truth <- c(
  "arm_id",
  "state_call",
  "masked"
)

missing_truth <- setdiff(
  required_truth,
  names(truth)
)

if (length(missing_truth)) {
  stop(
    "truth_arms.csv missing columns: ",
    paste(missing_truth, collapse = ", "),
    call. = FALSE
  )
}

truth[, arm_id := as.character(arm_id)]
truth[, state_call := tolower(trimws(as.character(state_call)))]
truth[, masked_bool := as_logical_strict(masked, "truth$masked")]

valid_truth_states <- c(
  "gain",
  "loss",
  "neutral",
  "mixed"
)

bad_truth_states <- setdiff(
  unique(
    truth[
      !is.na(state_call),
      state_call
    ]
  ),
  valid_truth_states
)

if (length(bad_truth_states)) {
  stop(
    "unexpected truth state_call values: ",
    paste(bad_truth_states, collapse = ", "),
    call. = FALSE
  )
}

truth_eval <- truth[
  masked_bool == FALSE &
    !is.na(state_call)
]

if (!nrow(truth_eval)) {
  stop("no unmasked truth arms available", call. = FALSE)
}

truth_eval[
  ,
  truth_direction := fifelse(
    state_call == "gain",
    "gain",
    fifelse(
      state_call == "loss",
      "loss",
      NA_character_
    )
  )
]

arm_signal[
  ,
  arm_id := arm
]

m <- merge(
  truth_eval,
  arm_signal,
  by = "arm_id",
  all.x = TRUE,
  sort = FALSE
)

m[
  ,
  arm_evaluable :=
    is.finite(tumor_median_signal) &
    is.finite(ref_median_signal) &
    n_tumor_cells > 0 &
    n_ref_cells >= MIN_REF_CELLS &
    n_features >= MIN_FEATURES_ARM
]

gain_truth <- m[
  state_call == "gain"
]

loss_truth <- m[
  state_call == "loss"
]

directional_truth <- m[
  state_call %in% c(
    "gain",
    "loss"
  )
]

evaluable_directional <- directional_truth[
  arm_evaluable == TRUE
]

gain_evaluable <- gain_truth[
  arm_evaluable == TRUE
]

loss_evaluable <- loss_truth[
  arm_evaluable == TRUE
]

safe_pct <- function(num, den) {
  if (den == 0L) return(NA_real_)
  100 * num / den
}

out <- data.table(
  caller = caller,
  caller_version = caller_version,
  dataset = ds,
  tumor_fraction = unique(sub$tumor_fraction),
  replicate = unique(sub$replicate),

  signal_source = signal_source,
  neutral_center = signal_center,
  reference_interval_alpha = REF_ALPHA,

  n_truth_arms_unmasked = nrow(m),
  n_truth_gain = nrow(gain_truth),
  n_truth_loss = nrow(loss_truth),
  n_truth_neutral = sum(m$state_call == "neutral"),
  n_truth_mixed = sum(m$state_call == "mixed"),

  n_directional_arms = nrow(directional_truth),
  n_directional_arms_evaluable = nrow(evaluable_directional),

  n_true_tumor_cells_manifest = length(tum),
  n_true_tumor_cells_scored = length(tum_signal_cells),
  tumor_cell_retention_pct = safe_pct(
    length(tum_signal_cells),
    length(tum)
  ),

  n_reference_B_input = length(ref),
  n_reference_B_scored = length(ref_signal_cells),
  reference_cell_retention_pct = safe_pct(
    length(ref_signal_cells),
    length(ref)
  ),

  gain_direction_recovery_pct = safe_pct(
    sum(
      gain_evaluable$caller_direction == "gain",
      na.rm = TRUE
    ),
    nrow(gain_evaluable)
  ),

  loss_direction_recovery_pct = safe_pct(
    sum(
      loss_evaluable$caller_direction == "loss",
      na.rm = TRUE
    ),
    nrow(loss_evaluable)
  ),

  directional_concordance_pct = safe_pct(
    sum(
      evaluable_directional$caller_direction ==
        evaluable_directional$truth_direction,
      na.rm = TRUE
    ),
    nrow(evaluable_directional)
  ),

  gain_recovery_ref95_pct = safe_pct(
    sum(
      gain_evaluable$caller_state_ref95 == "gain",
      na.rm = TRUE
    ),
    nrow(gain_evaluable)
  ),

  loss_recovery_ref95_pct = safe_pct(
    sum(
      loss_evaluable$caller_state_ref95 == "loss",
      na.rm = TRUE
    ),
    nrow(loss_evaluable)
  ),

  directional_concordance_ref95_pct = safe_pct(
    sum(
      evaluable_directional$caller_state_ref95 ==
        evaluable_directional$truth_direction,
      na.rm = TRUE
    ),
    nrow(evaluable_directional)
  )
)

m[
  ,
  `:=`(
    caller = caller,
    caller_version = caller_version,
    dataset = ds,
    tumor_fraction = unique(sub$tumor_fraction),
    replicate = unique(sub$replicate),
    signal_source = signal_source,
    neutral_center = signal_center
  )
]

setcolorder(
  m,
  c(
    "caller",
    "caller_version",
    "dataset",
    "tumor_fraction",
    "replicate",
    "arm_id",
    setdiff(
      names(m),
      c(
        "caller",
        "caller_version",
        "dataset",
        "tumor_fraction",
        "replicate",
        "arm_id"
      )
    )
  )
)

arm_cells_path <- file.path(
  wd,
  paste0(
    caller,
    "_",
    ds,
    "_tumor_arm_cell_signal.tsv"
  )
)

ref_cells_path <- file.path(
  wd,
  paste0(
    caller,
    "_",
    ds,
    "_reference_arm_cell_signal.tsv"
  )
)

arms_path <- file.path(
  OUTROOT,
  paste0(
    caller,
    "_",
    ds,
    "_arms.tsv"
  )
)

summary_path <- file.path(
  OUTROOT,
  paste0(
    caller,
    "_",
    ds,
    "_armrec.tsv"
  )
)

fwrite(
  tum_arm_cells,
  arm_cells_path,
  sep = "\t"
)

fwrite(
  ref_arm_cells,
  ref_cells_path,
  sep = "\t"
)

fwrite(
  m,
  arms_path,
  sep = "\t"
)

fwrite(
  out,
  summary_path,
  sep = "\t"
)

capture.output(
  sessionInfo(),
  file = file.path(
    wd,
    "sessionInfo.txt"
  )
)

writeLines(
  c(
    paste0("caller\t", caller),
    paste0("caller_version\t", caller_version),
    paste0("dataset\t", ds),
    paste0("tumor_fraction\t", unique(sub$tumor_fraction)),
    paste0("replicate\t", unique(sub$replicate)),
    paste0("signal_source\t", signal_source),
    paste0("neutral_center\t", signal_center),
    paste0("reference_interval_alpha\t", REF_ALPHA),
    paste0("min_reference_cells\t", MIN_REF_CELLS),
    paste0("min_features_per_arm\t", MIN_FEATURES_ARM),
    "arm_boundary_source\tUCSC_hg38_cytoBand_q-arm_acen_start"
  ),
  file.path(
    wd,
    "RUN_METADATA.tsv"
  )
)

writeLines(
  "DONE",
  file.path(
    wd,
    "DONE"
  )
)

cat(
  caller,
  ds,
  sprintf(
    ": direction gain=%.1f%% loss=%s concordance=%.1f%% | ref95 gain=%.1f%% loss=%s concordance=%.1f%% | tumor cells %d/%d\n",
    out$gain_direction_recovery_pct,
    ifelse(
      is.na(out$loss_direction_recovery_pct),
      "NA",
      sprintf("%.1f%%", out$loss_direction_recovery_pct)
    ),
    out$directional_concordance_pct,
    out$gain_recovery_ref95_pct,
    ifelse(
      is.na(out$loss_recovery_ref95_pct),
      "NA",
      sprintf("%.1f%%", out$loss_recovery_ref95_pct)
    ),
    out$directional_concordance_ref95_pct,
    out$n_true_tumor_cells_scored,
    out$n_true_tumor_cells_manifest
  )
)
