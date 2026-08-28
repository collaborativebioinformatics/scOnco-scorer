#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(numbat)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(jsonlite)
})

option_list <- list(
  make_option("--counts", type = "character", default = NA_character_,
              help = "Raw integer gene x run count matrix (.rds) [required]"),
  make_option("--df_allele", type = "character", default = NA_character_,
              help = "Numbat SMART-Seq allele-count TSV/TSV.GZ [required]"),
  make_option("--gtf", type = "character", default = NA_character_,
              help = "Exact GRCh38/Ensembl-v84 GTF used for alignment/counting [required]"),
  make_option("--outdir", type = "character", default = "out/numbat"),
  make_option("--ncores", type = "integer", default = 16L),
  make_option("--ref_frac", type = "double", default = 0.5),
  make_option("--seed", type = "integer", default = 100L),

  make_option("--init_k", type = "integer", default = 2L),
  make_option("--min_cells", type = "integer", default = 5L),
  make_option("--max_iter", type = "integer", default = 2L),

  make_option("--min_gtf_overlap", type = "double", default = 0.90),

  make_option("--segs_loh", type = "character", default = NA_character_),
  make_option("--segs_fix", type = "character", default = NA_character_),
  make_option("--segments_build", type = "character", default = NA_character_,
              help = "Required with --segs_loh/--segs_fix; must be GRCh38"),
  make_option("--segments_coordinates", type = "character", default = NA_character_,
              help = "Required with --segs_loh/--segs_fix; must be 1-based-closed"),

  make_option("--allow_missing_allele_units", action = "store_true", default = FALSE,
              help = "Allow target units with zero allele rows (default: fail)"),
  make_option("--overwrite", action = "store_true", default = FALSE,
              help = "Delete and recreate an existing Numbat output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

log_msg <- function(fmt, ...) {
  cat(sprintf(
    "%s [c1_numbat] %s\n",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    sprintf(fmt, ...)
  ))
}

abort <- function(fmt, ...) {
  stop(sprintf(fmt, ...), call. = FALSE)
}

is_missing <- function(x) {
  is.null(x) || length(x) != 1L || is.na(x) || !nzchar(x)
}

assert_file <- function(path, label) {
  if (is_missing(path)) abort("%s is required", label)
  if (!file.exists(path)) abort("%s does not exist: %s", label, path)
  info <- file.info(path)
  if (is.na(info$size) || info$size <= 0) abort("%s is empty: %s", label, path)
  invisible(path)
}

md5_or_na <- function(path) {
  if (is_missing(path) || !file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}

write_json_pretty <- function(x, path) {
  writeLines(
    jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA),
    con = path
  )
}

write_jsonl_one <- function(x, path) {
  line <- jsonlite::toJSON(x, auto_unbox = TRUE, pretty = FALSE, null = "null", digits = NA)
  writeLines(line, con = path, sep = "\n", useBytes = TRUE)
}

normalize_chrom <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x
}

assert_whole_nonnegative <- function(x, label, tol = 1e-8) {
  if (length(x) == 0L) return(invisible(TRUE))
  if (any(!is.finite(x))) abort("%s contains non-finite values", label)
  if (any(x < 0)) abort("%s contains negative values", label)
  if (any(abs(x - round(x)) > tol)) abort("%s contains non-integer values", label)
  invisible(TRUE)
}

collapse_sparse_rows <- function(x, ids) {
  if (length(ids) != nrow(x)) abort("internal error: row-ID vector length mismatch")
  if (anyNA(ids) || any(!nzchar(ids))) abort("gene identifiers contain NA/empty values")

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
  methods::as(y, "dgCMatrix")
}

read_exact_gene_gtf <- function(path) {
  assert_file(path, "--gtf")

  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    cmd <- sprintf(
      "gzip -dc -- %s | awk -F '\\t' '$3==\"gene\"'",
      shQuote(path)
    )
  } else {
    cmd <- sprintf(
      "awk -F '\\t' '$3==\"gene\"' %s",
      shQuote(path)
    )
  }

  gl <- data.table::fread(
    cmd = cmd,
    sep = "\t",
    header = FALSE,
    quote = "",
    fill = TRUE,
    showProgress = FALSE
  )

  if (nrow(gl) == 0L || ncol(gl) < 9L) {
    abort("could not parse gene records from GTF: %s", path)
  }

  gene_id <- stringr::str_match(gl[[9L]], 'gene_id "([^"]+)"')[, 2L]
  gene_id <- sub("\\.\\d+$", "", gene_id)

  g <- data.table::data.table(
    CHROM = normalize_chrom(gl[[1L]]),
    gene = gene_id,
    gene_start = suppressWarnings(as.integer(gl[[4L]])),
    gene_end = suppressWarnings(as.integer(gl[[5L]]))
  )

  g <- g[
    !is.na(gene) & nzchar(gene) &
      !is.na(gene_start) & !is.na(gene_end) &
      gene_start >= 1L & gene_end >= gene_start
  ]

  g <- g[CHROM %in% as.character(1:22)]

  if (nrow(g) == 0L) abort("no autosomal gene records remain after GTF parsing")

  cross_chr <- g[, .(n_chr = uniqueN(CHROM)), by = gene][n_chr > 1L]
  if (nrow(cross_chr) > 0L) {
    abort(
      "version-stripped GTF has %d gene IDs mapping to >1 autosome; first examples: %s",
      nrow(cross_chr),
      paste(head(cross_chr$gene, 10L), collapse = ", ")
    )
  }

  g <- g[
    ,
    .(
      CHROM = first(CHROM),
      gene_start = min(gene_start),
      gene_end = max(gene_end)
    ),
    by = gene
  ]
  g[, gene_length := as.integer(gene_end - gene_start + 1L)]

  expected <- colnames(numbat::gtf_hg38)
  required <- c("CHROM", "gene", "gene_length", "gene_start", "gene_end")

  if (!setequal(expected, required)) {
    abort(
      "unexpected Numbat GTF schema in installed package: [%s]",
      paste(expected, collapse = ", ")
    )
  }

  g <- as.data.frame(g)[, expected, drop = FALSE]
  rownames(g) <- g$gene
  g
}

prepare_segments <- function(path, mode, build, coordinates) {
  assert_file(path, paste0("--", mode))

  if (is_missing(build) || build != "GRCh38") {
    abort(
      "%s requires explicit --segments_build GRCh38; got: %s",
      mode, ifelse(is_missing(build), "<missing>", build)
    )
  }

  if (is_missing(coordinates) || coordinates != "1-based-closed") {
    abort(
      paste0(
        "%s requires explicit --segments_coordinates 1-based-closed. ",
        "Do not pass BED-style 0-based half-open coordinates without conversion."
      ),
      mode
    )
  }

  s <- data.table::fread(path, showProgress = FALSE)

  if (mode == "segs_loh") {
    req <- c("CHROM", "seg", "seg_start", "seg_end")
  } else if (mode == "segs_fix") {
    req <- c("CHROM", "seg", "seg_start", "seg_end", "cnv_state")
  } else {
    abort("internal error: unknown segment mode %s", mode)
  }

  miss <- setdiff(req, colnames(s))
  if (length(miss) > 0L) {
    abort("%s is missing required columns: %s", mode, paste(miss, collapse = ", "))
  }

  s[, CHROM := normalize_chrom(CHROM)]

  non_auto <- sum(!s$CHROM %in% as.character(1:22))
  if (non_auto > 0L) {
    log_msg("%s: dropping %d non-autosomal rows (Numbat CNV workflow uses chr1-22)", mode, non_auto)
  }
  s <- s[CHROM %in% as.character(1:22)]

  s[, seg_start := suppressWarnings(as.integer(seg_start))]
  s[, seg_end := suppressWarnings(as.integer(seg_end))]

  if (nrow(s) == 0L) abort("%s has no autosomal rows after normalization", mode)
  if (anyNA(s$seg_start) || anyNA(s$seg_end)) abort("%s has NA/non-integer coordinates", mode)
  if (any(s$seg_start < 1L) || any(s$seg_end < s$seg_start)) {
    abort("%s contains invalid 1-based closed intervals", mode)
  }

  seg_raw <- as.character(s$seg)
  if (anyNA(seg_raw) || any(!nzchar(seg_raw))) abort("%s contains NA/empty segment IDs", mode)
  if (anyDuplicated(seg_raw)) {
    s[, seg := paste0(CHROM, "_", seg_raw)]
  } else {
    s[, seg := seg_raw]
  }
  if (anyDuplicated(s$seg)) abort("%s segment IDs are not unique after normalization", mode)

  if (mode == "segs_fix") {
    s[, cnv_state := tolower(as.character(cnv_state))]
    valid <- c("neu", "amp", "del", "loh", "bamp", "bdel")
    bad <- setdiff(unique(s$cnv_state), valid)
    if (length(bad) > 0L) {
      abort(
        "invalid cnv_state value(s) in --segs_fix: %s; allowed: %s",
        paste(bad, collapse = ", "),
        paste(valid, collapse = ", ")
      )
    }
  }

  data.table::setorder(s, CHROM, seg_start, seg_end)
  as.data.frame(s)
}

assert_file(opt$counts, "--counts")
assert_file(opt$df_allele, "--df_allele")
assert_file(opt$gtf, "--gtf")

if (is.na(opt$ncores) || opt$ncores < 1L) abort("--ncores must be >= 1")
detected_cores <- parallel::detectCores(logical = TRUE)
if (!is.na(detected_cores) && opt$ncores > detected_cores) {
  log_msg("WARNING: --ncores=%d exceeds detected logical cores=%d", opt$ncores, detected_cores)
}

if (!is.finite(opt$ref_frac) || opt$ref_frac <= 0 || opt$ref_frac >= 1) {
  abort("--ref_frac must be strictly between 0 and 1")
}
if (is.na(opt$seed)) abort("--seed must be an integer")
if (is.na(opt$init_k) || opt$init_k < 1L) abort("--init_k must be >= 1")
if (is.na(opt$min_cells) || opt$min_cells < 1L) abort("--min_cells must be >= 1")
if (is.na(opt$max_iter) || opt$max_iter < 1L) abort("--max_iter must be >= 1")
if (!is.finite(opt$min_gtf_overlap) || opt$min_gtf_overlap <= 0 || opt$min_gtf_overlap > 1) {
  abort("--min_gtf_overlap must be in (0,1]")
}

has_loh <- !is_missing(opt$segs_loh)
has_fix <- !is_missing(opt$segs_fix)

if (has_loh && has_fix) {
  abort("provide only one of --segs_loh or --segs_fix, never both")
}

numbat_version <- utils::packageVersion("numbat")
if (numbat_version < "1.4.2") {
  abort(
    "installed numbat=%s is too old for this script's custom-GTF path; require >=1.4.2",
    as.character(numbat_version)
  )
}

if (dir.exists(opt$outdir)) {
  old_hits <- list.files(
    opt$outdir,
    pattern = "^(clone_post_|joint_post_|exp_post_|allele_post_|tree_final_|segs_consensus_|per_unit_posterior\\.csv)",
    full.names = TRUE
  )
  if (length(old_hits) > 0L && !isTRUE(opt$overwrite)) {
    abort(
      "outdir contains previous Numbat outputs (%d files). Use a new --outdir or --overwrite: %s",
      length(old_hits), opt$outdir
    )
  }
  if (isTRUE(opt$overwrite)) {
    unlink(opt$outdir, recursive = TRUE, force = TRUE)
  }
}
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

raw <- readRDS(opt$counts)
if (!(is.matrix(raw) || inherits(raw, "Matrix"))) {
  abort("--counts RDS must contain a matrix/Matrix object; got class: %s",
        paste(class(raw), collapse = ","))
}

m <- methods::as(Matrix::Matrix(raw, sparse = TRUE), "dgCMatrix")
rm(raw)

if (is.null(rownames(m)) || is.null(colnames(m))) {
  abort("count matrix must have both row names (genes) and column names (units)")
}
if (anyNA(rownames(m)) || any(!nzchar(rownames(m)))) abort("count matrix has NA/empty gene IDs")
if (anyNA(colnames(m)) || any(!nzchar(colnames(m)))) abort("count matrix has NA/empty unit IDs")
if (anyDuplicated(colnames(m))) abort("count matrix has duplicated unit names")

assert_whole_nonnegative(m@x, "count matrix")

gene_id_raw <- rownames(m)
gene_id <- sub("\\.\\d+$", "", gene_id_raw)
n_dup_rows <- length(gene_id) - length(unique(gene_id))

if (n_dup_rows > 0L) {
  log_msg("collapsing %d duplicate count rows created by Ensembl-version stripping", n_dup_rows)
}
m <- collapse_sparse_rows(m, gene_id)

if (anyDuplicated(rownames(m))) abort("internal error: duplicate genes remain after collapse")
if (any(Matrix::colSums(m) <= 0)) {
  bad <- colnames(m)[Matrix::colSums(m) <= 0]
  abort("zero-count target/reference units present: %s", paste(bad, collapse = ", "))
}

cn <- colnames(m)
is_A <- grepl("^C1_LLU_A(?:_|$)", cn, perl = TRUE)
is_B <- grepl("^C1_LLU_B(?:_|$)", cn, perl = TRUE)

if (any(is_A & is_B)) abort("a unit matched both A and B naming rules")
unknown <- cn[!(is_A | is_B)]
if (length(unknown) > 0L) {
  abort(
    "count matrix contains %d units that are neither C1_LLU_A nor C1_LLU_B; first: %s",
    length(unknown), paste(head(unknown, 10L), collapse = ", ")
  )
}

A <- m[, is_A, drop = FALSE]
B <- m[, is_B, drop = FALSE]

if (ncol(A) < 1L) abort("no C1_LLU_A units found")
if (ncol(B) < 3L) abort("need at least 3 C1_LLU_B units for a disjoint reference/scored split")

log_msg("loaded run-level count matrix: genes=%d, A=%d, B=%d", nrow(m), ncol(A), ncol(B))

gtf_custom <- read_exact_gene_gtf(opt$gtf)

genes_shared <- intersect(rownames(m), gtf_custom$gene)
overlap_frac <- length(genes_shared) / nrow(m)

log_msg(
  "count/GTF autosomal Ensembl overlap: %d/%d genes (%.2f%%)",
  length(genes_shared), nrow(m), 100 * overlap_frac
)

if (overlap_frac < opt$min_gtf_overlap) {
  abort(
    "count/GTF overlap %.3f is below --min_gtf_overlap %.3f; likely identifier/reference mismatch",
    overlap_frac, opt$min_gtf_overlap
  )
}

keep_genes <- rownames(m)[rownames(m) %in% genes_shared]
m <- m[keep_genes, , drop = FALSE]
A <- m[, colnames(A), drop = FALSE]
B <- m[, colnames(B), drop = FALSE]
gtf_custom <- gtf_custom[match(keep_genes, gtf_custom$gene), , drop = FALSE]

if (anyNA(gtf_custom$gene)) abort("internal error while ordering custom GTF")
if (!identical(rownames(m), gtf_custom$gene)) abort("count/GTF gene order mismatch after intersection")

saveRDS(gtf_custom, file.path(opt$outdir, "gtf_numbat_ensembl_v84_autosomes.rds"))

RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(opt$seed)

b_all <- sort(colnames(B))
n_ref <- floor(length(b_all) * opt$ref_frac)
n_score <- length(b_all) - n_ref

if (n_ref < 2L) abort("B reference subset would contain only %d unit(s)", n_ref)
if (n_score < 1L) abort("B scored subset would contain no units")

b_ref <- sort(sample(b_all, size = n_ref, replace = FALSE))
b_score <- setdiff(b_all, b_ref)

A_ids <- sort(colnames(A))

unit_roles <- data.table::rbindlist(list(
  data.table(unit = A_ids, source = "A", role = "target", true_label = "tumour"),
  data.table(unit = b_ref, source = "B", role = "reference", true_label = "normal"),
  data.table(unit = b_score, source = "B", role = "target", true_label = "normal")
))
data.table::fwrite(unit_roles, file.path(opt$outdir, "unit_roles.tsv"), sep = "\t")

ledger_entry <- list(
  step = "B_run_level_split",
  analysis_level = "run_level_technical_pilot",
  seed = opt$seed,
  rng_kind = "Mersenne-Twister/Inversion/Rejection",
  ref_frac = opt$ref_frac,
  n_A_target = length(A_ids),
  n_B_total = length(b_all),
  n_B_reference = length(b_ref),
  n_B_scored = length(b_score),
  B_reference = b_ref,
  B_scored = b_score
)
write_jsonl_one(ledger_entry, file.path(opt$outdir, "seed_ledger.jsonl"))

log_msg(
  "B split seed=%d: B_total=%d, B_reference=%d, B_scored=%d",
  opt$seed, length(b_all), length(b_ref), length(b_score)
)

count_mat <- methods::as(
  cbind(A[, A_ids, drop = FALSE], B[, b_score, drop = FALSE]),
  "dgCMatrix"
)

true_label <- setNames(
  c(rep("tumour", length(A_ids)), rep("normal", length(b_score))),
  c(A_ids, b_score)
)

target_role <- setNames(
  c(rep("A_target", length(A_ids)), rep("B_scored", length(b_score))),
  c(A_ids, b_score)
)

ref_mat <- B[, b_ref, drop = FALSE]
ref_annot <- data.frame(
  cell = b_ref,
  group = rep("HCC1395BL_reference", length(b_ref)),
  stringsAsFactors = FALSE
)

lambdas_ref <- numbat::aggregate_counts(
  count_mat = ref_mat,
  annot = ref_annot,
  normalized = TRUE,
  verbose = FALSE
)

if (!is.matrix(lambdas_ref)) lambdas_ref <- as.matrix(lambdas_ref)
if (anyNA(lambdas_ref)) abort("lambdas_ref contains NA")
if (anyDuplicated(rownames(lambdas_ref))) abort("lambdas_ref contains duplicated genes")
if (!all(rownames(count_mat) %in% rownames(lambdas_ref))) {
  abort("lambdas_ref is missing genes present in count_mat")
}
if (any(abs(colSums(lambdas_ref) - 1) > 1e-8)) {
  abort("lambdas_ref is not normalized to unit column sum")
}

saveRDS(lambdas_ref, file.path(opt$outdir, "lambdas_ref_B_only.rds"))

n_target <- ncol(count_mat)
log_msg(
  "Numbat target units=%d (A=%d + scored-B=%d); reference-B=%d",
  n_target, length(A_ids), length(b_score), length(b_ref)
)

if (opt$min_cells >= n_target) {
  abort("--min_cells=%d is >= target unit count=%d", opt$min_cells, n_target)
}
if (opt$init_k > n_target) {
  abort("--init_k=%d exceeds target unit count=%d", opt$init_k, n_target)
}

df_allele <- data.table::fread(opt$df_allele, showProgress = FALSE)

allele_req <- c("cell", "snp_id", "CHROM", "POS", "cM", "REF", "ALT", "AD", "DP", "GT")
allele_missing <- setdiff(allele_req, colnames(df_allele))
if (length(allele_missing) > 0L) {
  abort("df_allele missing required columns: %s", paste(allele_missing, collapse = ", "))
}

df_allele[, cell := as.character(cell)]
df_allele[, snp_id := as.character(snp_id)]
df_allele[, CHROM := normalize_chrom(CHROM)]
df_allele[, POS := suppressWarnings(as.integer(POS))]
df_allele[, AD := suppressWarnings(as.numeric(AD))]
df_allele[, DP := suppressWarnings(as.numeric(DP))]
df_allele[, GT := as.character(GT)]

if (anyNA(df_allele$cell) || any(!nzchar(df_allele$cell))) abort("df_allele has NA/empty cell/unit values")
if (anyNA(df_allele$snp_id) || any(!nzchar(df_allele$snp_id))) abort("df_allele has NA/empty snp_id values")
if (anyNA(df_allele$POS) || any(df_allele$POS < 1L)) abort("df_allele has invalid POS values")
if (anyNA(df_allele$GT) || any(!nzchar(df_allele$GT))) abort("df_allele has NA/empty GT values")

assert_whole_nonnegative(df_allele$AD, "df_allele$AD")
assert_whole_nonnegative(df_allele$DP, "df_allele$DP")
if (any(df_allele$AD > df_allele$DP)) abort("df_allele contains AD > DP")

n_nonauto <- sum(!df_allele$CHROM %in% as.character(1:22))
if (n_nonauto > 0L) {
  log_msg("df_allele: dropping %d non-autosomal rows", n_nonauto)
}
df_allele <- df_allele[CHROM %in% as.character(1:22)]

df_allele <- df_allele[cell %in% colnames(count_mat)]

if (nrow(df_allele) == 0L) {
  abort("no allele rows match A + scored-B target units")
}

dup_cell_snp <- df_allele[, .N, by = .(cell, snp_id)][N > 1L]
if (nrow(dup_cell_snp) > 0L) {
  abort(
    "df_allele has %d duplicated cell/snp_id pairs; first: %s/%s",
    nrow(dup_cell_snp), dup_cell_snp$cell[1L], dup_cell_snp$snp_id[1L]
  )
}

gt_conflict <- df_allele[
  !is.na(GT) & nzchar(GT),
  .(n_GT = uniqueN(GT)),
  by = snp_id
][n_GT > 1L]
if (nrow(gt_conflict) > 0L) {
  abort(
    "inconsistent phased genotype across units at %d SNPs; first: %s",
    nrow(gt_conflict), gt_conflict$snp_id[1L]
  )
}

allele_coverage <- df_allele[
  ,
  .(
    n_allele_rows = .N,
    n_snps = uniqueN(snp_id),
    total_DP = sum(DP),
    total_AD = sum(AD)
  ),
  by = cell
]

allele_coverage <- merge(
  data.table(cell = colnames(count_mat)),
  allele_coverage,
  by = "cell",
  all.x = TRUE,
  sort = FALSE
)
for (j in c("n_allele_rows", "n_snps", "total_DP", "total_AD")) {
  set(allele_coverage, which(is.na(allele_coverage[[j]])), j, 0)
}
allele_coverage[, true_label := true_label[cell]]
allele_coverage[, target_role := target_role[cell]]
data.table::fwrite(
  allele_coverage,
  file.path(opt$outdir, "allele_coverage_by_unit.tsv"),
  sep = "\t"
)

missing_allele_units <- allele_coverage[n_allele_rows == 0L, cell]
if (length(missing_allele_units) > 0L && !isTRUE(opt$allow_missing_allele_units)) {
  abort(
    paste0(
      "%d/%d target units have zero allele rows. ",
      "This technical allele-arm pilot is gated on non-empty allele evidence per target. ",
      "First missing: %s"
    ),
    length(missing_allele_units), n_target,
    paste(head(missing_allele_units, 10L), collapse = ", ")
  )
}
if (length(missing_allele_units) > 0L) {
  log_msg(
    "WARNING: allowing %d target units with zero allele rows because --allow_missing_allele_units was set",
    length(missing_allele_units)
  )
}

log_msg(
  "allele data after target/autosome filtering: rows=%d, SNPs=%d, represented_units=%d/%d",
  nrow(df_allele),
  data.table::uniqueN(df_allele$snp_id),
  data.table::uniqueN(df_allele$cell),
  n_target
)

seg_args <- list()
segment_mode <- NA_character_
segment_path <- NA_character_

if (has_loh) {
  seg_args$segs_loh <- prepare_segments(
    opt$segs_loh, "segs_loh", opt$segments_build, opt$segments_coordinates
  )
  segment_mode <- "segs_loh"
  segment_path <- opt$segs_loh
  log_msg("LOH mode: validated external segs_loh")
} else if (has_fix) {
  seg_args$segs_consensus_fix <- prepare_segments(
    opt$segs_fix, "segs_fix", opt$segments_build, opt$segments_coordinates
  )
  segment_mode <- "segs_consensus_fix"
  segment_path <- opt$segs_fix
  log_msg("LOH/CNV mode: validated external segs_consensus_fix")
} else {
  seg_args$call_clonal_loh <- TRUE
  segment_mode <- "call_clonal_loh"
  log_msg("LOH mode: fallback call_clonal_loh=TRUE")
}

meta <- list(
  analysis_level = "run_level_technical_pilot",
  reference_build = "GRCh38",
  package_versions = list(
    R = as.character(getRversion()),
    numbat = as.character(utils::packageVersion("numbat")),
    Matrix = as.character(utils::packageVersion("Matrix")),
    data_table = as.character(utils::packageVersion("data.table")),
    dplyr = as.character(utils::packageVersion("dplyr"))
  ),
  parameters = list(
    genome = "hg38",
    gamma = 5,
    ncores = opt$ncores,
    ref_frac = opt$ref_frac,
    seed = opt$seed,
    init_k = opt$init_k,
    min_cells = opt$min_cells,
    max_iter = opt$max_iter,
    min_gtf_overlap = opt$min_gtf_overlap,
    segment_mode = segment_mode,
    segments_build = if (is_missing(opt$segments_build)) NA_character_ else opt$segments_build,
    segments_coordinates = if (is_missing(opt$segments_coordinates)) NA_character_ else opt$segments_coordinates
  ),
  dimensions = list(
    genes_after_collapse = nrow(m),
    genes_used_by_custom_gtf = nrow(count_mat),
    A_run_units = ncol(A),
    B_run_units = ncol(B),
    B_reference_units = length(b_ref),
    B_scored_units = length(b_score),
    target_units = n_target,
    allele_rows = nrow(df_allele),
    allele_snps = data.table::uniqueN(df_allele$snp_id),
    allele_units = data.table::uniqueN(df_allele$cell)
  ),
  input_md5 = list(
    counts = md5_or_na(opt$counts),
    df_allele = md5_or_na(opt$df_allele),
    gtf = md5_or_na(opt$gtf),
    segments = md5_or_na(segment_path)
  ),
  input_paths = list(
    counts = normalizePath(opt$counts, mustWork = TRUE),
    df_allele = normalizePath(opt$df_allele, mustWork = TRUE),
    gtf = normalizePath(opt$gtf, mustWork = TRUE),
    segments = if (is_missing(segment_path)) NA_character_ else normalizePath(segment_path, mustWork = TRUE)
  )
)
write_json_pretty(meta, file.path(opt$outdir, "run_metadata.pre_numbat.json"))

RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(opt$seed)

run_args <- c(
  list(
    count_mat = count_mat,
    lambdas_ref = lambdas_ref,
    df_allele = as.data.frame(df_allele),
    genome = "hg38",
    gamma = 5,
    init_k = opt$init_k,
    min_cells = opt$min_cells,
    max_iter = opt$max_iter,
    ncores = opt$ncores,
    ncores_nni = opt$ncores,
    plot = TRUE,
    out_dir = opt$outdir
  ),
  seg_args
)

log_msg(
  paste0(
    "run_numbat: analysis=RUN-LEVEL TECHNICAL PILOT; genome=hg38; gamma=5; ",
    "targets=%d; init_k=%d; min_cells=%d; max_iter=%d; ncores=%d; mode=%s"
  ),
  n_target, opt$init_k, opt$min_cells, opt$max_iter, opt$ncores, segment_mode
)

status <- do.call(numbat::run_numbat, run_args)

if (!identical(status, "Success")) {
  abort("run_numbat did not complete successfully; returned: %s", paste(status, collapse = " "))
}

clone_files <- list.files(
  opt$outdir,
  pattern = "^clone_post_[0-9]+\\.tsv$",
  full.names = TRUE
)

if (length(clone_files) == 0L) {
  abort("run_numbat returned Success but no clone_post_<iteration>.tsv exists")
}

iter <- suppressWarnings(as.integer(sub(
  "^clone_post_([0-9]+)\\.tsv$", "\\1", basename(clone_files)
)))
iter <- iter[is.finite(iter)]
if (length(iter) == 0L) abort("could not parse Numbat iteration number from clone_post files")
i_final <- max(iter)

clone_path <- file.path(opt$outdir, sprintf("clone_post_%d.tsv", i_final))
tree_path <- file.path(opt$outdir, sprintf("tree_final_%d.rds", i_final))
mut_graph_path <- file.path(opt$outdir, sprintf("mut_graph_%d.rds", i_final))
joint_path <- file.path(opt$outdir, sprintf("joint_post_%d.tsv", i_final))
exp_path <- file.path(opt$outdir, sprintf("exp_post_%d.tsv", i_final))
allele_path <- file.path(opt$outdir, sprintf("allele_post_%d.tsv", i_final))
segs_path <- file.path(opt$outdir, sprintf("segs_consensus_%d.tsv", i_final))

for (p in c(clone_path, tree_path, mut_graph_path, joint_path, exp_path, allele_path, segs_path)) {
  assert_file(p, "required Numbat final-iteration output")
}

clone_post <- data.table::fread(clone_path, showProgress = FALSE)
post_req <- c("cell", "clone_opt", "p_cnv", "p_cnv_x", "p_cnv_y")
post_missing <- setdiff(post_req, colnames(clone_post))
if (length(post_missing) > 0L) {
  abort("clone_post_%d.tsv missing columns: %s", i_final, paste(post_missing, collapse = ", "))
}
if (anyDuplicated(clone_post$cell)) abort("clone_post_%d.tsv has duplicated cell/unit rows", i_final)

for (pcol in c("p_cnv", "p_cnv_x", "p_cnv_y")) {
  v <- clone_post[[pcol]]
  if (anyNA(v) || any(!is.finite(v))) abort("%s contains NA/non-finite posterior values", pcol)
  if (any(v < -1e-12 | v > 1 + 1e-12)) abort("%s contains values outside [0,1]", pcol)
}

missing_post <- setdiff(colnames(count_mat), clone_post$cell)
extra_post <- setdiff(clone_post$cell, colnames(count_mat))

if (length(missing_post) > 0L) {
  abort(
    "Numbat final clone posterior is missing %d target units; first: %s",
    length(missing_post), paste(head(missing_post, 10L), collapse = ", ")
  )
}
if (length(extra_post) > 0L) {
  abort(
    "Numbat final clone posterior contains %d unexpected units; first: %s",
    length(extra_post), paste(head(extra_post, 10L), collapse = ", ")
  )
}

clone_post <- clone_post[match(colnames(count_mat), cell)]

cp <- clone_post[
  ,
  .(
    unit = cell,
    clone_opt = clone_opt,
    p_cnv = p_cnv,
    p_cnv_x = p_cnv_x,
    p_cnv_y = p_cnv_y,
    compartment_opt = if ("compartment_opt" %in% names(clone_post)) compartment_opt else NA_character_
  )
]
cp[, true_label := true_label[unit]]
cp[, target_role := target_role[unit]]
cp[, analysis_level := "run_level_technical_pilot"]
cp[, numbat_iteration := i_final]

if (anyNA(cp$true_label) || anyNA(cp$target_role)) {
  abort("failed to attach known A/B labels to final posterior")
}

data.table::fwrite(
  cp,
  file.path(opt$outdir, "per_unit_posterior.csv"),
  sep = ","
)

file.copy(tree_path, file.path(opt$outdir, "gtree.rds"), overwrite = TRUE)
file.copy(mut_graph_path, file.path(opt$outdir, "mut_graph.rds"), overwrite = TRUE)
file.copy(segs_path, file.path(opt$outdir, "segs_consensus_final.tsv"), overwrite = TRUE)

summary_out <- list(
  status = "Success",
  analysis_level = "run_level_technical_pilot",
  final_iteration = i_final,
  n_posterior_units = nrow(cp),
  n_tumour_A = sum(cp$true_label == "tumour"),
  n_normal_B_scored = sum(cp$true_label == "normal"),
  median_p_cnv_A = unname(stats::median(cp$p_cnv[cp$true_label == "tumour"])),
  median_p_cnv_B = unname(stats::median(cp$p_cnv[cp$true_label == "normal"])),
  fraction_A_p_cnv_ge_0_5 = unname(mean(cp$p_cnv[cp$true_label == "tumour"] >= 0.5)),
  fraction_B_p_cnv_lt_0_5 = unname(mean(cp$p_cnv[cp$true_label == "normal"] < 0.5))
)
write_json_pretty(summary_out, file.path(opt$outdir, "run_summary.json"))

meta$output <- list(
  final_iteration = i_final,
  per_unit_posterior_md5 = md5_or_na(file.path(opt$outdir, "per_unit_posterior.csv")),
  tree_final_md5 = md5_or_na(tree_path),
  mut_graph_md5 = md5_or_na(mut_graph_path),
  segs_consensus_md5 = md5_or_na(segs_path)
)
write_json_pretty(meta, file.path(opt$outdir, "run_metadata.json"))

log_msg(
  paste0(
    "SUCCESS: wrote per_unit_posterior.csv (%d RUN-LEVEL units), ",
    "gtree.rds, mut_graph.rds, segs_consensus_final.tsv, run_summary.json"
  ),
  nrow(cp)
)
