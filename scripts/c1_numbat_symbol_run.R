#!/usr/bin/env Rscript
# C1 run-level Numbat: Ensembl->symbol conversion + run_numbat (numbat 1.4.2)
suppressPackageStartupMessages({library(Matrix); library(data.table); library(numbat); library(dplyr); library(stringr)})
set.seed(100)
m <- as(readRDS("out/c1_counts_runlevel.rds"), "dgCMatrix")
gl <- fread(cmd="grep -P \"\\tgene\\t\" ref/genes.gtf", sep="\t", header=FALSE)
gid <- sub("\\.\\d+$","", str_match(gl$V9, "gene_id \"([^\"]+)\"")[,2])
gsym <- str_match(gl$V9, "gene_name \"([^\"]+)\"")[,2]
map <- setNames(gsym, gid)
rn <- sub("\\.\\d+$","", rownames(m)); sym <- map[rn]
keep <- !is.na(sym) & sym!=""; m <- m[keep,]; rownames(m) <- sym[keep]
m <- as(rowsum(as.matrix(m), rownames(m)), "dgCMatrix")
A <- m[, grep("C1_LLU_A", colnames(m)), drop=FALSE]
B <- m[, grep("C1_LLU_B", colnames(m)), drop=FALSE]
set.seed(100); b_all <- sort(colnames(B))
b_ref <- sort(sample(b_all, floor(length(b_all)*0.5))); b_score <- setdiff(b_all, b_ref)
count_mat <- as(cbind(A, B[, b_score, drop=FALSE]), "dgCMatrix")
lambdas_ref <- Matrix::rowSums(B[, b_ref, drop=FALSE]); lambdas_ref <- lambdas_ref/sum(lambdas_ref)
true_label <- setNames(c(rep("tumour",ncol(A)), rep("normal",length(b_score))), colnames(count_mat))
d <- fread("out/allele/C1_HCC1395_allele_counts.tsv.gz"); d <- d[cell %in% colnames(count_mat)]
dir.create("out/numbat", showWarnings=FALSE, recursive=TRUE)
run_numbat(count_mat=count_mat, lambdas_ref=lambdas_ref, df_allele=as.data.frame(d),
           genome="hg38", gamma=5, init_k=2, min_cells=5, max_iter=2, ncores=16, ncores_nni=16,
           call_clonal_loh=TRUE, plot=TRUE, out_dir="out/numbat")
i <- max(as.integer(str_match(list.files("out/numbat","clone_post_\\d+.tsv"),"clone_post_(\\d+)")[,2]), na.rm=TRUE)
cp <- fread(sprintf("out/numbat/clone_post_%d.tsv", i)); cp$true_label <- true_label[cp$cell]
fwrite(cp[, .(unit=cell, clone_opt, p_cnv, p_cnv_x, p_cnv_y, true_label)], "out/numbat/per_unit_posterior.csv")
