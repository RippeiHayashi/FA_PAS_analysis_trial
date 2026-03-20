#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
samples_file <- args[2]
out_file <- args[3]

cts <- read_tsv(counts_file, show_col_types = FALSE)
samples <- read_tsv(samples_file, show_col_types = FALSE)

count_mat <- as.data.frame(cts)
rownames(count_mat) <- count_mat$pas_id
count_mat$pas_id <- NULL
count_mat <- as.matrix(count_mat)

samples <- samples[match(colnames(count_mat), samples$sample_id), ]
stopifnot(all(samples$sample_id == colnames(count_mat)))

group <- factor(samples$condition)
design <- model.matrix(~ group)

y <- DGEList(counts = count_mat, samples = samples)
keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)

res <- topTags(qlf, n = Inf)$table %>%
  tibble::rownames_to_column("pas_id") %>%
  as_tibble()

write_tsv(res, out_file)
