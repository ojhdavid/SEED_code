#!/usr/bin/env Rscript
if ("--help" %in% commandArgs(trailingOnly=TRUE)) {
  cat("Usage: Rscript visium/03_annotation/convert_luca_counts.R --counts-h5 FILE --metadata FILE.csv --output-rds FILE.rds\n")
  quit(status=0)
}
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(optparse)})
opt <- parse_args(OptionParser(option_list=list(
  make_option("--counts-h5", dest="counts_h5", type="character"),
  make_option("--metadata", type="character"),
  make_option("--output-rds", dest="output_rds", type="character"))))
if (is.null(opt$counts_h5) || is.null(opt$metadata) || is.null(opt$output_rds)) stop("All arguments required")
if (file.exists(opt$output_rds)) stop("Refusing to overwrite reference RDS")
counts <- Read10X_h5(opt$counts_h5, use.names=TRUE, unique.features=FALSE)
metadata <- read.csv(opt$metadata, stringsAsFactors=FALSE, check.names=FALSE)
if (!inherits(counts, "dgCMatrix") || anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) stop("Expected unique named gene-by-cell sparse matrix")
if (!identical(colnames(counts), metadata$cellID)) stop("Counts/metadata cell order mismatch")
saveRDS(counts, opt$output_rds)
cat("Reference RDS:", nrow(counts), "genes", ncol(counts), "cells\n")
