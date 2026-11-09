#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || '--help' %in% args) {
  cat('Usage: Rscript bulk_rnaseq/01_download_tcga.R download_dir output.rds\n')
  quit(status = if ('--help' %in% args) 0 else 1)
}
suppressPackageStartupMessages(library(TCGAbiolinks))
if (file.exists(args[2])) stop('Choose a new output file')
query <- GDCquery(
  project = 'TCGA-LUAD',
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  experimental.strategy = 'RNA-Seq',
  workflow.type = 'STAR - Counts',
  sample.type = c('Primary Tumor', 'Solid Tissue Normal')
)
GDCdownload(query, directory = args[1])
data <- GDCprepare(query, directory = args[1])
dir.create(dirname(args[2]), recursive = TRUE, showWarnings = FALSE)
saveRDS(data, args[2])
saveRDS(query, paste0(args[2], '.query.rds'))
writeLines(capture.output(sessionInfo()), paste0(args[2], '.sessionInfo.txt'))
