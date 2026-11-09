#!/usr/bin/env Rscript
# Convert prepared gene-symbol raw counts to the Seurat input for GSVA.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L || '--help' %in% args) {
  cat('Usage: Rscript 02_Lung_CYP24A1_Seurat.R counts.h5 metadata.csv output.rds\n')
  quit(status = if ('--help' %in% args) 0 else 1)
}
suppressPackageStartupMessages(library(Seurat))
if (file.exists(args[3])) stop('Choose a new output file')
counts <- Read10X_h5(args[1], use.names = TRUE, unique.features = FALSE)
metadata <- read.csv(args[2], row.names = 1, check.names = FALSE)
if (anyDuplicated(colnames(counts)) || anyDuplicated(rownames(counts)) ||
    anyDuplicated(rownames(metadata)) || !setequal(colnames(counts), rownames(metadata))) {
  stop('Unique gene/cell IDs and matching metadata are required')
}
metadata <- metadata[colnames(counts), , drop = FALSE]
object <- CreateSeuratObject(counts = counts, meta.data = metadata,
                             project = 'LuCA_Tumor', min.cells = 0, min.features = 0)
dir.create(dirname(args[3]), recursive = TRUE, showWarnings = FALSE)
saveRDS(object, args[3])
