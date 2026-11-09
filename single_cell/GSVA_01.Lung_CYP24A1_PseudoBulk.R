# R 4.3.1
# Seurat 5.0.1

library(Seurat)
library(tidyverse)
library(DESeq2)
library(optparse)
set.seed(1234)

option_list <- list(
  make_option(c("-o", "--obs_file"), type = "character",
              default = "../Data/LuCA_Ref_Tumor.obs.csv",
              help = "LuCA tumor cell metadata", metavar = "character"),
  make_option(c("-i", "--seurat_file"), type = "character",
              default = "../Data/LUAD_tumor_filtered.rds",
              help = "Seurat object with raw counts", metavar = "character"),
  make_option(c("-c", "--cluster_col"), type = "character", default = "leiden_res_1.1",
              help = "cluster column in metadata", metavar = "character"),
  make_option(c("-m", "--min_counts"), type = "numeric", default = 10,
              help = "minimum total counts per gene", metavar = "numeric"),
  make_option(c("-e", "--expr_file"), type = "character",
              default = "../Data/LuCA_Tumor_cluster_expr_vst.tsv",
              help = "output VST expression matrix (gene x cluster) for GSVA", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Input/output paths are resolved from the caller's working directory.




load_LuCA_counts <- function(obs_file, seurat_file, cluster_col) {
  LuCA_obs <- read_csv(obs_file)
  colnames(LuCA_obs)[1] <- "Cell_ID"

  seurat_obj <- readRDS(seurat_file)
  meta <- data.frame(Cell_ID = colnames(seurat_obj)) %>%
    left_join(LuCA_obs, by = "Cell_ID")
  rownames(meta) <- meta$Cell_ID
  seurat_obj@meta.data <- meta
  seurat_obj$cluster <- as.character(seurat_obj@meta.data[[cluster_col]])

  print(table(seurat_obj$cluster))
  return(seurat_obj)
}


write_gsva_input <- function(pseudobulk, min_counts, expr_file) {
  coldata <- data.frame(Cluster = colnames(pseudobulk), row.names = colnames(pseudobulk))

  # one profile per cluster (no replicates), so VST is used for normalization only
  dds <- DESeqDataSetFromMatrix(countData = pseudobulk, colData = coldata, design = ~ 1)
  dds <- dds[rowSums(counts(dds)) > min_counts, ]
  vsd <- vst(dds, blind = TRUE)

  # VST-normalized gene x cluster matrix for GSVA
  write.table(assay(vsd), file = expr_file, sep = "\t", quote = FALSE, col.names = NA)
}


PseudoBulk_cluster <- function(obs_file, seurat_file, cluster_col, min_counts, expr_file) {
  seurat_obj <- load_LuCA_counts(obs_file, seurat_file, cluster_col)

  # Pseudo-bulk: sum raw counts per cluster
  pseudobulk <- AggregateExpression(seurat_obj,
                                    group.by = "cluster",
                                    assays = "RNA",
                                    slot = "counts",
                                    return.seurat = FALSE)$RNA

  # column names are "19" or "g19" depending on the Seurat version -> "Cluster_19"
  colnames(pseudobulk) <- paste0("Cluster_", sub("^g", "", colnames(pseudobulk)))
  pseudobulk <- round(as.matrix(pseudobulk))

  write_gsva_input(pseudobulk, min_counts, expr_file)
  return(pseudobulk)
}


pseudobulk <- PseudoBulk_cluster(obs_file = opt$obs_file, seurat_file = opt$seurat_file,
                                 cluster_col = opt$cluster_col, min_counts = opt$min_counts,
                                 expr_file = opt$expr_file)
