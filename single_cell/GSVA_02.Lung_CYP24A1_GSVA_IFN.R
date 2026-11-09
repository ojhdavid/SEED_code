# R 4.3.1
# GSVA >= 1.50 (parameter-object interface)

library(GSVA)
library(org.Hs.eg.db)
library(dplyr)
library(tibble)
library(data.table)
library(pheatmap)
library(optparse)


option_list <- list(
  make_option(c("-e", "--expr_file"), type = "character",
              default = "../Data/LuCA_Tumor_cluster_expr_vst.tsv",
              help = "VST expression matrix (gene x cluster or gene x sample-cluster)", metavar = "character"),
  make_option(c("-g", "--geneset_dir"), type = "character",
              default = "../GSEA/IFN_11_Genesets",
              help = "directory of IFN gene set files", metavar = "character"),
  make_option(c("-k", "--target_cluster"), type = "character", default = "19",
              help = "cluster shown first in the heatmap", metavar = "character"),
  make_option(c("-f", "--figure_file"), type = "character",
              default = "../Figures/GSVA/IFN/GSVA_IFN_11_ISG_Genesets.pdf",
              help = "output heatmap", metavar = "character"),
  make_option(c("-W", "--width"), type = "numeric", default = 14,
              help = "PDF width (inches)", metavar = "numeric"),
  make_option(c("-H", "--height"), type = "numeric", default = 10,
              help = "PDF height (inches)", metavar = "numeric")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Input/output paths are resolved from the caller's working directory.

ISG_geneset <- c("IRF7", "SAMHD1", "ZNFX1", "STAT2", "UNC93B1", "TRIM22", "PARP10", "TRIM21",
                 "IFIT5", "HERC6", "PARP9", "DDX60", "HERC5", "IFI44L", "MX2", "EPSTI1", "IFI27",
                 "MX1", "IFIT2", "RSAD2", "CMPK2", "RIGI", "USP18", "IFIH1", "IFI35", "LAP3",
                 "HELZ2", "PARP14", "DTX3L", "SHISA5", "SMAD9", "TAP1", "OAS1", "OAS2", "IFIT1",
                 "IFI44", "PLSCR1", "IFI6", "ISG15", "IFITM1", "IFITM3", "LY6E", "IFIT3", "UBE2L6",
                 "OAS3", "STAT1", "HLA-B", "HLA-E", "BST2", "B2M")


get_cluster_id <- function(col_names) {
  # "g19" / "Cluster_19" (cluster-level) or "<sample>_19" (sample-cluster level) -> "19"
  ifelse(grepl("^g[0-9]+$", col_names),
         sub("^g", "", col_names),
         sub("^.*_([0-9]+)$", "\\1", col_names))
}


load_expression <- function(expr_file, target_cluster) {
  expr <- fread(expr_file, data.table = FALSE)
  colnames(expr)[1] <- "GENE"
  expr <- column_to_rownames(expr, var = "GENE")
  
  cluster_id <- get_cluster_id(colnames(expr))
  if (any(!grepl("^[0-9]+$", cluster_id))) {
    stop("Cluster ID could not be parsed from column names: ",
         paste(head(colnames(expr)[!grepl("^[0-9]+$", cluster_id)]), collapse = ", "))
  }
  cluster_num <- as.numeric(cluster_id)
  
  # cluster-level matrix: rename columns to Cluster_<id>
  if (all(grepl("^g[0-9]+$", colnames(expr)))) {
    colnames(expr) <- paste0("Cluster_", cluster_id)
  }
  
  # target cluster first, then the remaining clusters in numeric order
  is_target <- cluster_num == as.numeric(target_cluster)
  new_order <- order(!is_target, cluster_num)
  expr <- expr[, new_order]
  attr(expr, "cluster") <- cluster_id[new_order]
  expr
}


load_genesets <- function(geneset_dir) {
  paths <- list.files(geneset_dir, pattern = "\\.tsv$", full.names = TRUE)
  set_names <- basename(paths) %>%
    sub("\\.v[0-9.]+\\.Hs\\.tsv$", "", .) %>%
    gsub("[^A-Za-z0-9_]", "_", .)
  
  # gene symbols are stored in the header line of each file
  genesets <- lapply(paths, function(p) {
    colnames(read.table(p, sep = ",", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE))
  })
  names(genesets) <- set_names
  genesets
}


symbols_to_entrez <- function(symbols) {
  symbols <- unique(symbols)
  symbols <- symbols[!is.na(symbols) & symbols != ""]
  entrez <- unname(mapIds(org.Hs.eg.db, keys = symbols, keytype = "SYMBOL",
                          column = "ENTREZID", multiVals = "first"))
  unique(entrez[!is.na(entrez) & entrez != ""])
}


expr_to_entrez_matrix <- function(expr) {
  entrez_id <- mapIds(org.Hs.eg.db, keys = rownames(expr), keytype = "SYMBOL",
                      column = "ENTREZID", multiVals = "first")
  keep <- !is.na(entrez_id) & !duplicated(entrez_id)
  mat <- as.matrix(expr[keep, , drop = FALSE])
  rownames(mat) <- entrez_id[keep]
  mat
}


GSVA_IFN <- function(expr_file, geneset_dir, target_cluster, figure_file, width, height) {
  if (!dir.exists(dirname(figure_file))) {
    dir.create(dirname(figure_file), recursive = TRUE)
    cat("Created directory:", dirname(figure_file), "\n")
  } else {
    cat("Using directory:", dirname(figure_file), "\n")
  }
  
  expr <- load_expression(expr_file, target_cluster)
  cluster_id <- attr(expr, "cluster")
  expr_matrix <- expr_to_entrez_matrix(expr)
  
  genesets <- c(load_genesets(geneset_dir), list(ISG_geneset = ISG_geneset))
  genesets <- lapply(genesets, symbols_to_entrez)
  
  param <- gsvaParam(exprData = expr_matrix,
                     geneSets = genesets,
                     kcdf = "Gaussian",   # VST-transformed data
                     minSize = 3,
                     maxSize = 500,
                     maxDiff = TRUE)
  GSVA_results <- gsva(param, verbose = FALSE)
  print(dim(GSVA_results))
  
  # annotate columns by cluster when the matrix is at sample-cluster level
  annotation_col <- NA
  if (length(unique(cluster_id)) < ncol(GSVA_results)) {
    annotation_col <- data.frame(Cluster = factor(cluster_id, levels = unique(cluster_id)),
                                 row.names = colnames(GSVA_results))
  }
  
  # Fit the heatmap width to the number of columns.
  cell_width <- min(25, 180 * 72 / ncol(GSVA_results))
  
  # Write the heatmap PDF with pheatmap.
  pheatmap::pheatmap(GSVA_results,
                     filename = figure_file,
                     width = width,
                     height = height,
                     annotation_col = annotation_col,
                     show_colnames = cell_width >= 8,
                     fontsize_row = 8,
                     fontsize_col = 8,
                     cluster_cols = FALSE,
                     cluster_rows = FALSE,
                     cellwidth = cell_width,
                     cellheight = 25,
                     angle_col = "45")
  
  return(GSVA_results)
}


GSVA_results <- GSVA_IFN(expr_file = opt$expr_file, geneset_dir = opt$geneset_dir,
                         target_cluster = opt$target_cluster, figure_file = opt$figure_file,
                         width = opt$width, height = opt$height)
