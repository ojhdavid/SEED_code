library(Seurat)
library(ggplot2)
library(tidyverse)
library(data.table)
library(optparse)

option_list <- list(
  make_option(c("-m", "--motif_raw_file"), type = "character",
              default = "../Data/pySCENIC/RESULTS/LuCA_Ref_Tumor.motifs.csv",
              help = "raw pySCENIC motif table", metavar = "character"),
  make_option(c("-a", "--auc_file"), type = "character",
              default = "../Data/pySCENIC/RESULTS/LuCA_Ref_Tumor.auc.csv",
              help = "AUCell matrix from pySCENIC", metavar = "character"),
  make_option(c("-d", "--metadata_file"), type = "character",
              default = "../Data/LuCA_Ref_Tumor.obs.csv",
              help = "cell metadata with cluster labels", metavar = "character"),
  make_option(c("-o", "--motif_file"), type = "character",
              default = "../Data/pySCENIC/RESULTS/LuCA_Ref_Tumor.motifs_M.csv",
              help = "cleaned pySCENIC motif table", metavar = "character"),
  make_option(c("-g", "--figure_dir"), type = "character",
              default = "../Figures/pySCENIC",
              help = "output directory for figures", metavar = "character"),
  make_option(c("-t", "--result_dir"), type = "character",
              default = "../Data/pySCENIC/RESULTS",
              help = "output directory for tables", metavar = "character"),
  make_option(c("-c", "--cluster_col"), type = "character", default = "leiden_res_1.1",
              help = "cluster column in metadata", metavar = "character"),
  make_option(c("-k", "--target_cluster"), type = "character", default = "19",
              help = "cluster for top-marker DotPlot", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Input/output paths are resolved from the caller's working directory.


make_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    cat("Created directory:", path, "\n")
  } else {
    cat("Using directory:", path, "\n")
  }
}


fix_motif_header <- function(motif_raw_file, motif_file) {
  column_names <- c("TF", "MotifID", "Enrichment_AUC", "Enrichment_NES",
                    "Enrichment_MotifSimilarityQvalue", "Enrichment_OrthologousIdentity",
                    "Enrichment_Annotation", "Enrichment_Context",
                    "TargetGenes", "Enrichment_RankAtMax")
  # skip the 3-line pySCENIC header and keep all values as character
  motif <- fread(motif_raw_file, skip = 3, header = FALSE, sep = ",", colClasses = "character")
  stopifnot(ncol(motif) == length(column_names))
  setnames(motif, column_names)
  # TargetGenes contains commas, so quote every field
  fwrite(motif, file = motif_file, sep = ",", quote = TRUE)
}


dotplot_theme <- function(p) {
  p +
    scale_color_gradientn(colours = c("#2166AC", "#7EB8D6", "white", "#E38066", "#B2182B")) +
    scale_size(range = c(0.25, 10)) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 0.5, face = "bold", size = 12, colour = "black"),
      axis.text.y = element_text(face = "bold", size = 12, colour = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
}


LuCA_Regulon_analysis <- function(auc_file, metadata_file, motif_file,
                                  figure_dir, result_dir, cluster_col, target_cluster) {
  make_dir(figure_dir)
  make_dir(result_dir)

  # AUCell matrix (cells x regulons) -> Seurat object (regulons x cells)
  auc <- fread(auc_file)
  auc <- column_to_rownames(auc, var = "Cell")
  metadata <- fread(metadata_file)
  print(table(metadata[[cluster_col]]))

  seurat_obj <- CreateSeuratObject(counts = t(auc), project = "LuCA_Tumor")
  seurat_obj@assays$RNA$data <- seurat_obj@assays$RNA$counts

  # Add clinical metadata
  clin_meta <- metadata %>%
    dplyr::select(c("V1", "sample", "uicc_stage", "dataset", all_of(cluster_col)))
  colnames(clin_meta)[1] <- "Cell"
  obj_meta <- rownames_to_column(seurat_obj@meta.data, var = "Cell") %>%
    left_join(clin_meta, by = "Cell")
  rownames(obj_meta) <- obj_meta$Cell
  seurat_obj@meta.data <- obj_meta
  Idents(seurat_obj) <- as.factor(seurat_obj@meta.data[[cluster_col]])

  # Regulon markers per cluster
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(markers, file.path(result_dir, paste0(cluster_col, "_Regulon_markers.csv")), row.names = FALSE)

  # Top 30 regulons of the target cluster
  target_markers <- markers %>%
    filter(cluster == target_cluster) %>%
    arrange(desc(avg_log2FC))
  top_features <- head(target_markers$gene, 30)

  p1 <- dotplot_theme(DotPlot(seurat_obj, features = top_features)) +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5))
  ggsave(file.path(figure_dir, paste0("C", target_cluster, "_DotPlot.pdf")),
         plot = p1, width = 14, height = 8)

  # TF - TargetGene table
  motif <- fread(motif_file)
  LuCA_Tumor_motif_long_2 <- motif %>%
    mutate(TargetGene = str_extract_all(TargetGenes, "\\b[A-Za-z][A-Za-z0-9._-]*\\b")) %>%
    unnest(TargetGene) %>%
    dplyr::select(TF, TargetGene) %>%
    distinct()
  write.csv(LuCA_Tumor_motif_long_2, file.path(result_dir, "LuCA_Tumor_motif_long_2.csv"), row.names = FALSE)

  return(LuCA_Tumor_motif_long_2)
}


fix_motif_header(motif_raw_file = opt$motif_raw_file, motif_file = opt$motif_file)

LuCA_Tumor_motif_long_2 <- LuCA_Regulon_analysis(auc_file = opt$auc_file, metadata_file = opt$metadata_file,
                                                 motif_file = opt$motif_file,
                                                 figure_dir = opt$figure_dir, result_dir = opt$result_dir,
                                                 cluster_col = opt$cluster_col, target_cluster = opt$target_cluster)
