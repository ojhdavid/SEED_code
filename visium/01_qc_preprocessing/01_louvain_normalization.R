if ("--help" %in% commandArgs(trailingOnly = TRUE)) {
  cat("Usage: Rscript visium/01_qc_preprocessing/01_louvain_normalization.R --sample PUBLIC_SAMPLE_ID [--sample-map TSV] [--visium-root DIR] [--coordinate-root DIR]\n")
  quit(status = 0)
}
library(Seurat)
library(ggplot2)
library(readr)
library(optparse)

option_list <- list(
  make_option(c("-s", "--sample"), type = "character",
              help = "Public Visium sample ID", metavar = "character"),
  make_option(c("--sample-map"), dest = "sample_map", type = "character", default = "visium/resources/sample_id_mapping.tsv"),
  make_option(c("--visium-root"), dest = "visium_root", type = "character", default = "data/visium/Visium"),
  make_option(c("--coordinate-root"), dest = "coordinate_root", type = "character", default = "visium/resources/coordinates"),
  make_option(c("-f", "--nFeature_Spatial_min"), type = "numeric", default = 200,
              help = "minimum number of features", metavar = "numeric"),
  make_option(c("-c", "--nCount_Spatial_min"), type = "numeric", default = 3,
              help = "minimum number of counts", metavar = "numeric"),
  make_option(c("-p", "--percent_mito_max"), type = "numeric", default = 5,
              help = "maximum mitogene percentage", metavar = "numeric"),
  make_option(c("-r", "--resolution"), type = "numeric", default = 0.3,
              help = "resolution for clustering", metavar = "numeric")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

visium_version_check <- function(input_file_path) {
  version <- "Unknown"  # Default to "Unknown"
  if (file.exists(file.path(input_file_path, "spatial", "tissue_positions_list.csv"))) version <- "V1"
  if (file.exists(file.path(input_file_path, "spatial", "tissue_positions.csv"))) version <- "V2"
  version
}

resolve_sample <- function(sample, sample_map, coordinate_root) {
  mapping <- readr::read_tsv(sample_map, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  row <- mapping[mapping$public_sample_id == sample, , drop = FALSE]
  if (nrow(row) != 1) stop("Expected one mapping row for public sample ", sample)
  bounds <- NULL
  if (row$tissue[[1]] != "primary") {
    coordinates <- readr::read_csv(file.path(coordinate_root, row$coordinate_table[[1]]), show_col_types = FALSE)
    match <- coordinates[coordinates$sample_sub == row$coordinate_sample_id[[1]], , drop = FALSE]
    if (nrow(match) != 1 || match$sample[[1]] != row$space_ranger_sample[[1]]) {
      stop("Coordinate core mapping mismatch for ", sample)
    }
    bounds <- match[1, c("min_row", "max_row", "min_col", "max_col")]
  }
  list(run = row$space_ranger_sample[[1]], bounds = bounds)
}


Clustering_result_louvain <- function(sample, sample_map, visium_root, coordinate_root, nFeature_Spatial_min, nCount_Spatial_min,percent_mito_max, resolution){
  # sample = LUAD-SMC-51
  source <- resolve_sample(sample, sample_map, coordinate_root)
  input_file_path <- file.path(visium_root, source$run, "Spatial_matrix")
  safe_sample <- gsub("[^[:alnum:]_.]", "_", sample)
  download_path <- file.path("RDS_result", paste0("Louvain_", safe_sample, ".rds"))

  RDS_result <- "RDS_result"
  if (!dir.exists(RDS_result)) {
    dir.create(RDS_result, recursive = TRUE)
    cat("Created directory:", RDS_result)
  } else {
    cat("Using directory:", RDS_result)
  }

  seurat_obj <- Load10X_Spatial(
  data.dir= input_file_path,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice2",
  filter.matrix = TRUE)
  if (!is.null(source$bounds)) {
    xy <- seurat_obj@images[["slice2"]]@coordinates
    b <- source$bounds
    keep <- rownames(xy)[xy$imagerow >= b$min_row & xy$imagerow <= b$max_row &
                         xy$imagecol >= b$min_col & xy$imagecol <= b$max_col]
    if (!length(keep)) stop("No spots remain after mapped tissue-core crop: ", sample)
    seurat_obj <- subset(seurat_obj, cells = keep)
  }
  set.seed(100)
  seurat_obj@images[["slice2"]]@coordinates[["tissue"]] <- as.integer(seurat_obj@images[["slice2"]]@coordinates[["tissue"]])
  seurat_obj@images[["slice2"]]@coordinates[["row"]] <- as.integer(seurat_obj@images[["slice2"]]@coordinates[["row"]])
  seurat_obj@images[["slice2"]]@coordinates[["col"]] <- as.integer(seurat_obj@images[["slice2"]]@coordinates[["col"]])
  seurat_obj@images[["slice2"]]@coordinates[["imagerow"]] <- as.integer(seurat_obj@images[["slice2"]]@coordinates[["imagerow"]])
  seurat_obj@images[["slice2"]]@coordinates[["imagecol"]] <- as.integer(seurat_obj@images[["slice2"]]@coordinates[["imagecol"]])
  seurat_obj$percent_mito <- PercentageFeatureSet(seurat_obj, pattern="^MT.")
  #seurat_obj_filtering
  seurat_obj <- subset(seurat_obj, subset = nFeature_Spatial > nFeature_Spatial_min & nCount_Spatial > nCount_Spatial_min & percent_mito < percent_mito_max)
  # Normalization
  #builds regularized negative binomial models of gene expression in order to account for technical artifacts while preserving biological variance
  seurat_obj<- SCTransform(seurat_obj, assay = "Spatial", verbose=FALSE)
  seurat_obj <- RunPCA(seurat_obj, npcs = 20, verbose = F)
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  seurat_obj <- FindNeighbors(seurat_obj,reduction='pca', dims = 1:17)
  seurat_obj <- FindClusters(seurat_obj, verbose=FALSE, resolution = resolution)
  seurat_obj <- RunUMAP(seurat_obj, reduction='pca', dims = 1:17)
  saveRDS(seurat_obj, file = download_path)
  Dimplot <- "Dimplot"
  if (!dir.exists(Dimplot)) {
    dir.create(Dimplot, recursive = TRUE)
    cat("Created directory:", Dimplot)
  } else {
    cat("Using directory:", Dimplot)
  }
  version <-visium_version_check(input_file_path)
  print(version)
  Point_size_factor <- ifelse(version == "V1", 1, 3.5)
  colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == "seurat_clusters"] <- "Louvain"
  p1 <- SpatialDimPlot(seurat_obj, group.by = "Louvain", pt.size.factor = Point_size_factor, label = TRUE)
  ggsave(file.path("Dimplot", paste0("Louvain_", safe_sample, "_SpatialDimPlot.png")), plot = p1, width = 5, height = 5, units = "in")

}

Clustering_result_louvain(sample = opt$sample, sample_map = opt$sample_map,
                          visium_root = opt$visium_root, coordinate_root = opt$coordinate_root,
                          nFeature_Spatial_min = opt$nFeature_Spatial_min,
                          nCount_Spatial_min = opt$nCount_Spatial_min,
                          percent_mito_max = opt$percent_mito_max, resolution = opt$resolution)
