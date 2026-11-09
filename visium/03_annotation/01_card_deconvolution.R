#!/usr/bin/env Rscript
if ("--help" %in% commandArgs(trailingOnly = TRUE)) {
  cat("Usage: Rscript visium/03_annotation/01_card_deconvolution.R --sample PUBLIC_SAMPLE_ID --spatial-rds FILE --reference-counts FILE.rds --reference-metadata FILE.csv --output-dir DIR [--cell_type cell_type_major]\n")
  quit(status = 0)
}
suppressPackageStartupMessages({library(CARD); library(Seurat); library(Matrix); library(optparse); library(ggplot2)})
options(warn=1) # Show warnings immediately.
options <- list(
  make_option(c("--sample"), type="character"),
  make_option(c("--spatial-rds"), dest="spatial_rds", type="character"),
  make_option(c("--reference-counts"), dest="reference_counts", type="character"),
  make_option(c("--reference-metadata"), dest="reference_metadata", type="character"),
  make_option(c("--output-dir"), dest="output_dir", type="character"),
  make_option(c("--cell_type"), type="character", default="cell_type_major"),
  make_option(c("--reference"), type="character", default="core"))
opt <- parse_args(OptionParser(option_list=options))
for (key in c("sample", "spatial_rds", "reference_counts", "reference_metadata", "output_dir")) {
  if (is.null(opt[[key]]) || !nzchar(opt[[key]])) stop("Missing argument: ", key)
}
if (!grepl("^LUAD-(SMC|UHC)-", opt$sample)) stop("Expected a mapped public sample ID")
if (dir.exists(opt$output_dir) && length(list.files(opt$output_dir, all.files=FALSE))) stop("Refusing to reuse nonempty output directory")
dir.create(opt$output_dir, recursive=TRUE, showWarnings=FALSE)
safe_sample <- gsub("[^[:alnum:]_.]", "_", opt$sample)
seurat_object <- readRDS(opt$spatial_rds)
# CARD uses SCT corrected counts.
spatial_count <- seurat_object@assays$SCT@counts
coords <- GetTissueCoordinates(seurat_object)
if (!all(c("imagerow", "imagecol") %in% names(coords))) stop("Expected imagerow/imagecol tissue coordinates")
coords <- data.frame(y=coords$imagerow, x=coords$imagecol, row.names=rownames(coords))
if (!all(colnames(spatial_count) %in% rownames(coords))) stop("Spatial count/coordinate ID mismatch")
coords <- coords[colnames(spatial_count), , drop=FALSE]
reference_count <- readRDS(opt$reference_counts)
reference_meta <- read.csv(opt$reference_metadata, check.names=FALSE, stringsAsFactors=FALSE)
if (!all(c("cellID", "sampleInfo", opt$cell_type) %in% names(reference_meta))) stop("Missing reference metadata columns")
if (anyDuplicated(reference_meta$cellID) || anyDuplicated(colnames(reference_count))) stop("Duplicate reference cell IDs")
if (!setequal(reference_meta$cellID, colnames(reference_count))) stop("Reference count/metadata cell IDs differ")
reference_meta <- reference_meta[match(colnames(reference_count), reference_meta$cellID), , drop=FALSE]
rownames(reference_meta) <- reference_meta$cellID
if (anyDuplicated(rownames(reference_count))) stop("Reference gene symbols must be unique.")
for (matrix in list(reference_count, spatial_count)) {
  values <- if (inherits(matrix, "sparseMatrix")) matrix@x else as.vector(matrix)
  if (any(!is.finite(values)) || any(values < 0) || any(abs(values-round(values)) > 1e-6)) stop("Expected finite nonnegative count matrix")
}
if (anyNA(reference_meta[, c("sampleInfo", opt$cell_type)])) stop("Missing reference sample/cell-type label")
# Select cell types from the requested annotation column.
selected_types <- unique(as.character(reference_meta[[opt$cell_type]]))
input_summary <- data.frame(sample=opt$sample, spatial_spots=ncol(spatial_count), spatial_genes=nrow(spatial_count),
                         reference_cells=ncol(reference_count), reference_genes=nrow(reference_count),
                         common_genes=length(intersect(rownames(spatial_count), rownames(reference_count))),
                         cell_type_column=opt$cell_type, cell_types=length(selected_types),
                         matrix="SCT counts", minCountGene=100, minCountSpot=5)
write.table(input_summary, file.path(opt$output_dir, "input_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(data.frame(cell_type=selected_types), file.path(opt$output_dir, "selected_cell_types.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
card_object <- createCARDObject(sc_count=reference_count, sc_meta=reference_meta,
                               spatial_count=spatial_count, spatial_location=coords,
                               ct.varname=opt$cell_type, ct.select=selected_types,
                               sample.varname="sampleInfo", minCountGene=100, minCountSpot=5)
card_object <- CARD_deconvolution(CARD_object=card_object)
saveRDS(card_object, file.path(opt$output_dir, paste0("CARD_obj_", safe_sample, "_", opt$reference, ".rds")))
proportions <- as.data.frame(card_object@Proportion_CARD, check.names=FALSE)
if (any(!is.finite(as.matrix(proportions)))) stop("Nonfinite CARD proportions")
barcode <- rownames(proportions)
table <- data.frame(canonical_spot_id=paste0(opt$sample, "::", barcode), public_sample_id=opt$sample,
                    barcode=barcode, proportions, check.names=FALSE)
write.csv(table, file.path(opt$output_dir, paste0(safe_sample, "_CARD_proportions.csv")), row.names=FALSE)
# Plot transform: negate the first coordinate (y) only.
plot_coords <- card_object@spatial_location
plot_coords[, 1] <- -plot_coords[, 1]
plot <- CARD.visualize.prop(proportion=card_object@Proportion_CARD, spatial_location=plot_coords,
                           ct.visualize=colnames(card_object@Proportion_CARD),
                           colors=c("lightblue", "lightyellow", "red"), NumCols=4)
# CARD's native display rescales each cell type to its own min/max, not absolute fractions.
plot <- plot + facet_wrap(~Cell_Type, ncol=4, labeller=label_wrap_gen(width=18)) +
  labs(title=paste(opt$sample, "CARD coordinate-only display"), colour="Within-type\nmin-max") +
  theme(strip.text=element_text(size=9, face="bold"), plot.title=element_text(size=12),
        plot.background=element_rect(fill="white", colour=NA))
ggsave(file.path(opt$output_dir, paste0(safe_sample, "_", opt$reference, "_spatialPlot.png")),
       plot=plot, width=14, height=18, units="in", bg="white")
writeLines(c(paste("CARD", packageVersion("CARD")), "spatial_matrix=SCT counts", "ct.select=unique(reference_meta[[cell_type]])",
             "minCountGene=100; minCountSpot=5", "plot transform: negate y only; no histology overlay; within-type min-max display; table retains absolute fractions",
             paste("reference_counts", normalizePath(opt$reference_counts)), paste("reference_metadata", normalizePath(opt$reference_metadata))),
           file.path(opt$output_dir, "provenance.txt"))
cat("CARD completed:", nrow(proportions), "spots;", ncol(proportions), "cell types\n")
