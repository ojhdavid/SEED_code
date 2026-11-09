library(optparse)
library(readr)
library(dplyr)

option_list <- list(
  make_option(c("-s", "--sample"), type = "character", help = "Visium sample"),
  make_option(c("-a", "--annotation_csv"), type = "character", help = "Pathologist annotation CSV"),
  make_option(c("-c", "--card_rds"), type = "character", help = "CARD object RDS"),
  make_option(c("-o", "--output"), type = "character", default = NULL, help = "Merged CSV output")
)

opt <- parse_args(OptionParser(option_list = option_list))

read_pathology_annotation <- function(path, sample) {
  annotation <- read_csv(path, show_col_types = FALSE)
  barcode_col <- if ("Barcode" %in% colnames(annotation)) "Barcode" else colnames(annotation)[1]
  value_col <- paste0("Annotation_professor_", sample)
  if (!value_col %in% colnames(annotation)) {
    candidates <- grep("^Annotation", colnames(annotation), value = TRUE)
    if (length(candidates) == 0) {
      stop("No pathologist annotation column found.")
    }
    value_col <- candidates[1]
  }
  annotation %>%
    transmute(barcode = .data[[barcode_col]], pathology_annotation = .data[[value_col]])
}

read_card_argmax <- function(path) {
  card_obj <- readRDS(path)
  prop <- as.data.frame(card_obj@Proportion_CARD)
  prop$barcode <- rownames(prop)
  prop$card_max_proportion <- apply(prop[, setdiff(colnames(prop), "barcode"), drop = FALSE], 1, max)
  prop$card_argmax <- colnames(prop[, setdiff(colnames(prop), c("barcode", "card_max_proportion")), drop = FALSE])[
    apply(prop[, setdiff(colnames(prop), c("barcode", "card_max_proportion")), drop = FALSE], 1, which.max)
  ]
  prop %>% select(barcode, card_argmax, card_max_proportion, everything())
}

main <- function(sample, annotation_csv, card_rds, output) {
  if (is.null(sample) || !grepl("^LUAD-(SMC|UHC)-", sample)) {
    stop("--sample must be a mapped public sample ID, e.g. LUAD-SMC-103P")
  }
  if (is.null(annotation_csv) || !file.exists(annotation_csv)) stop("Missing --annotation_csv")
  if (is.null(card_rds) || !file.exists(card_rds)) stop("Missing --card_rds")
  if (is.null(output)) {
    dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
    output <- file.path("results/tables", paste0(sample, "_pathology_card_annotation.csv"))
  }
  merged <- read_pathology_annotation(annotation_csv, sample) %>%
    full_join(read_card_argmax(card_rds), by = "barcode")
  merged$public_sample_id <- sample
  merged$spot_id <- paste(sample, merged$barcode, sep = "::")
  if (anyNA(merged$barcode) || anyDuplicated(merged$spot_id)) {
    stop("Missing or duplicated spot barcode in merged CARD/pathology table")
  }
  write_csv(merged, output)
  print(table(merged$pathology_annotation, useNA = "ifany"))
}

main(opt$sample, opt$annotation_csv, opt$card_rds, opt$output)
