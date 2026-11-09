#!/usr/bin/env Rscript
# TCGA-LUAD expression, quartile-group Cox and Kaplan-Meier analyses.
# Input: saved GDCprepare SummarizedExperiment for tumor and normal samples.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || '--help' %in% args) {
  cat('Usage: Rscript bulk_rnaseq/02_tcga_survival.R GDCprepare_object.rds output_dir\n')
  quit(status = if ('--help' %in% args) 0 else 1)
}
if (file.exists(args[2])) stop('Choose a new output directory')
input_file <- normalizePath(args[1], mustWork = TRUE)
dir.create(args[2], recursive = TRUE)
output_dir <- normalizePath(args[2], mustWork = TRUE)
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(limma)
  library(org.Hs.eg.db)
  library(dplyr)
  library(survival)
  library(survminer)
  library(edgeR)
  library(glue)
  library(purrr)
  library(ggplot2)
})
writeLines(capture.output(sessionInfo()), file.path(output_dir, 'sessionInfo.txt'))
data <- readRDS(input_file)
if (!methods::is(data, 'SummarizedExperiment')) stop('Input must be a GDCprepare SummarizedExperiment')
required <- c('days_to_death', 'days_to_last_follow_up', 'vital_status', 'bcr_patient_barcode')
if (!all(required %in% names(colData(data)))) stop('Input lacks required clinical fields')
if (!'unstranded' %in% assayNames(data) || assayNames(data)[1] != 'unstranded') {
  stop('The first assay must be unstranded raw counts')
}
rna <- assay(data)
clinical <- colData(data) %>% as.data.frame()

rna <- rna[rowSums(rna == 0) <= (ncol(rna) / 2), ]

group <- factor(ifelse(substr(colnames(rna), 14, 14) == "0", "Tumor", "Normal"))
design <- model.matrix(~ group)
dge <- edgeR::DGEList(counts = rna)
dge <- edgeR::calcNormFactors(dge)
voom_res <- limma::voom(dge, design, plot = FALSE)
rna_vm <- voom_res$E

gene_ids <- gsub("\\..*", "", rownames(rna_vm))
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = gene_ids, column = "SYMBOL",
                       keytype = "ENSEMBL", multiVals = "first")
valid_idx <- !is.na(gene_symbols)
rna_vm <- rna_vm[valid_idx, ]
rownames(rna_vm) <- unname(gene_symbols[valid_idx])

clinical$OS <- ifelse(is.na(clinical$days_to_death),
                      clinical$days_to_last_follow_up,
                      clinical$days_to_death)
clinical$OS[clinical$OS == 0] <- NA
clinical$event <- ifelse(clinical$vital_status == "Alive", 0, 1)
clinical_df <- clinical[, c("bcr_patient_barcode", "OS", "event")]

run_univariate_cox_q <- function(gene_symbol) {
  if (!(gene_symbol %in% rownames(rna_vm))) return(NULL)
  expr <- rna_vm[gene_symbol, ]
  q1 <- quantile(expr, 0.25)
  q3 <- quantile(expr, 0.75)

  status <- rep(NA, length(expr))
  status[expr <= q1] <- "low"
  status[expr >= q3] <- "high"
  names(status) <- names(expr)
  status <- status[!is.na(status)]

  df <- data.frame(barcode = names(status), group = status)
  df <- merge(df, clinical_df, by.x = "barcode", by.y = "bcr_patient_barcode")
  if (length(unique(df$group)) < 2) return(NULL)

  df$group <- factor(df$group, levels = c("low", "high"))
  cox <- coxph(Surv(OS, event) ~ group, data = df)
  summ <- summary(cox)

  tibble(
    gene = gene_symbol,
    HR = round(summ$coefficients[1, "exp(coef)"], 3),
    lower_CI = round(summ$conf.int[,"lower .95"], 3),
    upper_CI = round(summ$conf.int[,"upper .95"], 3),
    pval = round(summ$coefficients[1, "Pr(>|z|)"], 5)
  )
}

tumor_barcodes <- colnames(rna_vm)[substr(colnames(rna_vm), 14, 14) == "0"]
normal_barcodes <- colnames(rna_vm)[substr(colnames(rna_vm), 14, 14) == "1"]

tumor_idx <- which(colnames(rna_vm) %in% tumor_barcodes)
normal_idx <- which(colnames(rna_vm) %in% normal_barcodes)

colnames(rna_vm) <- substr(colnames(rna_vm), 1, 12)

clinical_df$bcr_patient_barcode <- substr(clinical_df$bcr_patient_barcode, 1, 12)

genes_to_test_preselected<-c('TMSB4X', 'ZFP36L2', 'A2M', 'MCL1', 'AHNAK', 'PIM2', 'DERL3', 'TYMP', 'ANXA6', 'LTF', 'PLEKHO1', 'ZEB2', 'WIPF1', 'PLA2G7', 'P2RY10', 'SEC11C', 'SLC8A1', 'TXNDC5', 'CPVL', 'CRIP1', 'COTL1', 'JAML', 'IGFLR1', 'PLEKHO2', 'LAX1', 'KLF4', 'TXNIP', 'COLEC12', 'CCL18', 'LIMD2', 'HCLS1', 'ITM2A', 'RARRES1', 'CALHM6', 'FYB1', 'LIX1L', 'PRKACB', 'TSC22D3', 'ADGRA2', 'ELMO1', 'STK17B', 'CXorf21', 'SLCO2B1', 'BOC', 'LTB', 'FAM49A', 'CYP1B1', 'SFMBT2', 'NNMT', 'PSAT1', 'SPHK1', 'MICAL1', 'SNX20', 'EVI2A', 'NRCAM', 'ARHGAP30', 'ATP2A3', 'NPC2', 'LHFPL6', 'CALD1', 'NXPE3', 'OLFML1', 'PLXNC1', 'F2R', 'PCDH18', 'ABCD2', 'DENND5B', 'ALDH3A1', 'QKI', 'AQP1', 'MPP1', 'WARS', 'CHST11', 'ABI3BP', 'FRZB', 'AKAP12', 'CPNE5', 'TGFB3', 'RCSD1', 'PLEK', 'APBB1IP', 'CADM1', 'FGF7', 'MDFIC', 'STK4', 'GNG2', 'UBE2J1', 'NME2', 'TCF4', 'ANO1', 'CYP24A1', 'BTN3A2', 'ABLIM3')

genes_to_test_fixed <- genes_to_test_preselected
genes_to_test_fixed[genes_to_test_fixed == "FAM49A"] <- "CYRIA"
genes_to_test_fixed[genes_to_test_fixed == "CXorf21"] <- "TASL"
genes_to_test_fixed[genes_to_test_fixed == "WARS"]   <- "WARS1"

cox_results_selected <- purrr::map_dfr(genes_to_test_fixed, run_univariate_cox_q)

significant_genes <- cox_results_selected %>% filter(pval < 0.05) %>% pull(gene)

sig_cox_res <- cox_results_selected %>% filter(pval < 0.05)

sig_cox_res <- sig_cox_res %>%
  mutate(gene = factor(gene, levels = gene[order(HR)]))

write.csv(cox_results_selected, file.path(output_dir, 'cox_results_93.csv'), row.names = FALSE)
write.csv(sig_cox_res, file.path(output_dir, 'significant_genes.csv'), row.names = FALSE)
write.csv(data.frame(gene = genes_to_test_fixed, present = genes_to_test_fixed %in% rownames(rna_vm)),
          file.path(output_dir, 'gene_coverage.csv'), row.names = FALSE)
write.csv(data.frame(expression_columns = ncol(rna_vm), clinical_rows = nrow(clinical_df),
                     unique_expression_patients = length(unique(colnames(rna_vm))),
                     unique_clinical_patients = length(unique(clinical_df$bcr_patient_barcode)),
                     tested_genes = nrow(cox_results_selected), significant_genes = nrow(sig_cox_res)),
          file.path(output_dir, 'analysis_summary.csv'), row.names = FALSE)
saveRDS(list(rna_vm = rna_vm, clinical_df = clinical_df, tumor_idx = tumor_idx,
             cox_results = cox_results_selected), file.path(output_dir, 'analysis_objects.rds'))
if (anyDuplicated(clinical_df$bcr_patient_barcode)) {
  warning('Clinical IDs are repeated: the barcode merge can contribute multiple rows per patient. See analysis_summary.csv.')
}

forest <- ggplot(sig_cox_res, aes(x = gene, y = HR)) +
  geom_point(color = "#2E86AB", size = 3) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title = "Univariate Cox Regression (Q1 vs Q3)",
       x = "Gene", y = "Hazard Ratio (HR)") +
  theme_bw(base_size = 13)
ggsave(file.path(output_dir, 'cox_forest.pdf'), forest, width = 7, height = max(5, nrow(sig_cox_res)*0.24))

plot_km_gene_q <- function(gene) {
  if (!(gene %in% rownames(rna_vm))) {
    message(glue("⚠️ Gene {gene} not found."))
    return(NULL)
  }

  expr_values <- rna_vm[gene, tumor_idx]
  quantiles <- quantile(expr_values, probs = c(0.25, 0.75), na.rm = TRUE)
  group <- ifelse(expr_values <= quantiles[1], "Low",
                  ifelse(expr_values >= quantiles[2], "High", NA))

  df <- data.frame(
    barcode = substr(names(expr_values), 1, 12),
    group = factor(group, levels = c("Low", "High"))
  )
  df <- df[!is.na(df$group), ]

  df <- merge(df, clinical_df, by.x = "barcode", by.y = "bcr_patient_barcode")

  if (nrow(df) < 10) {
    message(glue("⚠️ Not enough samples for gene {gene} after filtering."))
    return(NULL)
  }

  fit <- survfit(Surv(OS, event) ~ group, data = df)

  ggsurvplot(fit,
             data = df,
             pval = TRUE,
             conf.int = FALSE ,
             risk.table = TRUE,
             xlim = c(0, 1825),
             break.time.by = 365,
             title = paste("Kaplan-Meier for", gene, "(Q1 vs Q3, 5-year)"),
             xlab = "Time (days, up to 5 years)",
             palette = c("#E64B35", "#4DBBD5"))
}

km <- plot_km_gene_q('CYP24A1')
if (!is.null(km)) {
  pdf(file.path(output_dir, 'CYP24A1_KM.pdf'), width = 8, height = 8)
  print(km)
  dev.off()
}
writeLines(capture.output(sessionInfo()), file.path(output_dir, 'sessionInfo.txt'))
cat('Tested:', nrow(cox_results_selected), 'Significant:', nrow(sig_cox_res), '\n')
