# TCGA-LUAD bulk RNA-seq

Survival analysis of 93 SEED-associated genes using TCGA-LUAD RNA-seq and
clinical data, with univariate Cox regression and Kaplan–Meier plots.
Expression groups are defined by the lower and upper quartiles.

## Run

Run from the repository root. R dependencies are listed in each script.
The analysis uses a saved GDCprepare RDS object.

```bash
Rscript bulk_rnaseq/01_download_tcga.R /data/GDC /data/TCGA_LUAD.rds
Rscript bulk_rnaseq/02_tcga_survival.R /data/TCGA_LUAD.rds /outputs/TCGA
```

Outputs include Cox regression tables, selected genes (nominal `p < 0.05`)
and survival plots.
