# CARD and pathology annotation

CARD uses SCT counts from [QC](../01_qc_preprocessing/README.md) and the
[LuCA core reference](https://datasets.cellxgene.cziscience.com/dd538ee7-f5e4-49e9-9f1e-2a1ea5246cf4.h5ad).
Download that asset as `/inputs/luca_core.h5ad`; its SHA-256 is
`c5c78b19c066b32de0422d306a16fcadc3cd23e6ec0fb238611962b9693e2746`.

Run from the repository root. The Python adapter uses the main Visium Python
environment; R conversion and inference use `seed-visium-card`:

```bash
conda env create -f visium/environment-card.yml
conda run -n seed-visium-card Rscript visium/03_annotation/install_card.R
/new/env/visium/bin/python -s visium/03_annotation/prepare_luca_reference.py \
  --input-h5ad /inputs/luca_core.h5ad --output-dir /outputs/card_reference
conda run -n seed-visium-card Rscript visium/03_annotation/convert_luca_counts.R \
  --counts-h5 /outputs/card_reference/CZI_core_reference_counts_10x.h5 \
  --metadata /outputs/card_reference/CZI_core_reference_celltype_all.csv \
  --output-rds /outputs/card_reference/counts.rds
conda run -n seed-visium-card Rscript visium/03_annotation/01_card_deconvolution.R \
  --sample LUAD-SMC-103P --spatial-rds RDS_result/Louvain_LUAD_SMC_103P.rds \
  --reference-counts /outputs/card_reference/counts.rds \
  --reference-metadata /outputs/card_reference/CZI_core_reference_celltype_all.csv \
  --cell_type cell_type_major --output-dir /outputs/card
```

Outputs are a CARD RDS, proportions CSV and coordinate plots. Plot colors are
scaled per cell type; CSV values retain absolute fractions. Use
`02_merge_pathology_card.R --help` for joining pathology annotations to CARD results.
