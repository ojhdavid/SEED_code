# Seurat QC

Space Ranger counts and mapped TMA crops → QC/normalization → SCT Seurat RDS.
Run from the repository root:

```bash
conda env create -f visium/environment-r.yml
conda run -n seed-visium-r Rscript visium/01_qc_preprocessing/01_louvain_normalization.R \
  --sample LUAD-SMC-103P --visium-root /inputs/Visium --coordinate-root visium/resources/coordinates
```

`coordinate-root` contains the CSV filenames in the [sample map](../resources/sample_id_mapping.tsv).
Output: `RDS_result/Louvain_LUAD_SMC_103P.rds`, used by CARD. See `--help` for options.
