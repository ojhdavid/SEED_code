# Single-cell analyses

Start with the tumor-only LuCA H5AD containing raw UMI counts in `layers['count']`
and `sample`, `dataset`, `uicc_stage` metadata. Create the Python environment
with `conda env create -f single_cell/environment.yml`.

```bash
conda activate seed-singlecell
python single_cell/01_Lung_CYP24A1_BBKNN.py \
  --input /data/LuCA_ref_only_tumor_all_stage.h5ad \
  --cluster-metadata /data/LuCA_Ref_Tumor.obs.csv --output-dir /outputs/LuCA
Rscript single_cell/02_Lung_CYP24A1_Seurat.R \
  /outputs/LuCA/LuCA_Ref_Tumor_counts.h5 /outputs/LuCA/LuCA_Ref_Tumor.obs.csv \
  /outputs/LuCA/LUAD_tumor_filtered.rds
```

Preprocessing uses 2,000 HVGs, dataset-based BBKNN, ridge regression and
Leiden clustering (resolution 1.1). Use `--gene-map` for a saved gene mapping.
With `--cluster-metadata`, supplied labels are used downstream;
`--prepare-only` exports counts and metadata without reclustering.

Run `GSVA_01`–`02` with the exported Seurat RDS and metadata CSV. Run
`pySCENIC_01`–`03` with the exported counts H5AD, specifying
`--gene_id_type symbol` in `pySCENIC_01`, and use the same metadata CSV for
regulon plots. The default target is cluster 19 of `leiden_res_1.1`.
Gene sets, pySCENIC motif/ranking resources and raw data are external inputs.
GSVA requires Seurat, DESeq2 and GSVA; regulon inference requires pySCENIC.
