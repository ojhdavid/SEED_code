# Spatial domains

SpaGCN/Leiden analysis of mapped Space Ranger samples and TMA crops.

```bash
conda env create -f visium/02_spatial_domains/environment.yml
conda activate seed-spagcn
PYTHONNOUSERSITE=1 PIP_USER=0 python -m pip install --no-user SpaGCN==1.2.7
python visium/02_spatial_domains/01_spagcn_leiden.py \
  --sample LUAD-SMC-103P --method SpaGCN --resolution 0.3
```

Run from the repository root; set input paths with the options in `--help`.
The example writes `h5ad_result/SpaGCN_LUAD_SMC_103P.h5ad`.
