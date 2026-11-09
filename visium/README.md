# Visium

The SMC-103 workflow computes case embeddings, SEED groups, a gene signature,
spatial/UMAP plots and meta-scores. It requires the full SMC input cohort;
prepare the [input files](../docs/inputs.md) before running.

## Setup and run

Use Python **3.10.12**. The installer reads the
[package lock](06_downstream/environment-raw-signature.lock.txt); the runner
sets the numerical runtime options. Run from the repository root:

```bash
BASE_PYTHON=/path/to/python3.10 BOOTSTRAP_PYTHON=/path/to/python-with-pip \
  bash visium/06_downstream/install_raw_signature_environment.sh /new/env/visium
cp config/visium.example.json config/visium.local.json
# Set visium_root and output; use a new output directory.
/new/env/visium/bin/python -s visium/run_alk103_workflow.py \
  --config config/visium.local.json
```

Main results are `analysis/derived_signature.tsv`,
`analysis/SEED_early_membership.tsv`, the `case/` H5AD, and figures/scores in
`primary_spatial/`, `seed_umap/`, `early_umap/` and `scoring/`.
Execution logs are in `logs/`.

## Modules

| Directory | Purpose |
|---|---|
| [01 QC](01_qc_preprocessing/README.md) | Seurat QC and SCT counts for CARD |
| [02 Spatial domains](02_spatial_domains/README.md) | SpaGCN/Leiden domains |
| [03 Annotation](03_annotation/README.md) | CARD proportions and pathology merge |
| [04 Integration](04_integration/README.md) | Whole-cohort and matched-case UMAPs |
| [05 SEED](05_seed_definition/README.md) | SEED definition and case thresholds |
| [06 Downstream](06_downstream/README.md) | Signature derivation and scoring |

Integration reads Space Ranger counts and manual annotations directly. QC/CARD
and SpaGCN are companion analyses. Each module lists its environment and entry point.
