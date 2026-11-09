# Integration

- `01_bbknn/01_build_bbknn_umap.py`: whole-cohort uncorrected/BBKNN UMAPs (3,000 HVGs).
- [Matched-case workflow](02_concat_seed/README.md): primary/LN UMAPs for SEED selection (20,000 HVGs).

Both consume Space Ranger counts, mapped crop coordinates and manual annotations.
Use the companion environment for these standalone commands:

```bash
conda env create -f visium/environment.yml
conda activate seed-visium
python visium/04_integration/01_bbknn/01_build_bbknn_umap.py \
  --visium-root /inputs/Visium --coordinate-table visium/resources/coordinates/TMA_coordinates_except_german.csv \
  --manual-annotation-root visium/resources/annotations
```

Set output paths with `--help`. `02_relabel_stored_umaps.py` maps IDs in stored
whole-cohort objects without recalculating coordinates.
