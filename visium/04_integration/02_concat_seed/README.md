# Matched-case UMAPs

Use the [integration environment](../README.md) and run from the repository root:

```bash
python visium/04_integration/02_concat_seed/01_build_concat_case_umaps.py \
  --visium-root /inputs/Visium --manual-annotation-root visium/resources/annotations \
  --coordinate-tables visium/resources/coordinates/TMA_coordinates_except_german.csv \
    visium/resources/coordinates/TMA_coordinates_germany.csv \
  --cases LUAD-SMC-103 LUAD-UHC-04 --output-dir /outputs/case_objects \
  --figure-dir /outputs/case_umaps --metadata-tsv /outputs/case_umaps.tsv
```

Public case IDs accept spaces or commas. The script writes per-case H5AD objects
and sample-colored UMAPs for [SEED selection](../../05_seed_definition/README.md).

An optional concat-specific environment is provided in [environment-concat.yml](environment-concat.yml).
