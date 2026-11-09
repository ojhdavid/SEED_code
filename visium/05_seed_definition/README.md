# SEED definition

Primary SEED spots are malignant primary-tumor spots whose minimum Euclidean
distance to qualifying LN anchors in the matched-case UMAP is **strictly below**
the [case threshold](../resources/case_specific_seed_thresholds.tsv).
Anchors use `Malignant_cell_border_mixed`; candidates use `Malignant_cell`, with
Hot/Cold labels combined. **UHC-04 = 1.0**; SMC-51 has no qualifying LN anchor.

`is_primary_seed` identifies primary SEED, `is_ln_anchor` identifies LN anchors,
and `is_seed` denotes the combined analysis group, including nearby LN malignant
spots. The selector filters anchors by label; `--anchor-tissue lymph_node`
adds an explicit tissue restriction when required by the input annotations.

```bash
python visium/05_seed_definition/01_define_seed.py \
  --input /inputs/case_umap.h5ad --cases LUAD-UHC-04 \
  --case-key public_case_id --sample-key public_sample_id \
  --output-h5ad /outputs/UHC04_seed.h5ad --output-csv /outputs/UHC04_assignments.csv \
  --output-summary /outputs/UHC04_summary.json
```

Use the [integration environment](../04_integration/README.md). Input requires
unique sample-qualified IDs, `X_umap`, annotations, public case/sample IDs and
`tissue` (`primary`/`lymph_node`). Plotting and external early-label attachment
scripts provide `--help`; early display labels leave SEED membership intact.
