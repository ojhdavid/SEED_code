# Signature and meta-score

Scripts: `01_derive_signature.py` (shared derivation module),
`02_plot_metascore.py`, `03_raw_to_signature.py`, and `04_export_calculated_case.py`.
Use the [main workflow](../README.md) for SMC-103 derivation and execution order.
Wilcoxon-ranked genes are filtered at top 260 using the union of the five
[exclusion categories](resources/signature_exclusions.tsv), as specified in
[signature_config.json](resources/signature_config.json).

To score an existing case object in the same environment:

```bash
/new/env/visium/bin/python -s visium/06_downstream/02_plot_metascore.py \
  --input /outputs/run/case/case_for_plots_and_scoring.h5ad \
  --signature /outputs/run/analysis/derived_signature.tsv \
  --score-mode integrated-x --output-dir /outputs/run/scoring
```

`integrated-x` scores the stored expression matrix. `lognorm-unscaled` instead
normalizes raw counts to 10,000 and applies `log1p`; it requires `--counts-layer`
or `--use-raw-counts`. Use `--use-paper-signature` instead of `--signature` to
score the bundled [93-gene signature](resources/luad_smc_103_seed_meta_signature.tsv).
