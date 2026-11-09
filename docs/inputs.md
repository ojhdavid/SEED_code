# Visium inputs

The SMC-103 runner uses primary samples ALK51/59/101/102/103 and their TMA cores.
Spot annotations, SMC/UHC crop coordinates and SMC-103 early-group barcodes are
included under [`visium/resources/`](../visium/resources/). Supply the Space Ranger
counts and spatial files. The [sample map](../visium/resources/sample_id_mapping.tsv)
links source/public IDs and selects analysis samples.

Set `visium_root` and `output` in a copy of
[visium.example.json](../config/visium.example.json); the other paths use bundled files:

| Key | Required content |
|---|---|
| `visium_root` | Run directories `ALK51`, `ALK59`, `ALK101`, `ALK102`, `ALK103`, `TMA1`, `TMA2`; each contains `Spatial_matrix/filtered_feature_bc_matrix.h5` and `Spatial_matrix/spatial/` positions, images and scalefactors |
| `coordinate_table` | Bundled SMC crop CSV: `ALK`, `sample`, `sample_sub`, `min_row`, `max_row`, `min_col`, `max_col` |
| `annotation_root` | Bundled `Annotation_professor_<source>.csv` files: barcode index and one annotation column |
| `trajectory_csv` | Bundled SMC-103 early-group barcode list, column `x` |
| `output` | A new writable directory; relative paths resolve against the configuration file |

Counts, annotations and spatial metadata must use matching barcodes.
UHC crops use the bundled `coordinates/TMA_coordinates_germany.csv`.
