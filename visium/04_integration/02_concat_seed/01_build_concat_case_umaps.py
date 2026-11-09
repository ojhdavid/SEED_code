#!/usr/bin/env python
"""Build matched primary/LN case UMAPs with public identifiers.

Use 20,000 HVGs, case-level neighbors and sample-to-public ID mapping.
The final figures show the combined case UMAP colored by public sample ID.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import platform
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import seaborn as sns
from anndata import AnnData, concat
from scipy.spatial.distance import cdist


RANDOM_SEED = 42
ANCHOR_LABEL = "Malignant_cell_border_mixed"
CANDIDATE_LABEL = "Malignant_cell"
MALIGNANT_LABEL_REPLACEMENTS = {
    "Malignant_cell(Hot)": "Malignant_cell",
    "Malignant_cell(Cold)": "Malignant_cell",
}
MALIGNANT_LABELS = {
    "Malignant_cell", "Malignant_cell_border_mixed", "Malignant_cell_mucin",
    "Malignant_cell_mixed_with_immune_cell", "Malignant_cell_mixed_with_stromal",
    "Malignant_cell(Hot)", "Malignant_cell(Cold)", "Malignant_mucin",
}


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str | None:
    path = Path(path)
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def array_sha256(value: np.ndarray) -> str:
    array = np.ascontiguousarray(value)
    digest = hashlib.sha256()
    digest.update(str(array.shape).encode())
    digest.update(str(array.dtype).encode())
    digest.update(array.tobytes())
    return digest.hexdigest()


def filename_safe(value: object) -> str:
    return "".join(ch if ch.isalnum() or ch in "._" else "_" for ch in str(value))


def read_mapping(path: str | Path) -> pd.DataFrame:
    mapping = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {
        "source_case_id", "source_sample_id", "public_case_id", "public_sample_id",
        "cohort", "tissue", "core_id", "include_in_concat", "exclusion_reason",
        "space_ranger_sample", "coordinate_sample_id", "annotation_source_id",
        "coordinate_table",
    }
    missing = sorted(required.difference(mapping.columns))
    if missing:
        raise ValueError(f"Sample mapping is missing columns: {missing}")
    for key in ("source_sample_id", "public_sample_id"):
        if mapping[key].duplicated().any():
            values = mapping.loc[mapping[key].duplicated(False), key].tolist()
            raise ValueError(f"Mapping column {key} is not one-to-one: {values}")
    mapping["include_in_concat"] = mapping["include_in_concat"].str.lower().map(
        {"true": True, "false": False}
    )
    if mapping["include_in_concat"].isna().any():
        raise ValueError("include_in_concat must contain only true/false")
    return mapping


def read_coordinate_tables(paths: list[str]) -> dict[str, pd.DataFrame]:
    tables: dict[str, pd.DataFrame] = {}
    for value in paths:
        path = Path(value)
        table = pd.read_csv(path)
        required = {"sample", "sample_sub", "min_row", "max_row", "min_col", "max_col"}
        if not required.issubset(table.columns):
            raise ValueError(f"{path} lacks required coordinate columns")
        tables[path.name] = table
    return tables


def build_sample_rows(mapping: pd.DataFrame, coordinate_tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for _, record in mapping.loc[mapping["include_in_concat"]].iterrows():
        row: dict[str, object] = record.to_dict()
        if record["tissue"] == "primary":
            row.update({key: np.nan for key in ["min_row", "max_row", "min_col", "max_col"]})
        else:
            table_name = record["coordinate_table"]
            if table_name not in coordinate_tables:
                raise FileNotFoundError(f"Coordinate table not supplied: {table_name}")
            table = coordinate_tables[table_name]
            match = table.loc[table["sample_sub"].astype(str).eq(record["coordinate_sample_id"])]
            if len(match) != 1:
                raise ValueError(
                    f"Expected one coordinate row for {record['source_sample_id']}; found {len(match)}"
                )
            coordinate = match.iloc[0]
            if str(coordinate["sample"]) != record["space_ranger_sample"]:
                raise ValueError(f"Space Ranger source mismatch for {record['source_sample_id']}")
            row.update({key: coordinate[key] for key in ["min_row", "max_row", "min_col", "max_col"]})
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["public_case_id", "source_sample_id"]).reset_index(drop=True)


def crop(adata: AnnData, row: pd.Series) -> AnnData:
    values = [row[key] for key in ["min_row", "max_row", "min_col", "max_col"]]
    if all(pd.isna(value) for value in values):
        return adata
    if any(pd.isna(value) for value in values):
        raise ValueError(f"Incomplete crop coordinates for {row['source_sample_id']}: {values}")
    min_row, max_row, min_col, max_col = values
    keep = (
        (adata.obsm["spatial"][:, 1] >= min_row)
        & (adata.obsm["spatial"][:, 1] <= max_row)
        & (adata.obsm["spatial"][:, 0] >= min_col)
        & (adata.obsm["spatial"][:, 0] <= max_col)
    )
    return adata[keep, :].copy()


def preprocess_sample(adata: AnnData, n_top_genes: int) -> AnnData:
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    q1 = adata.obs.total_counts.quantile(0.25)
    q3 = adata.obs.total_counts.quantile(0.75)
    adata = adata[adata.obs.total_counts <= q3 + 1.5 * (q3 - q1), :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=n_top_genes
    )
    adata = adata[:, adata.var["highly_variable"]]
    sc.pp.scale(adata)
    sc.tl.pca(adata, svd_solver="arpack", random_state=RANDOM_SEED)
    return adata


def load_sample(row: pd.Series, args: argparse.Namespace) -> AnnData:
    source_sample = row["source_sample_id"]
    print(f"[load] {source_sample} -> {row['public_sample_id']}")
    adata = sc.read_visium(
        str(Path(args.visium_root) / row["space_ranger_sample"] / "Spatial_matrix")
    )
    input_spots_before_crop = int(adata.n_obs)
    adata = crop(adata, row)
    input_spots_after_crop = int(adata.n_obs)
    adata = preprocess_sample(adata, args.n_top_genes)
    annotation_path = Path(args.manual_annotation_root) / f"Annotation_professor_{row['annotation_source_id']}.csv"
    annotation = pd.read_csv(annotation_path, index_col=0).dropna()
    if annotation.shape[1] != 1:
        raise ValueError(f"Expected one annotation column in {annotation_path}")
    annotation.columns = ["annotation"]
    annotation = annotation.loc[~annotation.index.duplicated(keep="first")]
    adata.obs = adata.obs.merge(annotation, left_index=True, right_index=True, how="left")
    adata = adata[adata.obs["annotation"].notna() & adata.obs["annotation"].astype(str).ne("Exclude")].copy()
    adata.obs["source_obs_name"] = adata.obs_names.astype(str)
    adata.obs["source_case_id"] = row["source_case_id"]
    adata.obs["source_sample_id"] = source_sample
    adata.obs["public_case_id"] = row["public_case_id"]
    adata.obs["public_sample_id"] = row["public_sample_id"]
    adata.obs["case"] = row["public_case_id"]
    adata.obs["sample"] = row["public_sample_id"]
    adata.obs["cohort"] = row["cohort"]
    adata.obs["tissue"] = row["tissue"]
    adata.obs["core_id"] = row["core_id"]
    adata.obs["canonical_spot_id"] = (
        adata.obs["public_sample_id"].astype(str) + "::" + adata.obs["source_obs_name"].astype(str)
    )
    adata.obs["malignant_status"] = np.where(
        adata.obs["annotation"].astype(str).isin(MALIGNANT_LABELS), "Malignant", "Non-malignant"
    )
    adata.obs_names = adata.obs["canonical_spot_id"].to_numpy()
    adata.uns["sample_input_provenance"] = {
        "source_sample_id": source_sample,
        "public_sample_id": row["public_sample_id"],
        "spots_before_crop": input_spots_before_crop,
        "spots_after_crop": input_spots_after_crop,
        "spots_after_annotation_filter": int(adata.n_obs),
    }
    return adata


def intersect_genes(samples: dict[str, AnnData]) -> dict[str, AnnData]:
    keys = list(samples)
    common = set(samples[keys[0]].var_names)
    for key in keys[1:]:
        common.intersection_update(samples[key].var_names)
    ordered = [gene for gene in samples[keys[0]].var_names if gene in common]
    return {key: value[:, ordered].copy() for key, value in samples.items()}


def build_case(public_case_id: str, rows: pd.DataFrame, args: argparse.Namespace) -> tuple[AnnData, list[str]]:
    case_rows = rows.loc[rows["public_case_id"].eq(public_case_id)]
    if case_rows.empty:
        raise ValueError(f"No included sample rows for {public_case_id}")
    samples = {row["source_sample_id"]: load_sample(row, args) for _, row in case_rows.iterrows()}
    dropped = sorted(key for key, value in samples.items() if value.n_obs < args.min_cells_per_sample)
    samples = {key: value for key, value in samples.items() if value.n_obs >= args.min_cells_per_sample}
    if len(samples) < 2:
        raise ValueError(f"{public_case_id}: fewer than two samples remain after min-cells filtering")
    hvg_counts = {key: int(value.n_vars) for key, value in samples.items()}
    input_spot_counts = {
        key: dict(value.uns["sample_input_provenance"]) for key, value in samples.items()
    }
    samples = intersect_genes(samples)
    common_gene_count = next(iter(samples.values())).n_vars

    case_adata = concat(list(samples.values()), join="inner", label="_concat_key", keys=list(samples), index_unique=None)
    # Build BBKNN, then use ordinary neighbors for the case embedding.
    sc.external.pp.bbknn(case_adata, batch_key="sample")
    sc.pp.neighbors(case_adata, n_neighbors=8, random_state=RANDOM_SEED)
    sc.tl.leiden(case_adata, resolution=0.3, random_state=RANDOM_SEED)
    sc.tl.umap(case_adata, random_state=RANDOM_SEED)
    for source_sample_id, sample_adata in samples.items():
        mask = case_adata.obs["source_sample_id"].astype(str).eq(source_sample_id).to_numpy()
        sample_adata.obsm["X_umap"] = np.asarray(case_adata.obsm["X_umap"])[mask].copy()
        sample_adata.obs["leiden"] = case_adata.obs.loc[mask, "leiden"].astype(str).to_numpy()

    adata = concat(list(samples.values()), join="inner", label="_concat_key", keys=list(samples), index_unique=None)
    if not adata.obs_names.equals(case_adata.obs_names):
        raise RuntimeError(f"{public_case_id}: final concat changed spot order")
    for key in ("X_pca", "X_umap"):
        if not np.array_equal(np.asarray(adata.obsm[key]), np.asarray(case_adata.obsm[key])):
            raise RuntimeError(f"{public_case_id}: final concat changed {key}")
    for key in ("connectivities", "distances"):
        adata.obsp[key] = case_adata.obsp[key].copy()
    adata.uns["neighbors"] = copy.deepcopy(case_adata.uns["neighbors"])
    # Final concat retains the copied UMAP; cohort integration calls BBKNN again.
    sc.external.pp.bbknn(adata, batch_key="sample")
    adata.uns["matched_case_umap_parameters"] = {
        "public_case_id": public_case_id,
        "source_case_id": sorted(adata.obs["source_case_id"].astype(str).unique()),
        "n_top_genes": int(args.n_top_genes),
        "hvg_counts_before_intersection": hvg_counts,
        "sample_input_spot_counts": input_spot_counts,
        "common_gene_count": int(common_gene_count),
        "n_neighbors": 8, "leiden_resolution": 0.3, "random_seed": RANDOM_SEED,
        "bbknn_call_omitted": False,
        "case_graph_sequence": "BBKNN then ordinary neighbors overwrites graph",
        "final_concat_graph": "case-local BBKNN; no UMAP recomputation",
        "sample_level_pca_retained": True,
        "case_umap_copied_to_samples_before_final_concat": True,
        "umap_recomputed_after_final_concat": False,
        "dropped_samples_below_min_cells": dropped,
        "x_pca_sha256": array_sha256(np.asarray(adata.obsm["X_pca"])),
        "x_umap_sha256": array_sha256(np.asarray(adata.obsm["X_umap"])),
    }
    return adata, dropped


def plot_combined_sample_umap(adata: AnnData, path: Path, title: str, dpi: int) -> None:
    xy = np.asarray(adata.obsm["X_umap"])
    values = adata.obs["public_sample_id"].astype(str)
    categories = sorted(values.unique())
    colors = dict(zip(categories, sns.color_palette("tab20", len(categories)).as_hex()))
    fig, ax = plt.subplots(figsize=(7.2, 5.7))
    for category in categories:
        mask = values.eq(category).to_numpy()
        ax.scatter(xy[mask, 0], xy[mask, 1], s=9, c=colors[category], linewidths=0,
                   rasterized=True, label=f"{category} (n={mask.sum():,})")
    ax.set_title(title)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_aspect("equal", adjustable="datalim")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def seed_summary(adata: AnnData, radius: float) -> dict[str, object]:
    """Report the distance-selected analysis group; do not alter the embedding."""
    labels = adata.obs["annotation"].astype(str).replace(MALIGNANT_LABEL_REPLACEMENTS)
    anchors = adata.obs_names[labels.eq(ANCHOR_LABEL)]
    candidates = adata.obs_names[labels.eq(CANDIDATE_LABEL)]
    nearby = []
    if len(anchors) and len(candidates):
        xy = pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names)
        distances = cdist(xy.loc[candidates], xy.loc[anchors])
        nearby = candidates[distances.min(axis=1) < radius].tolist()
    combined = adata.obs_names.isin(list(anchors) + nearby)
    primary = adata.obs["tissue"].astype(str).eq("primary").to_numpy()
    ln = adata.obs["tissue"].astype(str).eq("lymph_node").to_numpy()
    return {
        "radius": radius, "operator": "<", "n_anchors": len(anchors),
        "n_ln_anchors": int((adata.obs_names.isin(anchors) & ln).sum()),
        "n_nearby_primary": int((adata.obs_names.isin(nearby) & primary).sum()),
        "n_nearby_LN": int((adata.obs_names.isin(nearby) & ln).sum()),
        "n_primary_seed": int((combined & primary).sum()),
        "n_LN_seed": int((combined & ln).sum()),
        "n_seed_total": int(combined.sum()),
        "primary_seed_spot_ids": adata.obs_names[combined & primary].tolist(),
        "membership_stage": "distance-selected membership",
    }


def resolve_cases(values: str | list[str], mapping: pd.DataFrame) -> list[str]:
    tokens = [values] if isinstance(values, str) else values
    requested = [value.strip() for token in tokens for value in token.split(",") if value.strip()]
    aliases = dict(zip(mapping["source_case_id"], mapping["public_case_id"]))
    public = set(mapping["public_case_id"])
    resolved = [value if value in public else aliases.get(value, value) for value in requested]
    unknown = sorted(set(resolved).difference(public))
    if unknown:
        raise ValueError(f"Unknown case IDs: {unknown}")
    return list(dict.fromkeys(resolved))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-map", default="visium/resources/sample_id_mapping.tsv")
    parser.add_argument("--coordinate-tables", nargs="+", default=[
        "visium/resources/coordinates/TMA_coordinates_except_german.csv",
        "visium/resources/coordinates/TMA_coordinates_germany.csv",
    ])
    parser.add_argument("--visium-root", default="data/visium/Visium")
    parser.add_argument("--manual-annotation-root", default="visium/resources/annotations")
    parser.add_argument("--output-dir", default="results/concat_seed/objects")
    parser.add_argument("--figure-dir", default="results/concat_seed/figures")
    parser.add_argument("--metadata-tsv", default="results/concat_seed/matched_case_umap_metadata.tsv")
    parser.add_argument("--cases", nargs="+", default=[
        "LUAD-SMC-51", "LUAD-SMC-59", "LUAD-SMC-101", "LUAD-SMC-102",
        "LUAD-SMC-103", "LUAD-UHC-02", "LUAD-UHC-04", "LUAD-UHC-24",
    ])
    parser.add_argument("--n-top-genes", type=int, default=20000)
    parser.add_argument("--min-cells-per-sample", type=int, default=10)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    np.random.seed(RANDOM_SEED)
    mapping = read_mapping(args.sample_map)
    coordinate_tables = read_coordinate_tables(args.coordinate_tables)
    cases = resolve_cases(args.cases, mapping)
    rows = build_sample_rows(mapping.loc[mapping["public_case_id"].isin(cases)], coordinate_tables)
    output_dir, figure_dir = Path(args.output_dir), Path(args.figure_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    Path(args.metadata_tsv).parent.mkdir(parents=True, exist_ok=True)
    metadata_rows, figure_rows, seed_summaries = [], [], {}
    for case_id in cases:
        stem = filename_safe(case_id)
        output = output_dir / f"{stem}_case_level_umap.h5ad"
        if output.exists() and not args.overwrite:
            raise FileExistsError(f"Refusing to overwrite {output}; pass --overwrite")
        adata, dropped = build_case(case_id, rows, args)
        adata.write_h5ad(output)
        figure = figure_dir / f"{stem}_combined_umap_by_sample.png"
        plot_combined_sample_umap(adata, figure, f"{case_id}: matched primary/LN samples", args.dpi)
        figure_rows.append({"public_case_id": case_id, "figure_type": "combined_umap_by_sample", "path": str(figure)})
        for sample_id, count in adata.obs["public_sample_id"].astype(str).value_counts().sort_index().items():
            source = adata.obs.loc[adata.obs["public_sample_id"].astype(str).eq(sample_id), "source_sample_id"].iloc[0]
            metadata_rows.append({"public_case_id": case_id, "public_sample_id": sample_id,
                                  "source_sample_id": source, "n_spots": int(count),
                                  "n_genes": int(adata.n_vars), "h5ad": str(output),
                                  "dropped_below_min_cells": source in dropped})
        if case_id == "LUAD-SMC-103":
            seed_summaries[case_id] = seed_summary(adata, radius=1.5)

    pd.DataFrame(metadata_rows).to_csv(args.metadata_tsv, sep="\t", index=False)
    pd.DataFrame(figure_rows).to_csv(figure_dir / "matched_case_workflow_output_manifest.tsv", sep="\t", index=False)
    manifest = {
        "script": "visium/04_integration/02_concat_seed/01_build_concat_case_umaps.py",
        "cases": cases, "random_seed": RANDOM_SEED, "n_top_genes": args.n_top_genes,
        "default_n_top_genes": 20000,
        "sample_mapping": {"path": args.sample_map, "sha256": sha256_file(args.sample_map)},
        "coordinate_tables": [{"path": path, "sha256": sha256_file(path)} for path in args.coordinate_tables],
        "umap_graph": "sc.pp.neighbors(n_neighbors=8, random_state=42) after BBKNN",
        "final_graph": "case-local BBKNN after final concat; copied UMAP retained",
        "final_figures": "combined UMAP by public sample ID only",
        "seed_summary": seed_summaries,
        "software": {"python": platform.python_version(), "scanpy": sc.__version__,
                     "scipy": scipy.__version__, "numpy": np.__version__, "pandas": pd.__version__},
    }
    (figure_dir / "run_manifest_matched_case_umaps.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"[done] cases={len(cases)} figures={len(figure_rows)}")


if __name__ == "__main__":
    main()
