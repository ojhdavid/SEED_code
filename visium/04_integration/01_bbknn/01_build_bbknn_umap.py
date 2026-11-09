#!/usr/bin/env python
"""Build whole-cohort Visium UMAP embeddings.

Purpose
-------
Create two cohort-wide embeddings from a single common expression space:
1. uncorrected PCA-neighbors-UMAP
2. BBKNN-adjusted PCA-UMAP

These embeddings are intended for cohort-wide visualization of sample, tissue,
case, and annotation structure. They are not used to define SEED.

Inputs
------
- TMA coordinate table
- Space Ranger Visium matrices
- manual annotation CSVs

Outputs
-------
- results/intermediate/joint_cohort_umap_uncorrected.h5ad
- results/intermediate/joint_cohort_umap_bbknn.h5ad
- results/figures/joint_cohort_umaps/*.png
- results/figures/joint_cohort_umaps/joint_workflow_output_manifest.tsv
- results/figures/joint_cohort_umaps/run_manifest_joint_cohort_umaps.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import re
import shutil
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


SEED = 42
STRICT_MALIGNANT_LABELS = [
    "Malignant_cell",
    "Malignant_cell_border_mixed",
    "Malignant_cell_mucin",
    "Malignant_cell(Hot)",
    "Malignant_cell(Cold)",
    "Malignant_mucin",
]
MIXED_TUMOR_ASSOCIATED_LABELS = [
    "Malignant_cell_mixed_with_stromal",
    "Malignant_cell_mixed_with_immune_cell",
]
TUMOR_ANNOTATION_GROUP_ORDER = [
    "Other non-malignant",
    "Mixed tumor-associated",
    "Strict malignant",
]
TUMOR_GROUP_PALETTE = {
    "Other non-malignant": "#D9D9D9",
    "Mixed tumor-associated": "#E67E22",
    "Strict malignant": "#A71930",
}
STRICT_BINARY_PALETTE = {
    "Other non-malignant": "#D9D9D9",
    "Strict malignant": "#A71930",
}
INCLUSIVE_BINARY_PALETTE = {
    "Other non-tumor": "#D9D9D9",
    "Tumor-associated": "#A71930",
}


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str | None:
    path = Path(path)
    if not path.exists() or not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_barcode(sample: str, barcode: str) -> str:
    barcode = str(barcode)
    return barcode if barcode.startswith(f"{sample}_") else f"{sample}_{barcode}"


def array_sha256(value: np.ndarray) -> str:
    array = np.ascontiguousarray(value)
    digest = hashlib.sha256()
    digest.update(str(array.shape).encode())
    digest.update(str(array.dtype).encode())
    digest.update(array.tobytes())
    return digest.hexdigest()


def apply_public_ids(adata: AnnData, mapping_path: str | Path) -> AnnData:
    """Map sample and case metadata to public IDs."""
    mapping = pd.read_csv(mapping_path, sep="\t", dtype=str, keep_default_na=False)
    required = {"source_case_id", "source_sample_id", "public_case_id", "public_sample_id"}
    missing = sorted(required.difference(mapping.columns))
    if missing:
        raise ValueError(f"Sample mapping is missing columns: {missing}")
    if mapping["source_sample_id"].duplicated().any() or mapping["public_sample_id"].duplicated().any():
        raise ValueError("Sample mapping must be one-to-one")
    sample_map = dict(zip(mapping["source_sample_id"], mapping["public_sample_id"]))
    case_map = dict(zip(mapping["source_case_id"], mapping["public_case_id"]))
    source_samples = adata.obs["sample"].astype(str)
    source_cases = adata.obs["case"].astype(str)
    missing_samples = sorted(set(source_samples).difference(sample_map))
    missing_cases = sorted(set(source_cases).difference(case_map))
    if missing_samples or missing_cases:
        raise ValueError(f"Unmapped source IDs: samples={missing_samples}, cases={missing_cases}")
    before = {key: array_sha256(np.asarray(adata.obsm[key])) for key in ["X_pca", "X_umap"]}
    adata.obs["source_obs_name"] = adata.obs["raw_barcode"].astype(str)
    adata.obs["source_sample_id"] = source_samples.to_numpy()
    adata.obs["source_case_id"] = source_cases.to_numpy()
    adata.obs["public_sample_id"] = source_samples.map(sample_map).to_numpy()
    adata.obs["public_case_id"] = source_cases.map(case_map).to_numpy()
    adata.obs["sample"] = adata.obs["public_sample_id"].to_numpy()
    adata.obs["case"] = adata.obs["public_case_id"].to_numpy()
    adata.obs["canonical_spot_id"] = (
        adata.obs["public_sample_id"].astype(str) + "::" + adata.obs["source_obs_name"].astype(str)
    )
    adata.obs_names = adata.obs["canonical_spot_id"].to_numpy()
    after = {key: array_sha256(np.asarray(adata.obsm[key])) for key in ["X_pca", "X_umap"]}
    if before != after:
        raise RuntimeError("Public-ID relabeling changed stored PCA or UMAP coordinates")
    adata.uns["public_id_mapping_summary"] = {
        "mapping_path": str(mapping_path),
        "coordinate_hashes_before": before,
        "coordinate_hashes_after": after,
        "coordinates_unchanged": True,
    }
    return adata


def assert_can_write(path: Path, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing output without --overwrite: {path}")


def assign_tumor_annotation_group(annotation: str) -> str:
    if annotation in STRICT_MALIGNANT_LABELS:
        return "Strict malignant"
    if annotation in MIXED_TUMOR_ASSOCIATED_LABELS:
        return "Mixed tumor-associated"
    return "Other non-malignant"


def add_tumor_annotation_columns(adata: AnnData) -> AnnData:
    annotations = adata.obs["annotation"].astype(str)
    group = annotations.map(assign_tumor_annotation_group)
    adata.obs["tumor_annotation_group"] = pd.Categorical(
        group,
        categories=TUMOR_ANNOTATION_GROUP_ORDER,
        ordered=True,
    )
    strict = annotations.isin(STRICT_MALIGNANT_LABELS)
    mixed = annotations.isin(MIXED_TUMOR_ASSOCIATED_LABELS)
    adata.obs["strict_malignant_status_excluding_mixed"] = pd.Categorical(
        np.where(strict, "Strict malignant", np.where(mixed, "Mixed tumor-associated excluded", "Other non-malignant")),
        categories=["Other non-malignant", "Strict malignant", "Mixed tumor-associated excluded"],
        ordered=True,
    )
    adata.obs["inclusive_tumor_associated_status"] = pd.Categorical(
        np.where(strict | mixed, "Tumor-associated", "Other non-tumor"),
        categories=["Other non-tumor", "Tumor-associated"],
        ordered=True,
    )
    return adata


def crop(adata: AnnData, min_row, max_row, min_col, max_col) -> AnnData:
    keep = adata.obs[
        (adata.obsm["spatial"][:, 1] >= min_row)
        & (adata.obsm["spatial"][:, 1] <= max_row)
        & (adata.obsm["spatial"][:, 0] >= min_col)
        & (adata.obsm["spatial"][:, 0] <= max_col)
    ].index
    return adata[keep, :].copy()


def apply_crop_if_needed(adata: AnnData, row: pd.Series) -> AnnData:
    values = [row.get("min_row", np.nan), row.get("max_row", np.nan), row.get("min_col", np.nan), row.get("max_col", np.nan)]
    if all(pd.isna(v) for v in values):
        return adata
    if any(pd.isna(v) for v in values):
        raise ValueError(f"Incomplete crop coordinates for {row.get('sample_sub')}: {values}")
    return crop(adata, *values)


def annotation_path(sample_sub: str, root: str | Path) -> Path:
    key = sample_sub[:-1] if sample_sub.endswith("P") else sample_sub
    path = Path(root) / f"Annotation_professor_{key}.csv"
    if not path.exists():
        raise FileNotFoundError(f"Manual annotation not found for {sample_sub}: {path}")
    return path


def build_sample_rows(coordinate_table: str | Path, mapping_path: str | Path) -> pd.DataFrame:
    raw = pd.read_csv(coordinate_table)
    mapping = pd.read_csv(mapping_path, sep="\t", dtype=str, keep_default_na=False)
    required = {"source_case_id", "source_sample_id", "cohort", "tissue", "include_in_concat"}
    missing = sorted(required.difference(mapping.columns))
    if missing:
        raise ValueError(f"Sample mapping is missing columns: {missing}")
    if mapping["source_sample_id"].duplicated().any():
        raise ValueError("Sample map is not one-to-one")
    included = mapping["include_in_concat"].str.lower().map({"true": True, "false": False})
    if included.isna().any():
        raise ValueError("include_in_concat must contain only true/false")
    primary_mapping = mapping.loc[
        mapping["cohort"].eq("SMC") & mapping["tissue"].eq("primary")
    ].copy()
    primary = pd.DataFrame(
        {
            "sample": primary_mapping["source_case_id"],
            "ALK": primary_mapping["source_case_id"],
            "sample_sub": primary_mapping["source_sample_id"],
        }
    )
    rows = pd.merge(raw.T, primary.T, left_index=True, right_index=True, how="outer").T
    rows = rows.reset_index(drop=True)
    rows = rows.loc[rows["ALK"].notna()].copy()
    unknown = sorted(set(rows["sample_sub"].astype(str)) - set(mapping["source_sample_id"]))
    if unknown:
        raise ValueError("Samples absent from mapping: " + ", ".join(unknown))
    rows = rows.loc[rows["sample_sub"].isin(mapping.loc[included, "source_sample_id"])].copy()
    rows = rows.loc[~rows["sample_sub"].astype(str).str.contains("T2", na=False)].copy()
    return rows.sort_values(["ALK", "sample_sub"]).reset_index(drop=True)


def qc_filter_counts(adata: AnnData) -> AnnData:
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    q1 = adata.obs.total_counts.quantile(0.25)
    q3 = adata.obs.total_counts.quantile(0.75)
    iqr = q3 - q1
    adata = adata[adata.obs.total_counts <= (q3 + 1.5 * iqr), :].copy()
    return adata[adata.obs.pct_counts_mt < 5, :].copy()


def load_sample(row: pd.Series, args: argparse.Namespace) -> AnnData:
    sample = str(row["sample"])
    sample_sub = str(row["sample_sub"])
    print(f"[load] {sample_sub}")
    adata = sc.read_visium(f"{args.visium_root}/{sample}/Spatial_matrix/")
    adata = apply_crop_if_needed(adata, row)
    adata = qc_filter_counts(adata)
    annotation = pd.read_csv(annotation_path(sample_sub, args.manual_annotation_root), index_col=0)
    annotation = annotation.dropna()
    annotation.columns = ["annotation"]
    annotation = annotation[~annotation.index.duplicated(keep="first")]
    adata.obs = adata.obs.merge(annotation, left_index=True, right_index=True, how="left")
    adata.obs["sample"] = sample_sub
    adata.obs["case"] = str(row["ALK"])
    adata.obs["source_sample"] = sample
    adata.obs["raw_barcode"] = adata.obs_names.astype(str)
    adata.obs["canonical_barcode"] = [normalize_barcode(sample_sub, bc) for bc in adata.obs["raw_barcode"].astype(str)]
    adata.obs["tissue"] = np.where(sample_sub.endswith("P"), "primary", "lymph_node")
    adata = adata[adata.obs["annotation"].notna() & (adata.obs["annotation"].astype(str) != "Exclude")].copy()
    adata.obs_names = adata.obs["canonical_barcode"].astype(str).to_numpy()
    return adata


def deterministic_gene_intersection(sample_adata: dict[str, AnnData]) -> dict[str, AnnData]:
    keys = list(sample_adata.keys())
    common = set(sample_adata[keys[0]].var_names)
    for key in keys[1:]:
        common &= set(sample_adata[key].var_names)
    ordered = [gene for gene in sample_adata[keys[0]].var_names if gene in common]
    return {key: adata[:, ordered].copy() for key, adata in sample_adata.items()}


def build_joint_object(args: argparse.Namespace, rows: pd.DataFrame) -> AnnData:
    sample_adata = {str(row["sample_sub"]): load_sample(row, args) for _, row in rows.iterrows()}
    sample_adata = deterministic_gene_intersection(sample_adata)
    adata = concat(list(sample_adata.values()), join="inner", label="sample", keys=list(sample_adata.keys()), index_unique=None)
    adata.layers["counts"] = adata.X.copy()
    adata = add_tumor_annotation_columns(adata)
    return adata


def cohort_preprocess(adata: AnnData, args: argparse.Namespace) -> AnnData:
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_top_genes, flavor="seurat")
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=args.n_pcs, svd_solver="arpack", random_state=SEED)
    adata.uns["cohort_umap_parameters"] = {
        "n_top_genes": args.n_top_genes,
        "n_pcs": args.n_pcs,
        "n_neighbors": args.n_neighbors,
        "random_seed": SEED,
        "deterministic_gene_order": True,
        "sample_level_pca_concatenated": False,
    }
    return adata


def build_uncorrected(adata: AnnData, args: argparse.Namespace) -> AnnData:
    out = adata.copy()
    sc.pp.neighbors(out, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs, use_rep="X_pca", random_state=SEED)
    sc.tl.umap(out, random_state=SEED)
    out.uns["embedding_method"] = {"name": "uncorrected whole-cohort PCA-neighbors-UMAP", "batch_key": None}
    return out


def build_bbknn(adata: AnnData) -> AnnData:
    out = adata.copy()
    sc.external.pp.bbknn(out, batch_key="sample", use_rep="X_pca")
    sc.tl.umap(out, random_state=SEED)
    out.uns["embedding_method"] = {"name": "BBKNN-adjusted whole-cohort PCA-UMAP", "batch_key": "sample"}
    return out


def palette(values: pd.Series, key: str) -> dict[str, str]:
    if key == "tumor_annotation_group":
        return TUMOR_GROUP_PALETTE
    if key == "strict_malignant_status_excluding_mixed":
        return {**STRICT_BINARY_PALETTE, "Mixed tumor-associated excluded": "#FFFFFF"}
    if key == "inclusive_tumor_associated_status":
        return INCLUSIVE_BINARY_PALETTE
    cats = sorted(pd.unique(values.astype(str)))
    return dict(zip(cats, sns.color_palette("tab20", n_colors=max(1, len(cats))).as_hex()))


def plot_umap(adata: AnnData, color_key: str, path: Path, title: str, point_size: float = 5) -> None:
    xy = np.asarray(adata.obsm["X_umap"])[:, :2]
    values = adata.obs[color_key].astype(str)
    pal = palette(values, color_key)
    fig, ax = plt.subplots(figsize=(10, 6))
    if color_key == "tumor_annotation_group":
        categories = TUMOR_ANNOTATION_GROUP_ORDER
    elif color_key == "strict_malignant_status_excluding_mixed":
        categories = ["Other non-malignant", "Strict malignant"]
        values = values.replace("Mixed tumor-associated excluded", np.nan)
    elif color_key == "inclusive_tumor_associated_status":
        categories = ["Other non-tumor", "Tumor-associated"]
    else:
        categories = list(pd.unique(values))
    for category in categories:
        mask = values.eq(category).to_numpy()
        if mask.sum() == 0:
            continue
        alpha = 0.40 if category in {"Other non-malignant", "Other non-tumor"} else 0.90
        size = point_size if category in {"Other non-malignant", "Other non-tumor"} else point_size * 1.4
        ax.scatter(xy[mask, 0], xy[mask, 1], s=point_size, color=pal.get(category, "gray"), alpha=0.85, linewidths=0, rasterized=True, label=f"{category} (n={mask.sum():,})")
        ax.collections[-1].set_alpha(alpha)
        ax.collections[-1].set_sizes(np.array([size]))
    ax.set_title(title)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_aspect("equal", adjustable="datalim")
    if values.nunique() <= 30:
        ax.legend(
            loc="center left", bbox_to_anchor=(0.72, 0.5), bbox_transform=fig.transFigure,
            frameon=False, fontsize=7,
        )
    fig.subplots_adjust(left=0.09, right=0.70, bottom=0.11, top=0.91)
    fig.savefig(path, dpi=300)
    plt.close(fig)


def make_figures(embeddings: dict[str, AnnData], figure_dir: Path) -> pd.DataFrame:
    figure_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    # The retained uncorrected views are deliberately limited to case and
    # sample. Selected BBKNN views use the same stored embedding.
    requested = {
        "uncorrected": ["case", "sample"],
        "bbknn": ["case", "sample", "tumor_annotation_group"],
    }
    for name, keys in requested.items():
        adata = embeddings[name]
        for key in keys:
            path = figure_dir / f"{name}_umap_by_{key}.png"
            title = f"{'BBKNN' if name == 'bbknn' else 'Uncorrected'} whole-cohort UMAP by {key}"
            plot_umap(adata, key, path, title)
            rows.append({"embedding": name, "figure_type": f"umap_by_{key}", "path": str(path)})
    return pd.DataFrame(rows)


def summarize_group_counts(adata: AnnData, group_cols: list[str], scope_values: tuple[str, ...]) -> pd.DataFrame:
    rows = []
    df = adata.obs.copy()
    if group_cols:
        grouped = df.groupby(group_cols, sort=False, observed=True)
        for keys, sub in grouped:
            if not isinstance(keys, tuple):
                keys = (keys,)
            denom = len(sub)
            counts = sub["tumor_annotation_group"].value_counts().reindex(TUMOR_ANNOTATION_GROUP_ORDER, fill_value=0)
            for group in TUMOR_ANNOTATION_GROUP_ORDER:
                rows.append(
                    {
                        **{col: str(val) for col, val in zip(group_cols, keys)},
                        "scope": scope_values[0],
                        "tumor_annotation_group": group,
                        "spot_count": int(counts[group]),
                        "denominator": int(denom),
                        "spot_fraction": float(counts[group] / denom) if denom else np.nan,
                        "spot_percentage": float(100 * counts[group] / denom) if denom else np.nan,
                    }
                )
    else:
        denom = len(df)
        counts = df["tumor_annotation_group"].value_counts().reindex(TUMOR_ANNOTATION_GROUP_ORDER, fill_value=0)
        for group in TUMOR_ANNOTATION_GROUP_ORDER:
            rows.append(
                {
                    "scope": scope_values[0],
                    "sample": "",
                    "case": "",
                    "tissue": "",
                    "tumor_annotation_group": group,
                    "spot_count": int(counts[group]),
                    "denominator": int(denom),
                    "spot_fraction": float(counts[group] / denom) if denom else np.nan,
                    "spot_percentage": float(100 * counts[group] / denom) if denom else np.nan,
                }
            )
    return pd.DataFrame(rows)


def write_count_tables(adata: AnnData, figure_dir: Path) -> dict[str, str]:
    sample = summarize_group_counts(adata, ["sample", "case", "tissue"], ("sample",))
    case = summarize_group_counts(adata, ["case"], ("case",))
    pooled_frames = []
    for scope, sub in [
        ("whole_cohort_pooled", adata),
        ("primary_samples_pooled", adata[adata.obs["tissue"].astype(str).eq("primary")].copy()),
        ("lymph_node_samples_pooled", adata[adata.obs["tissue"].astype(str).eq("lymph_node")].copy()),
    ]:
        pooled_frames.append(summarize_group_counts(sub, [], (scope,)))
    pooled = pd.concat(pooled_frames, ignore_index=True)
    paths = {
        "sample": str(figure_dir / "tumor_annotation_group_counts_by_sample.tsv"),
        "case": str(figure_dir / "tumor_annotation_group_counts_by_case.tsv"),
        "pooled": str(figure_dir / "tumor_annotation_group_counts_pooled.tsv"),
    }
    required_cols = ["scope", "sample", "case", "tissue", "tumor_annotation_group", "spot_count", "denominator", "spot_fraction", "spot_percentage"]
    for frame in [sample, case, pooled]:
        for col in required_cols:
            if col not in frame.columns:
                frame[col] = ""
    sample = sample[required_cols]
    case = case[required_cols]
    pooled = pooled[required_cols]
    sample.to_csv(paths["sample"], sep="\t", index=False)
    case.to_csv(paths["case"], sep="\t", index=False)
    pooled.to_csv(paths["pooled"], sep="\t", index=False)
    return paths


def annotation_reason(annotation: str, group: str) -> str:
    if annotation == "Malignant_cell_border_mixed":
        return "Strict malignant: this label denotes malignant spots in matched lymph-node samples."
    if annotation in STRICT_MALIGNANT_LABELS:
        return "Strict malignant annotation label."
    if annotation in MIXED_TUMOR_ASSOCIATED_LABELS:
        return "Mixed tumor-associated label; separated from strict malignant and not assigned to non-malignant."
    return "Other non-malignant: not an explicit malignant or mixed tumor-associated label."


def write_annotation_mapping_table(adata: AnnData, figure_dir: Path) -> str:
    rows = []
    obs = adata.obs.copy()
    for ann, sub in obs.groupby("annotation", sort=True, observed=True):
        group = assign_tumor_annotation_group(str(ann))
        rows.append(
            {
                "original_annotation": str(ann),
                "total_count": int(len(sub)),
                "primary_count": int(sub["tissue"].astype(str).eq("primary").sum()),
                "lymph_node_count": int(sub["tissue"].astype(str).eq("lymph_node").sum()),
                "samples_containing_label": ",".join(sorted(sub["sample"].astype(str).unique())),
                "assigned_tumor_annotation_group": group,
                "reason_for_assignment": annotation_reason(str(ann), group),
            }
        )
    out = figure_dir / "tumor_annotation_label_mapping.tsv"
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    return str(out)


def image_nonblank(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        arr = plt.imread(path)
        return bool(np.nanstd(arr) > 0)
    except Exception:
        return False


def write_input_checks(uncorrected: AnnData, bbknn: AnnData, figure_manifest: pd.DataFrame, archive_manifest: pd.DataFrame, figure_dir: Path, args: argparse.Namespace) -> str:
    obs = uncorrected.obs
    labels_present = set(obs["annotation"].astype(str))
    absent = sorted((set(STRICT_MALIGNANT_LABELS) | set(MIXED_TUMOR_ASSOCIATED_LABELS)) - labels_present)
    border = obs[obs["annotation"].astype(str).eq("Malignant_cell_border_mixed")]
    mixed = obs["annotation"].astype(str).isin(MIXED_TUMOR_ASSOCIATED_LABELS)
    checks = []
    def add(check: str, status: str, details: str) -> None:
        checks.append({"check": check, "status": status, "details": details})

    add("Malignant_cell_border_mixed assigned to Strict malignant", "PASS" if assign_tumor_annotation_group("Malignant_cell_border_mixed") == "Strict malignant" else "FAIL", "Required lymph-node malignant label classification.")
    add("Malignant_cell_mixed_with_stromal assigned to Mixed tumor-associated", "PASS" if assign_tumor_annotation_group("Malignant_cell_mixed_with_stromal") == "Mixed tumor-associated" else "FAIL", "Required mixed tumor-associated label classification.")
    add("Malignant_cell_mixed_with_immune_cell assigned to Mixed tumor-associated", "PASS" if assign_tumor_annotation_group("Malignant_cell_mixed_with_immune_cell") == "Mixed tumor-associated" else "FAIL", "Required mixed tumor-associated label classification.")
    add("Mixed tumor-associated spots not assigned to non-malignant", "PASS" if (obs.loc[mixed, "tumor_annotation_group"].astype(str) == "Mixed tumor-associated").all() else "FAIL", f"Mixed tumor-associated spot count={int(mixed.sum())}.")
    add("Three categories are mutually exclusive", "PASS" if obs["tumor_annotation_group"].notna().all() else "FAIL", "Stored as one ordered categorical value per retained spot.")
    add("Every retained annotated spot assigned exactly one group", "PASS" if obs["tumor_annotation_group"].notna().sum() == uncorrected.n_obs else "FAIL", f"Assigned={int(obs['tumor_annotation_group'].notna().sum())}; total={uncorrected.n_obs}.")
    counts_sum = int(obs["tumor_annotation_group"].value_counts().sum())
    add("Category counts sum to retained spots", "PASS" if counts_sum == uncorrected.n_obs else "FAIL", f"Category sum={counts_sum}; total={uncorrected.n_obs}.")
    add("Malignant_cell_border_mixed distribution recorded", "PASS", f"Total={len(border)}; primary={int(border['tissue'].astype(str).eq('primary').sum())}; lymph_node={int(border['tissue'].astype(str).eq('lymph_node').sum())}; samples={','.join(sorted(border['sample'].astype(str).unique()))}.")
    add("Uncorrected object contains tumor_annotation_group", "PASS" if "tumor_annotation_group" in uncorrected.obs.columns else "FAIL", "Categorical column in uncorrected h5ad.")
    add("BBKNN object contains tumor_annotation_group", "PASS" if "tumor_annotation_group" in bbknn.obs.columns else "FAIL", "Categorical column in BBKNN h5ad.")
    add("UMAP coordinates present", "PASS" if "X_umap" in uncorrected.obsm and "X_umap" in bbknn.obsm else "FAIL", "Uncorrected and BBKNN objects contain X_umap.")
    generated_figs = [Path(p) for p in figure_manifest["path"] if str(p).endswith(".png")]
    nonblank = [image_nonblank(p) for p in generated_figs]
    add("All generated PNG figures are nonblank", "PASS" if all(nonblank) else "FAIL", f"Nonblank={sum(nonblank)}/{len(nonblank)}.")
    add("Absent requested labels recorded", "PASS", "Absent labels: " + (", ".join(absent) if absent else "none"))
    out = figure_dir / "input_checks.tsv"
    pd.DataFrame(checks).to_csv(out, sep="\t", index=False)
    return str(out)


def archive_previous_outputs(args: argparse.Namespace) -> pd.DataFrame:
    archive_dir = Path(args.archive_dir)
    archived = []
    if archive_dir.exists():
        raise FileExistsError(f"Archive directory already exists; refusing to overwrite: {archive_dir}")
    archive_dir.mkdir(parents=True, exist_ok=False)
    for path in [Path(args.output_uncorrected), Path(args.output_bbknn)]:
        if path.exists():
            dest = archive_dir / path.name
            shutil.copy2(path, dest)
            archived.append({"source_path": str(path), "archive_path": str(dest), "type": "h5ad"})
    figure_dir = Path(args.figure_dir)
    if figure_dir.exists():
        fig_archive = archive_dir / "joint_cohort_umaps"
        fig_archive.mkdir(parents=True, exist_ok=False)
        for path in sorted(figure_dir.iterdir()):
            if path.is_file():
                dest = fig_archive / path.name
                shutil.copy2(path, dest)
                archived.append({"source_path": str(path), "archive_path": str(dest), "type": "figure_or_manifest"})
    manifest = pd.DataFrame(archived)
    manifest.to_csv(archive_dir / "archive_manifest.tsv", sep="\t", index=False)
    return manifest


def write_workflow_summary(adata: AnnData, figure_dir: Path, archive_manifest: pd.DataFrame) -> str:
    counts = adata.obs["tumor_annotation_group"].value_counts().reindex(TUMOR_ANNOTATION_GROUP_ORDER, fill_value=0)
    total = int(adata.n_obs)
    text = [
        "# Whole-cohort UMAP summary", "",
        "Pathology annotations are grouped as Strict malignant, Mixed tumor-associated, and Other non-malignant.",
        "`Malignant_cell_border_mixed` belongs to Strict malignant. Tumor-stromal and tumor-immune mixed labels belong to Mixed tumor-associated.",
        "", "The workflow writes uncorrected and BBKNN embeddings, annotation plots and group-count tables.", "",
    ]
    for group in TUMOR_ANNOTATION_GROUP_ORDER:
        text.append(f"- {group}: {int(counts[group]):,} / {total:,} ({100 * counts[group] / total:.2f}%)")
    out = figure_dir / "workflow_summary.md"
    out.write_text("\n".join(text) + "\n", encoding="utf-8")
    return str(out)


def write_manifest(args: argparse.Namespace, rows: pd.DataFrame, adata: AnnData, figure_manifest: pd.DataFrame, archive_manifest: pd.DataFrame, count_paths: dict[str, str], annotation_mapping_path: str, checks_path: str, workflow_summary_path: str) -> None:
    labels_present = set(adata.obs["annotation"].astype(str))
    manifest = {
        "script": "visium/04_integration/01_bbknn/01_build_bbknn_umap.py",
        "random_seed": SEED,
        "parameters": {"n_top_genes": args.n_top_genes, "n_pcs": args.n_pcs, "n_neighbors": args.n_neighbors},
        "annotation_definitions": {
            "strict_malignant_labels": STRICT_MALIGNANT_LABELS,
            "mixed_tumor_associated_labels": MIXED_TUMOR_ASSOCIATED_LABELS,
            "tumor_annotation_group_order": TUMOR_ANNOTATION_GROUP_ORDER,
            "tumor_group_palette": TUMOR_GROUP_PALETTE,
            "absent_requested_labels": sorted((set(STRICT_MALIGNANT_LABELS) | set(MIXED_TUMOR_ASSOCIATED_LABELS)) - labels_present),
        },
        "inputs": {
            "coordinate_table": {"path": args.coordinate_table, "sha256": sha256_file(args.coordinate_table)},
            "visium_root": args.visium_root,
            "manual_annotation_root": args.manual_annotation_root,
        },
        "sample_count": int(rows.shape[0]),
        "spot_count": int(adata.n_obs),
        "hvg_count": int(adata.n_vars),
        "outputs": {
            "uncorrected_h5ad": args.output_uncorrected,
            "bbknn_h5ad": args.output_bbknn,
            "figure_dir": args.figure_dir,
            "figure_count": int(len(figure_manifest)),
            "count_tables": count_paths,
            "annotation_mapping": annotation_mapping_path,
            "input_checks": checks_path,
            "workflow_summary": workflow_summary_path,
        },
        "archive": {
            "archive_dir": args.archive_dir,
            "archived_file_count": int(len(archive_manifest)),
        },
        "software": {
            "python": platform.python_version(),
            "scanpy": sc.__version__,
            "scipy": scipy.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
    }
    Path(args.figure_dir).mkdir(parents=True, exist_ok=True)
    (Path(args.figure_dir) / "run_manifest_joint_cohort_umaps.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build uncorrected and BBKNN whole-cohort Visium UMAPs.")
    parser.add_argument("--coordinate-table", default="visium/resources/coordinates/TMA_coordinates_except_german.csv")
    parser.add_argument("--sample-map", default="visium/resources/sample_id_mapping.tsv")
    parser.add_argument("--visium-root", default="data/visium/Visium")
    parser.add_argument("--manual-annotation-root", default="visium/resources/annotations")
    parser.add_argument("--output-uncorrected", default="results/intermediate/joint_cohort_umap_uncorrected.h5ad")
    parser.add_argument("--output-bbknn", default="results/intermediate/joint_cohort_umap_bbknn.h5ad")
    parser.add_argument("--figure-dir", default="results/figures/joint_cohort_umaps")
    parser.add_argument("--archive-dir", default="results/qc/joint_cohort_umap_output_backup")
    parser.add_argument("--n-top-genes", type=int, default=3000)
    parser.add_argument("--n-pcs", type=int, default=50)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    np.random.seed(SEED)
    sc.settings.verbosity = 2
    archive_manifest = archive_previous_outputs(args) if args.overwrite else pd.DataFrame()
    for path in [Path(args.output_uncorrected), Path(args.output_bbknn)]:
        path.parent.mkdir(parents=True, exist_ok=True)
        assert_can_write(path, args.overwrite)
    figure_dir = Path(args.figure_dir)
    figure_dir.mkdir(parents=True, exist_ok=True)

    rows = build_sample_rows(args.coordinate_table, args.sample_map)
    rows.to_csv(figure_dir / "joint_cohort_input_samples.tsv", sep="\t", index=False)
    adata = cohort_preprocess(build_joint_object(args, rows), args)
    uncorrected = build_uncorrected(adata, args)
    bbknn = build_bbknn(adata)
    uncorrected = apply_public_ids(uncorrected, args.sample_map)
    bbknn = apply_public_ids(bbknn, args.sample_map)

    uncorrected.write_h5ad(args.output_uncorrected)
    bbknn.write_h5ad(args.output_bbknn)
    figure_manifest = make_figures({"uncorrected": uncorrected, "bbknn": bbknn}, figure_dir)
    figure_manifest.to_csv(figure_dir / "joint_workflow_output_manifest.tsv", sep="\t", index=False)
    count_paths = write_count_tables(uncorrected, figure_dir)
    annotation_mapping_path = write_annotation_mapping_table(uncorrected, figure_dir)
    checks_path = write_input_checks(uncorrected, bbknn, figure_manifest, archive_manifest, figure_dir, args)
    workflow_summary_path = write_workflow_summary(uncorrected, figure_dir, archive_manifest)
    write_manifest(args, rows, uncorrected, figure_manifest, archive_manifest, count_paths, annotation_mapping_path, checks_path, workflow_summary_path)
    print(f"[done] spots={uncorrected.n_obs:,} HVGs={uncorrected.n_vars:,} figures={len(figure_manifest)}")


if __name__ == "__main__":
    main()
