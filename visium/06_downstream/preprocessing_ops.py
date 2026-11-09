"""Per-sample QC, case preparation and BBKNN integration."""
from __future__ import annotations
import hashlib
import math
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import concat

SEED = 42

def sha256_bytes(*parts: bytes) -> str:
    digest = hashlib.sha256()
    for part in parts:
        digest.update(part)
    return digest.hexdigest()

def array_digest(value) -> str:
    arr = np.ascontiguousarray(np.asarray(value))
    return sha256_bytes(str(arr.dtype).encode(), str(arr.shape).encode(), arr.tobytes())

def crop(adata, min_row, max_row, min_col, max_col):
    cropped_cells = adata.obs[
        (adata.obsm["spatial"][:, 1] >= min_row)
        & (adata.obsm["spatial"][:, 1] <= max_row)
        & (adata.obsm["spatial"][:, 0] >= min_col)
        & (adata.obsm["spatial"][:, 0] <= max_col)
    ].index
    return adata[cropped_cells, :].copy()

def preprocessing(adata, n_top_genes=20000):
    """Filter, normalize, select HVGs, scale and calculate PCA."""
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    q1 = adata.obs.total_counts.quantile(0.25)
    q3 = adata.obs.total_counts.quantile(0.75)
    iqr = q3 - q1
    adata = adata[adata.obs.total_counts <= (q3 + 1.5 * iqr), :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    n_vars_before_hvg = int(adata.n_vars)
    sc.pp.highly_variable_genes(
        adata,
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
        n_top_genes=n_top_genes,
    )
    n_hvg_true = int(adata.var.highly_variable.sum())
    hvg_names = adata.var_names[adata.var.highly_variable].astype(str)
    hvg_digest = sha256_bytes("\n".join(hvg_names).encode())
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata)
    sc.tl.pca(adata, svd_solver="arpack", random_state=SEED)
    audit = {
        "n_vars_before_hvg_subset": n_vars_before_hvg,
        "n_highly_variable_true_before_subset": n_hvg_true,
        "n_vars_after_hvg_subset": int(adata.n_vars),
        "all_retained_var_highly_variable": bool(adata.var.highly_variable.all()),
        "hvg_name_order_sha256": hvg_digest,
        "X_pca_shape": list(adata.obsm["X_pca"].shape),
        "X_pca_sha256": array_digest(adata.obsm["X_pca"]),
    }
    return adata, audit

def make_adata(row, visium_root: Path, annotation_root: Path, n_top_genes: int):
    sample = str(row["sample"])
    sample_sub = str(row["sample_sub"])
    min_row, max_row = row["min_row"], row["max_row"]
    min_col, max_col = row["min_col"], row["max_col"]
    adata = sc.read_visium(f"{visium_root}/{sample}/Spatial_matrix/")
    if not (
        math.isnan(min_row)
        and math.isnan(max_row)
        and math.isnan(min_col)
        and math.isnan(max_col)
    ):
        adata = crop(adata, min_row, max_row, min_col, max_col)
    adata, audit = preprocessing(adata, n_top_genes=n_top_genes)
    key = sample_sub[:-1] if sample_sub[-1] == "P" else sample_sub
    annotation = pd.read_csv(
        annotation_root / f"Annotation_professor_{key}.csv", index_col=0
    )
    annotation = annotation.dropna()
    annotation.columns = ["annotation"]
    annotation = annotation[~annotation.index.duplicated(keep="first")]
    adata.obs = adata.obs.merge(
        annotation, left_index=True, right_index=True, how="left"
    )
    adata.obs["sample"] = sample_sub
    adata = adata[
        adata.obs["annotation"].notna() & (adata.obs["annotation"] != "Exclude")
    ].copy()
    audit.update(
        {
            "sample": sample,
            "sample_sub": sample_sub,
            "n_obs_after_annotation_filter": int(adata.n_obs),
            "obs_names_unique": bool(adata.obs_names.is_unique),
        }
    )
    return adata, audit

def adata_intersection(dict_adata):
    """Select the set intersection of genes across samples."""
    list_adata = [adatum.var_names for adatum in dict_adata.values()]
    list_var_names = list(set.intersection(*map(set, list_adata)))
    return {
        sample_sub: adatum[:, list_var_names]
        for sample_sub, adatum in dict_adata.items()
    }

def integration(dict_adata):
    """Concatenate sample matrices and calculate the BBKNN graph."""
    dict_adata = {k: v for k, v in dict_adata.items() if v.shape[0] >= 10}
    adata = concat(
        list(dict_adata.values()),
        join="inner",
        label="sample",
        keys=dict_adata.keys(),
    )
    sc.external.pp.bbknn(adata, batch_key="sample")
    return adata
