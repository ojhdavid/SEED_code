#!/usr/bin/env python
"""Map cohort sample IDs and plot UMAPs."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns


def array_hash(value) -> str:
    array = np.ascontiguousarray(value)
    digest = hashlib.sha256()
    digest.update(str(array.shape).encode())
    digest.update(str(array.dtype).encode())
    digest.update(array.tobytes())
    return digest.hexdigest()


def relabel(adata, mapping: pd.DataFrame) -> dict[str, object]:
    sample_map = dict(zip(mapping["source_sample_id"], mapping["public_sample_id"]))
    case_map = dict(zip(mapping["source_case_id"], mapping["public_case_id"]))
    source_sample = adata.obs["sample"].astype(str)
    source_case = adata.obs["case"].astype(str)
    if not set(source_sample).issubset(sample_map) or not set(source_case).issubset(case_map):
        raise ValueError("Stored object contains IDs absent from the explicit mapping table")
    before = {key: array_hash(adata.obsm[key]) for key in ["X_pca", "X_umap"]}
    source_obs = adata.obs["raw_barcode"].astype(str) if "raw_barcode" in adata.obs else pd.Series(
        [str(value).split("_", 1)[-1] for value in adata.obs_names], index=adata.obs_names
    )
    adata.obs["source_obs_name"] = source_obs.to_numpy()
    adata.obs["source_sample_id"] = source_sample.to_numpy()
    adata.obs["source_case_id"] = source_case.to_numpy()
    adata.obs["public_sample_id"] = source_sample.map(sample_map).to_numpy()
    adata.obs["public_case_id"] = source_case.map(case_map).to_numpy()
    adata.obs["sample"] = adata.obs["public_sample_id"].to_numpy()
    adata.obs["case"] = adata.obs["public_case_id"].to_numpy()
    adata.obs["canonical_spot_id"] = (
        adata.obs["public_sample_id"].astype(str) + "::" + adata.obs["source_obs_name"].astype(str)
    )
    if adata.obs["canonical_spot_id"].duplicated().any():
        raise ValueError("Public canonical spot IDs are not unique")
    adata.obs_names = adata.obs["canonical_spot_id"].to_numpy()
    after = {key: array_hash(adata.obsm[key]) for key in ["X_pca", "X_umap"]}
    if before != after:
        raise RuntimeError("Coordinate hash changed during metadata-only relabeling")
    return {"before": before, "after": after, "unchanged": True, "n_spots": adata.n_obs}


def plot(adata, key: str, output: Path, title: str, dpi: int) -> None:
    xy = np.asarray(adata.obsm["X_umap"])
    values = adata.obs[key].astype(str)
    categories = sorted(values.unique())
    colors = dict(zip(categories, sns.color_palette("tab20", len(categories)).as_hex()))
    fig, ax = plt.subplots(figsize=(10, 6))
    for category in categories:
        mask = values.eq(category).to_numpy()
        ax.scatter(xy[mask, 0], xy[mask, 1], s=5, color=colors[category], linewidths=0,
                   rasterized=True, label=f"{category} (n={mask.sum():,})")
    ax.set_title(title)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_aspect("equal", adjustable="datalim")
    fig.subplots_adjust(left=0.09, right=0.70, bottom=0.11, top=0.91)
    ax.legend(
        loc="center left", bbox_to_anchor=(0.72, 0.5), bbox_transform=fig.transFigure,
        frameon=False, fontsize=7,
    )
    fig.savefig(output, dpi=dpi)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--uncorrected-input", required=True)
    parser.add_argument("--bbknn-input", required=True)
    parser.add_argument("--sample-map", default="visium/resources/sample_id_mapping.tsv")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--object-dir", required=True)
    parser.add_argument("--audit-json", required=True)
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    mapping = pd.read_csv(args.sample_map, sep="\t", dtype=str, keep_default_na=False)
    output_dir, object_dir = Path(args.output_dir), Path(args.object_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    object_dir.mkdir(parents=True, exist_ok=True)
    audit = {}
    objects = {
        "uncorrected": sc.read_h5ad(args.uncorrected_input),
        "bbknn": sc.read_h5ad(args.bbknn_input),
    }
    for name, adata in objects.items():
        audit[name] = relabel(adata, mapping)
        adata.write_h5ad(object_dir / f"joint_umap_{name}_public_ids.h5ad")
    plot(objects["uncorrected"], "public_case_id", output_dir / "uncorrected_umap_by_case.png",
         "Uncorrected whole-cohort UMAP by case", args.dpi)
    plot(objects["uncorrected"], "public_sample_id", output_dir / "uncorrected_umap_by_sample.png",
         "Uncorrected whole-cohort UMAP by sample", args.dpi)
    plot(objects["bbknn"], "public_case_id", output_dir / "bbknn_umap_by_case.png",
         "BBKNN whole-cohort UMAP by case", args.dpi)
    plot(objects["bbknn"], "public_sample_id", output_dir / "bbknn_umap_by_sample.png",
         "BBKNN whole-cohort UMAP by sample", args.dpi)
    Path(args.audit_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.audit_json).write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(f"[done] objects={len(objects)} figures=4 output={output_dir}")


if __name__ == "__main__":
    main()
