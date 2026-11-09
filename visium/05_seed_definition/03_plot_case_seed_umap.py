#!/usr/bin/env python3
"""Plot case-level SEED assignments using X_umap."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc


COLORS = {"Other malignant spots": "#ffe082", "SEED": "#d7191c"}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def filename_safe(value: str) -> str:
    return "".join(char if char.isalnum() or char in "._" else "_" for char in value)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--case-id", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--formats", default="png,pdf,svg")
    parser.add_argument("--point-size", type=float, default=12.0)
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)
    required = {"case", "is_malignant_for_seed", "is_seed", "seed_group"}
    missing = sorted(required.difference(adata.obs.columns))
    if missing:
        raise KeyError(f"Missing required assignment columns: {missing}")
    if "X_umap" not in adata.obsm:
        raise KeyError("Input object has no stored X_umap coordinates")
    case_values = adata.obs["case"].astype(str)
    if set(case_values) != {args.case_id}:
        raise ValueError(f"Input cases {sorted(set(case_values))} do not equal {args.case_id}")
    xy = np.asarray(adata.obsm["X_umap"])
    if xy.ndim != 2 or xy.shape[1] < 2 or not np.isfinite(xy[:, :2]).all():
        raise ValueError("Stored UMAP coordinates are missing or non-finite")

    malignant = adata.obs["is_malignant_for_seed"].astype(bool).to_numpy()
    seed = adata.obs["is_seed"].astype(bool).to_numpy()
    shown = malignant | seed
    labels = np.full(adata.n_obs, "not_shown", dtype=object)
    labels[malignant & ~seed] = "Other malignant spots"
    labels[seed] = "SEED"

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{filename_safe(args.case_id)}_case_level_seed_umap"
    outputs = []
    fig, ax = plt.subplots(figsize=(8, 6))
    for label in ["Other malignant spots", "SEED"]:
        mask = labels == label
        ax.scatter(
            xy[mask, 0], xy[mask, 1], s=args.point_size,
            color=COLORS[label], linewidths=0, rasterized=True,
            label=f"{label} (n={int(mask.sum()):,})",
        )
    primary_count = int((seed & adata.obs["tissue"].astype(str).eq("primary").to_numpy()).sum())
    ln_count = int((seed & adata.obs["tissue"].astype(str).eq("lymph_node").to_numpy()).sum())
    radii_label = ", ".join(f"{float(r):g}" for r in adata.obs["seed_radius"].dropna().unique())
    plot_label = f"Primary SEED and LN analysis spots\nPrimary={primary_count}; LN={ln_count}; total={int(seed.sum())}; strict distance < {radii_label}"
    ax.set_title(f"{args.case_id}: {plot_label}")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_aspect("equal", adjustable="datalim")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    fig.tight_layout()
    for extension in [item.strip().lower() for item in args.formats.split(",") if item.strip()]:
        path = output_dir / f"{stem}.{extension}"
        fig.savefig(path, dpi=args.dpi, bbox_inches="tight")
        outputs.append({"path": str(path), "sha256": file_sha256(path)})
    plt.close(fig)

    metadata = {
        "case_id": args.case_id,
        "input_object": args.input,
        "plot_label": plot_label,
        "stored_umap_used": True,
        "pca_recomputed": False,
        "umap_recomputed": False,
        "seed_recomputed": False,
        "n_spots": int(adata.n_obs),
        "n_shown": int(shown.sum()),
        "n_other_malignant": int((labels == "Other malignant spots").sum()),
        "n_seed": int((labels == "SEED").sum()),
        "n_primary_seed": primary_count, "n_LN_seed": ln_count,
        "membership_stage": "distance-selected membership",
        "outputs": outputs,
    }
    metadata_path = output_dir / f"{stem}.metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
