#!/usr/bin/env python
"""Plot SEED and pseudotime-early labels using X_umap."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc


FIGURE3_ORDER = ["Malignant cells", "Malignant pseudotime early", "SEED"]
FIGURE3_COLORS = {
    "Malignant cells": "#ffe082",
    "Malignant pseudotime early": "#2166ac",
    "SEED": "#d7191c",
}
EARLY_ORDER = ["Other malignant spots", "Malignant pseudotime early"]
EARLY_COLORS = {
    "Other malignant spots": "#ffe082",
    "Malignant pseudotime early": "#2166ac",
}


def safe_name(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._" else "_" for ch in value)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-dir", default="results/seed/umap")
    parser.add_argument("--case-id", default="LUAD-SMC-103")
    parser.add_argument("--case-key", default="public_case_id")
    parser.add_argument("--point-size", type=float, default=15.0)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--formats", default="png,pdf,svg")
    return parser.parse_args()


def plot_stored(
    xy: np.ndarray,
    labels: np.ndarray,
    order: list[str],
    colors: dict[str, str],
    title: str,
    stem: Path,
    formats: list[str],
    point_size: float,
    dpi: int,
) -> list[Path]:
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    for category in order:
        mask = labels == category
        if not mask.any():
            continue
        size = point_size * (1.25 if category == "SEED" else 1.0)
        ax.scatter(xy[mask, 0], xy[mask, 1], c=colors[category], s=size,
                   alpha=0.95, linewidths=0, rasterized=True,
                   label=f"{category} (n={int(mask.sum()):,})")
    ax.set_title(title)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_aspect("equal", adjustable="datalim")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=8)
    fig.tight_layout()
    written = []
    for extension in formats:
        target = stem.with_suffix(f".{extension}")
        fig.savefig(target, dpi=dpi if extension == "png" else None, bbox_inches="tight")
        written.append(target)
    plt.close(fig)
    return written


def main() -> None:
    args = parse_args()
    formats = [value.strip().lower() for value in args.formats.split(",") if value.strip()]
    if not formats or not set(formats).issubset({"png", "pdf", "svg"}):
        raise ValueError("--formats must be a comma-separated subset of png,pdf,svg")
    adata = sc.read_h5ad(args.input)
    required = {args.case_key, "figure3_group", "pseudotime_early_group", "figure3_display_scope", "is_seed", "is_pseudotime_early"}
    missing = sorted(required.difference(adata.obs.columns))
    if missing:
        raise KeyError(f"Stored display assignments are missing: {missing}")
    if "X_umap" not in adata.obsm:
        raise KeyError("Input lacks stored X_umap")
    # Check that display labels match the input assignments.
    provenance = adata.uns.get("seed_display_provenance")
    if provenance is None:
        raise ValueError("Display provenance missing: rerun 04_attach_pseudotime_early.py after SEED definition")
    def digest_ids(mask):
        return hashlib.sha256("\n".join(sorted(adata.obs_names[mask].astype(str))).encode()).hexdigest()
    seed_mask = adata.obs["is_seed"].astype(bool).to_numpy()
    early_mask = adata.obs["is_pseudotime_early"].astype(bool).to_numpy()
    scope_mask = adata.obs["figure3_display_scope"].astype(bool).to_numpy()
    for key, selection in [("seed_id_sha256", seed_mask), ("early_id_sha256", early_mask), ("scope_id_sha256", scope_mask)]:
        if provenance.get(key) != digest_ids(selection):
            raise ValueError("Stale display provenance: " + key)
    if provenance.get("x_umap_sha256") != hashlib.sha256(np.ascontiguousarray(adata.obsm["X_umap"]).tobytes()).hexdigest():
        raise ValueError("UMAP changed after display assignments were generated")
    expected = np.full(adata.n_obs, "not_shown", dtype=object)
    expected[scope_mask] = "Malignant cells"
    expected[seed_mask] = "SEED"
    expected[early_mask] = "Malignant pseudotime early"
    if not np.array_equal(expected, adata.obs["figure3_group"].astype(str).to_numpy()):
        raise ValueError("Stored figure3_group disagrees with early-last display precedence")
    case_mask = adata.obs[args.case_key].astype(str).eq(args.case_id).to_numpy()
    scope = adata.obs["figure3_display_scope"].astype(bool).to_numpy()
    mask = case_mask & scope
    if not mask.any():
        raise ValueError(f"No malignant display spots found for {args.case_id}")
    xy = np.asarray(adata.obsm["X_umap"])[mask, :2]
    if not np.isfinite(xy).all():
        raise ValueError("Stored X_umap contains non-finite values")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = safe_name(args.case_id)
    figure3_labels = adata.obs.loc[mask, "figure3_group"].astype(str).to_numpy()
    early_labels = adata.obs.loc[mask, "pseudotime_early_group"].astype(str).to_numpy()
    radii_label = ", ".join(f"{float(r):g}" for r in adata.obs.loc[mask, "seed_radius"].dropna().unique())
    written = []
    written.extend(plot_stored(
        xy, figure3_labels, FIGURE3_ORDER, FIGURE3_COLORS,
        f"{args.case_id}: SEED and pseudotime-early groups\n"
        f"strict UMAP distance < {radii_label}; combined membership={int(seed_mask.sum())}, shown={int((seed_mask & ~early_mask).sum())}",
        output_dir / f"{prefix}_combined_pseudotime_seed_umap",
        formats, args.point_size, args.dpi,
    ))
    written.extend(plot_stored(
        xy, early_labels, EARLY_ORDER, EARLY_COLORS,
        f"{args.case_id}: pseudotime early",
        output_dir / f"{prefix}_pseudotime_early_umap",
        formats, args.point_size, args.dpi,
    ))
    print("Stored figure3 counts:", dict(zip(*np.unique(figure3_labels, return_counts=True))))
    print("Stored early-only counts:", dict(zip(*np.unique(early_labels, return_counts=True))))
    for path in written:
        print(path)


if __name__ == "__main__":
    main()
