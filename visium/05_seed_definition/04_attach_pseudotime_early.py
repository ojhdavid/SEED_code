#!/usr/bin/env python
"""Attach pseudotime-early labels from a trajectory table to a SEED object."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


LABEL_REPLACEMENTS = {
    "Malignant_cell(Hot)": "Malignant_cell",
    "Malignant_cell(Cold)": "Malignant_cell",
}


def sorted_id_sha256(values: pd.Series | pd.Index | list[str]) -> str:
    text = "\n".join(sorted(map(str, values)))
    return hashlib.sha256(text.encode()).hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="SEED-annotated case h5ad")
    parser.add_argument("--trajectory-csv", required=True)
    parser.add_argument("--trajectory-id-column", default="x")
    parser.add_argument("--output-h5ad", required=True)
    parser.add_argument("--output-mapping", required=True)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--source-obs-key", default="source_obs_name")
    parser.add_argument("--sample-key", default="public_sample_id")
    parser.add_argument("--primary-sample", default=None, help="Optional explicit restriction; the bare-barcode early lookup has no tissue restriction.")
    parser.add_argument("--case-id", default="LUAD-SMC-103")
    parser.add_argument("--annotation-key", default="annotation")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    adata = sc.read_h5ad(args.input)
    required = {
        args.source_obs_key, args.sample_key, args.annotation_key,
        "is_seed", "is_ln_anchor", "seed_group", "case", "tissue",
    }
    missing = sorted(required.difference(adata.obs.columns))
    if missing:
        raise KeyError(f"Input object lacks required obs columns: {missing}")
    if "X_umap" not in adata.obsm:
        raise KeyError("Input object lacks stored X_umap coordinates")
    trajectory = pd.read_csv(args.trajectory_csv)
    if args.trajectory_id_column not in trajectory:
        raise KeyError(f"Trajectory ID column not found: {args.trajectory_id_column}")
    input_ids = pd.Index(trajectory[args.trajectory_id_column].dropna().astype(str).drop_duplicates())

    obs = adata.obs.copy()
    source_ids = obs[args.source_obs_key].astype(str)
    if set(obs["case"].astype(str)) != {args.case_id}:
        raise ValueError("Attach the case-specific trajectory to exactly the requested case")
    primary = obs["tissue"].astype(str).eq("primary")
    annotation = obs[args.annotation_key].astype(str).replace(LABEL_REPLACEMENTS)
    id_match = source_ids.isin(input_ids)
    # The set intersection uses bare source barcodes in the malignant pool.
    early = annotation.eq("Malignant_cell") & id_match
    if args.primary_sample:
        early &= obs[args.sample_key].astype(str).eq(args.primary_sample)
    seed = obs["is_seed"].astype(bool)
    malignant = annotation.isin(["Malignant_cell", "Malignant_cell_border_mixed"])

    # Preserve the distance-selected membership. Display order is
    # default -> metastatic_seeds -> malignant_pseudotime_early (early wins).
    display_group = np.full(adata.n_obs, "not_in_malignant_scope", dtype=object)
    display_group[malignant.to_numpy()] = "Malignant_cells"
    display_group[seed.to_numpy()] = "metastatic_seeds"
    display_group[early.to_numpy()] = "malignant_pseudotime_early"
    figure_group = np.full(adata.n_obs, "not_shown", dtype=object)
    figure_group[malignant.to_numpy()] = "Malignant cells"
    figure_group[seed.to_numpy()] = "SEED"
    figure_group[early.to_numpy()] = "Malignant pseudotime early"
    adata.obs["group"] = display_group
    adata.obs["is_seed_display"] = (seed & ~early).to_numpy()
    early_group = np.where(early, "Malignant pseudotime early", "Other malignant spots")
    adata.obs["is_pseudotime_early"] = early.to_numpy()
    adata.obs["figure3_group"] = pd.Categorical(
        figure_group,
        categories=["Malignant cells", "Malignant pseudotime early", "SEED", "not_shown"],
        ordered=True,
    )
    adata.obs["pseudotime_early_group"] = pd.Categorical(
        early_group,
        categories=["Other malignant spots", "Malignant pseudotime early"],
        ordered=True,
    )
    adata.obs["figure3_display_scope"] = malignant.to_numpy()
    adata.uns["seed_display_provenance"] = {
        "seed_id_sha256": sorted_id_sha256(adata.obs_names[seed]),
        "early_id_sha256": sorted_id_sha256(adata.obs_names[early]),
        "scope_id_sha256": sorted_id_sha256(adata.obs_names[malignant]),
        "x_umap_sha256": hashlib.sha256(np.ascontiguousarray(adata.obsm["X_umap"]).tobytes()).hexdigest(),
        "category_precedence": "malignant_pseudotime_early overwrites metastatic_seeds",
        "trajectory_csv": str(args.trajectory_csv),
        "trajectory_sha256": hashlib.sha256(Path(args.trajectory_csv).read_bytes()).hexdigest(),
        "source_barcode_matching": "set intersection, no deduplication of object rows",
    }

    lookup = pd.DataFrame({
        "canonical_spot_id": adata.obs_names.astype(str),
        "source_obs_name": source_ids.to_numpy(),
        "public_sample_id": obs[args.sample_key].astype(str).to_numpy(),
        "annotation": annotation.to_numpy(),
        "is_seed": seed.to_numpy(),
        "is_pseudotime_early": early.to_numpy(),
    })
    rows = []
    for trajectory_id in input_ids:
        matched = lookup.loc[lookup["source_obs_name"].eq(trajectory_id)]
        if matched.empty:
            rows.append({
                "trajectory_id": trajectory_id, "match_status": "unmatched",
                "canonical_spot_id": "", "public_sample_id": "", "annotation": "",
                "is_seed": False, "is_pseudotime_early": False,
            })
        else:
            for _, row in matched.iterrows():
                rows.append({"trajectory_id": trajectory_id, "match_status": "matched", **row.to_dict()})
    audit = pd.DataFrame(rows)
    Path(args.output_h5ad).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output_mapping).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output_summary).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output_h5ad)
    audit.to_csv(args.output_mapping, sep="\t", index=False)

    early_ids = adata.obs_names[early]
    summary = {
        "source_trajectory_file": args.trajectory_csv,
        "trajectory_id_column": args.trajectory_id_column,
        "n_input_ids": int(len(input_ids)),
        "n_matched_ids": int(audit.loc[audit["match_status"].eq("matched"), "trajectory_id"].nunique()),
        "n_unmatched_ids": int(audit.loc[audit["match_status"].eq("unmatched"), "trajectory_id"].nunique()),
        "n_pseudotime_early": int(early.sum()),
        "early_primary": int((early & primary).sum()),
        "early_lymph_node": int((early & ~primary).sum()),
        "overlap_with_seed": int((early & seed).sum()),
        "category_precedence": "malignant_pseudotime_early overwrites metastatic_seeds",
        "membership_unchanged": bool(np.array_equal(seed.to_numpy(), adata.obs["is_seed"].to_numpy())),
        "early_sample_restriction": args.primary_sample,
        "n_seed_before_early": int(seed.sum()),
        "n_seed_display_after_early": int((seed & ~early).sum()),
        "primary_seed_display_after_early": int((seed & ~early & primary).sum()),
        "LN_seed_display_after_early": int((seed & ~early & obs["tissue"].astype(str).eq("lymph_node")).sum()),
        "overlap_seed_early_primary": int((seed & early & primary).sum()),
        "overlap_seed_early_LN": int((seed & early & obs["tissue"].astype(str).eq("lymph_node")).sum()),
        "bare_barcode_duplicate_rows": int(source_ids.duplicated(False).sum()),
        "sorted_public_early_id_sha256": sorted_id_sha256(early_ids),
        "x_umap_recomputed": False,
    }
    Path(args.output_summary).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"[done] matched_ids={summary['n_matched_ids']} early_spots={summary['n_pseudotime_early']}")


if __name__ == "__main__":
    main()
