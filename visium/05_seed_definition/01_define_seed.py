"""Define SEED spots from stored case UMAP coordinates.

Case thresholds are resolved by public case ID from a tracked TSV resource.
CLI overrides are permitted only when they are recorded as overrides in the
output metadata. Unknown cases never receive an implicit fallback threshold.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial.distance import cdist


STALE_DISPLAY_COLUMNS = [
    "group", "figure3_group", "figure3_display_scope", "pseudotime_early_group",
    "is_pseudotime_early", "is_seed_display",
]

DEFAULT_MALIGNANT_LABELS = [
    "Malignant_cell",
    "Malignant_cell_border_mixed",
]

MALIGNANT_LABEL_REPLACEMENTS = {
    "Malignant_cell(Hot)": "Malignant_cell",
    "Malignant_cell(Cold)": "Malignant_cell",
}

BROAD_LABEL_REPLACEMENTS = {
    **MALIGNANT_LABEL_REPLACEMENTS,
    "Malignant_cell_mucin": "Malignant_cell",
    "Malignant_mucin": "Malignant_cell",
    "Malignant_inflammatory_mixed": "Malignant_cell",
}


def parse_key_values(value: str | None) -> dict[str, float]:
    if not value:
        return {}
    parsed: dict[str, float] = {}
    for item in value.split(","):
        key, radius = item.split("=", 1)
        parsed[key.strip()] = float(radius)
    return parsed


def load_threshold_config(path: str | Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {
        "public_case_id", "threshold", "comparison_operator",
        "provenance_status", "case_status", "notes",
    }
    missing = sorted(required.difference(table.columns))
    if missing:
        raise ValueError(f"Threshold configuration is missing columns: {missing}")
    if table["public_case_id"].duplicated().any():
        duplicated = sorted(table.loc[table["public_case_id"].duplicated(False), "public_case_id"].unique())
        raise ValueError(f"Threshold configuration has duplicate case IDs: {duplicated}")
    unsupported = sorted(set(table["comparison_operator"]) - {"<"})
    if unsupported:
        raise ValueError(f"Only strict '<' comparison is supported; found: {unsupported}")
    unsupported_status = sorted(set(table["case_status"]) - {"DEFINED", "ANCHOR_ABSENT"})
    if unsupported_status:
        raise ValueError(f"Unsupported case_status values: {unsupported_status}")
    return table.set_index("public_case_id", drop=False)


def resolve_thresholds(
    cases: list[str],
    table: pd.DataFrame,
    overrides: dict[str, float],
) -> tuple[dict[str, float], set[str], dict[str, dict[str, object]]]:
    unexpected_overrides = sorted(set(overrides).difference(cases))
    if unexpected_overrides:
        raise ValueError(f"Threshold overrides were supplied for unrequested cases: {unexpected_overrides}")
    radii: dict[str, float] = {}
    anchor_absent_cases: set[str] = set()
    resolution: dict[str, dict[str, object]] = {}
    for case in cases:
        configured = table.loc[case] if case in table.index else None
        configured_threshold = "" if configured is None else configured["threshold"].strip()
        notes = "" if configured is None else configured["notes"]
        anchor_absent = configured is not None and configured["case_status"] == "ANCHOR_ABSENT"
        if anchor_absent:
            if configured_threshold:
                raise ValueError(f"{case} is ANCHOR_ABSENT but has a configured threshold")
            if case in overrides:
                raise ValueError(f"{case} is configured as anchor-absent; a threshold override is not allowed")
            anchor_absent_cases.add(case)
            resolution[case] = {
                "configured_threshold": None,
                "resolved_threshold": None,
                "comparison_operator": "<",
                "threshold_source": "anchor_absent_config",
                "provenance_status": configured["provenance_status"],
                "status": "NOT DEFINED — NO QUALIFYING LN ANCHOR",
                "notes": notes,
            }
            continue
        if case in overrides:
            radius = float(overrides[case])
            source = "cli_override"
        elif configured_threshold:
            radius = float(configured_threshold)
            source = "tracked_configuration"
        else:
            raise ValueError(
                f"No configured threshold exists for {case}. Add a tracked configuration row "
                "or supply an explicit --radii override."
            )
        radii[case] = radius
        resolution[case] = {
            "configured_threshold": float(configured_threshold) if configured_threshold else None,
            "resolved_threshold": radius,
            "comparison_operator": "<",
            "threshold_source": source,
            "provenance_status": "USER_OVERRIDE" if source == "cli_override" else configured["provenance_status"],
            "status": "DEFINED",
            "notes": notes,
        }
    return radii, anchor_absent_cases, resolution


def parse_replacements(value: str | None, mode: str) -> dict[str, str]:
    if mode == "none":
        replacements: dict[str, str] = {}
    elif mode == "broad":
        replacements = BROAD_LABEL_REPLACEMENTS.copy()
    else:
        replacements = MALIGNANT_LABEL_REPLACEMENTS.copy()
    if not value:
        return replacements
    for item in value.split(","):
        source, target = item.split("=", 1)
        replacements[source.strip()] = target.strip()
    return replacements


def infer_case(obs: pd.DataFrame, sample_key: str, explicit_case_key: str | None) -> pd.Series:
    if explicit_case_key and explicit_case_key in obs:
        return obs[explicit_case_key].astype(str)
    raise KeyError(
        "A public case column is required; implicit case parsing from sample IDs is disabled."
    )


def add_seed_labels(
    adata,
    annotation_key: str,
    sample_key: str,
    case_key: str | None,
    target_labels: list[str],
    anchor_tissue: str | None,
    malignant_labels: list[str],
    nearby_labels: list[str],
    label_replacements: dict[str, str],
    radii: dict[str, float],
    anchor_absent_cases: set[str],
) -> pd.DataFrame:
    if "X_umap" not in adata.obsm:
        raise KeyError("Integrated object must contain adata.obsm['X_umap'].")
    if annotation_key not in adata.obs:
        raise KeyError(f"Missing annotation column: {annotation_key}")
    if sample_key not in adata.obs:
        raise KeyError(f"Missing sample column: {sample_key}")

    obs = adata.obs.drop(columns=STALE_DISPLAY_COLUMNS, errors="ignore").copy()
    obs["seed_case"] = infer_case(obs, sample_key, case_key)
    labels = obs[annotation_key].astype(str).replace(label_replacements)
    obs["seed_annotation"] = labels
    is_anchor = labels.isin(target_labels)
    if anchor_tissue and "tissue" in obs:
        is_anchor = is_anchor & obs["tissue"].astype(str).eq(anchor_tissue)
    is_malignant = labels.isin(malignant_labels) | is_anchor
    is_nearby_candidate = labels.isin(nearby_labels) & is_malignant

    seed = pd.Series(False, index=obs.index)
    min_distance = pd.Series(np.nan, index=obs.index, dtype=float)
    seed_radius = pd.Series(np.nan, index=obs.index, dtype=float)
    umap = pd.DataFrame(adata.obsm["X_umap"], index=obs.index, columns=["UMAP1", "UMAP2"])

    for case, case_obs in obs.groupby("seed_case", sort=False):
        case_index = case_obs.index
        anchor_index = case_index[is_anchor.loc[case_index]]
        if case in anchor_absent_cases:
            if len(anchor_index):
                raise ValueError(
                    f"{case} is configured as anchor-absent but {len(anchor_index)} qualifying anchors were found"
                )
            continue
        radius = radii.get(case)
        if radius is None:
            raise ValueError(f"No resolved threshold exists for {case}")
        # The analysis group includes nearby Malignant_cell spots from both tissues.
        candidate_index = case_index[is_nearby_candidate.loc[case_index] | is_anchor.loc[case_index]]
        if len(anchor_index) == 0 or len(candidate_index) == 0:
            continue

        distances = cdist(umap.loc[candidate_index].to_numpy(), umap.loc[anchor_index].to_numpy())
        nearest = distances.min(axis=1)
        min_distance.loc[candidate_index] = nearest
        seed_radius.loc[candidate_index] = radius
        seed.loc[candidate_index[nearest < radius]] = True
        seed.loc[anchor_index] = True

    obs["is_malignant_for_seed"] = is_malignant
    obs["is_nearby_candidate_for_seed"] = is_nearby_candidate
    obs["is_seed_anchor"] = is_anchor
    obs["is_ln_anchor"] = is_anchor & obs["tissue"].astype(str).eq("lymph_node")
    obs["seed_umap_min_distance"] = min_distance
    obs["seed_radius"] = seed_radius
    is_primary = obs["tissue"].astype(str).eq("primary") if "tissue" in obs else pd.Series(True, index=obs.index)
    is_lymph_node = obs["tissue"].astype(str).eq("lymph_node") if "tissue" in obs else pd.Series(False, index=obs.index)
    obs["seed_group"] = np.select(
        [
            is_anchor & is_primary,
            is_anchor & is_lymph_node,
            is_anchor & ~is_primary & ~is_lymph_node,
            seed & ~is_anchor & is_primary,
            is_malignant & ~seed & is_primary,
            seed & ~is_anchor & is_lymph_node,
            is_malignant & ~seed & is_lymph_node,
        ],
        ["primary_anchor", "ln_anchor", "other_tissue_anchor", "primary_seed", "other_malignant", "ln_near_anchor_malignant", "ln_other_malignant"],
        default="other",
    )
    obs["is_seed"] = seed
    obs["is_primary_seed"] = seed & is_primary
    obs["is_metastatic_seed"] = seed  # Distance-selected membership, independent of display labels.
    obs["seed_membership_rule"] = "anchors_and_nearby_malignant_spots"
    obs["analysis_group"] = np.select(
        [
            seed & is_malignant,
            is_malignant & ~seed,
        ],
        ["combined_cells", "other_malignant_cells"],
        default="other",
    )
    return obs


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Apply public-case SEED thresholds to stored case UMAP coordinates."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-h5ad", default="results/intermediate/visium_seed_annotated.h5ad")
    parser.add_argument("--output-csv", default="results/tables/visium_seed_assignments.csv")
    parser.add_argument("--output-summary", default="results/tables/visium_seed_summary.json")
    parser.add_argument("--annotation-key", default="annotation")
    parser.add_argument("--sample-key", default="sample")
    parser.add_argument("--case-key", default="case")
    parser.add_argument(
        "--threshold-config",
        default="visium/resources/case_specific_seed_thresholds.tsv",
        help="Tracked public-case threshold table.",
    )
    parser.add_argument(
        "--samples",
        default="all",
        help="Optional explicit sample subset; default retains all samples of the selected case.",
    )
    parser.add_argument(
        "--cases",
        default="LUAD-SMC-103",
        help="Comma-separated cases to retain; use 'all' to keep every case.",
    )
    parser.add_argument(
        "--target-labels",
        default="Malignant_cell_border_mixed",
        help="Comma-separated anchor labels; the default does not restrict tissue.",
    )
    parser.add_argument("--anchor-tissue", default=None, help="Optional explicit restriction; the default has no tissue restriction.")
    parser.add_argument("--malignant-labels", default=",".join(DEFAULT_MALIGNANT_LABELS))
    parser.add_argument("--nearby-labels", default="Malignant_cell", help="Labels eligible as nearby malignant cells in primary AND LN.")
    parser.add_argument(
        "--label-replacement-mode",
        choices=["hot-cold", "broad", "none"],
        default="hot-cold",
        help="Hot-cold mode maps Hot/Cold malignant labels to Malignant_cell.",
    )
    parser.add_argument("--label-replacements", default=None, help="Comma-separated source=target replacements appended to mode defaults.")
    parser.add_argument(
        "--radii",
        default=None,
        help="Explicit public-case overrides; each is recorded as cli_override in output metadata.",
    )
    args = parser.parse_args()

    malignant_labels = [x.strip() for x in args.malignant_labels.split(",") if x.strip()]
    nearby_labels = [x.strip() for x in args.nearby_labels.split(",") if x.strip()]
    label_replacements = parse_replacements(args.label_replacements, args.label_replacement_mode)

    adata = sc.read_h5ad(args.input)
    if not adata.obs_names.is_unique:
        raise ValueError("Input obs_names must be sample-qualified and unique.")
    if "canonical_spot_id" in adata.obs:
        canonical = pd.Index(adata.obs["canonical_spot_id"].astype(str))
        if not canonical.equals(pd.Index(adata.obs_names.astype(str))):
            raise ValueError("obs_names and obs['canonical_spot_id'] disagree.")

    requested_samples = [x.strip() for x in args.samples.split(",") if x.strip()]
    if requested_samples and requested_samples != ["all"]:
        adata = adata[
            adata.obs[args.sample_key].astype(str).isin(requested_samples)
        ].copy()
        missing_samples = sorted(set(requested_samples) - set(adata.obs[args.sample_key].astype(str)))
        if missing_samples:
            raise ValueError(f"Requested samples missing from input: {missing_samples}")

    requested_cases = [x.strip() for x in args.cases.split(",") if x.strip()]
    if requested_cases and requested_cases != ["all"]:
        case_values = infer_case(adata.obs, args.sample_key, args.case_key)
        adata = adata[case_values.isin(requested_cases)].copy()
        retained_cases = set(infer_case(adata.obs, args.sample_key, args.case_key))
        missing_cases = sorted(set(requested_cases) - retained_cases)
        if missing_cases:
            raise ValueError(f"Requested cases missing from input: {missing_cases}")

    retained_cases = sorted(infer_case(adata.obs, args.sample_key, args.case_key).unique().tolist())
    threshold_table = load_threshold_config(args.threshold_config)
    radii, anchor_absent_cases, threshold_resolution = resolve_thresholds(
        retained_cases, threshold_table, parse_key_values(args.radii)
    )

    # Derived display groups from an earlier assignment must never survive rerunning SEED.
    removed_display_columns = [c for c in STALE_DISPLAY_COLUMNS if c in adata.obs]
    adata.obs = adata.obs.drop(columns=STALE_DISPLAY_COLUMNS, errors="ignore")
    adata.uns.pop("seed_display_provenance", None)
    obs = add_seed_labels(
        adata,
        annotation_key=args.annotation_key,
        sample_key=args.sample_key,
        case_key=args.case_key,
        target_labels=[x.strip() for x in args.target_labels.split(",") if x.strip()],
        anchor_tissue=args.anchor_tissue,
        malignant_labels=malignant_labels,
        nearby_labels=nearby_labels,
        label_replacements=label_replacements,
        radii=radii,
        anchor_absent_cases=anchor_absent_cases,
    )
    obs["seed_threshold_source"] = obs["seed_case"].map(
        {case: item["threshold_source"] for case, item in threshold_resolution.items()}
    )
    obs["seed_threshold_provenance_status"] = obs["seed_case"].map(
        {case: item["provenance_status"] for case, item in threshold_resolution.items()}
    )
    for column in [
        "seed_case",
        "seed_annotation",
        "is_malignant_for_seed",
        "is_nearby_candidate_for_seed",
        "is_seed_anchor",
        "is_ln_anchor",
        "is_metastatic_seed",
        "seed_membership_rule",
        "seed_umap_min_distance",
        "seed_radius",
        "seed_group",
        "is_seed",
        "is_primary_seed",
        "analysis_group",
        "seed_threshold_source",
        "seed_threshold_provenance_status",
    ]:
        adata.obs[column] = obs[column]

    Path(args.output_h5ad).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output_csv).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output_summary).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output_h5ad)
    if "canonical_spot_id" in obs:
        output = obs.reset_index(drop=True)
    else:
        output = obs.reset_index(names="canonical_spot_id")
    output.to_csv(args.output_csv, index=False)
    summary = {
        "input": args.input,
        "samples": sorted(obs[args.sample_key].astype(str).unique().tolist()),
        "cases": sorted(obs["seed_case"].astype(str).unique().tolist()),
        "radius_operator": "<",
        "threshold_config": args.threshold_config,
        "threshold_resolution": threshold_resolution,
        "threshold_override_requested": bool(args.radii),
        "label_replacements": label_replacements,
        "eligible_labels": DEFAULT_MALIGNANT_LABELS,
        "n_spots": int(len(obs)),
        "n_ln_anchors": int(obs["is_ln_anchor"].sum()),
        "n_primary_seed": int((obs["is_seed"] & obs["tissue"].astype(str).eq("primary")).sum()),
        "n_ln_seed": int((obs["is_seed"] & obs["tissue"].astype(str).eq("lymph_node")).sum()),
        "n_anchors": int(obs["is_seed_anchor"].sum()),
        "n_primary_anchors": int((obs["is_seed_anchor"] & obs["tissue"].astype(str).eq("primary")).sum()),
        "n_nearby_primary_malignant": int(obs["seed_group"].eq("primary_seed").sum()),
        "n_nearby_ln_malignant": int(obs["seed_group"].eq("ln_near_anchor_malignant").sum()),
        "anchor_tissue_restriction": args.anchor_tissue,
        "candidate_tissue_restriction": None,
        "membership_stage": "distance-selected membership",
        "removed_stale_display_columns": removed_display_columns,
        "n_seed_total": int(obs["is_seed"].sum()),
        "spot_order_sha256": hashlib.sha256(
            "\n".join(obs.index.astype(str)).encode()
        ).hexdigest(),
    }
    Path(args.output_summary).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(obs["seed_group"].value_counts(dropna=False).to_string())
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
