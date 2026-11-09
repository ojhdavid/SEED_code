#!/usr/bin/env python
"""Calculate and plot SEED meta-scores."""

from __future__ import annotations

import argparse
import copy
import hashlib
import inspect
import json
import platform
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy import sparse


RANDOM_SEED = 42
SCORE_MODES = ["integrated-x", "lognorm-unscaled"]
CANONICAL_SIGNATURE = Path(__file__).resolve().parent / "resources/luad_smc_103_seed_meta_signature.tsv"


def safe_name(value: object) -> str:
    return "".join(ch if ch.isalnum() or ch in "._" else "_" for ch in str(value))


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_signature(path: str | Path) -> tuple[list[str], pd.DataFrame]:
    table = pd.read_csv(path, sep="\t")
    if "gene" not in table:
        raise KeyError("Signature TSV must contain a gene column")
    genes = table["gene"].dropna().astype(str).tolist()
    if len(genes) != len(set(genes)):
        raise ValueError("Signature genes must be unique")
    if not genes or any(gene.strip() != gene or not gene for gene in genes):
        raise ValueError("Signature gene symbols must be nonempty and whitespace-free")
    return genes, table


def array_hash(values) -> str:
    values = np.ascontiguousarray(values)
    return hashlib.sha256(values.tobytes()).hexdigest()


def assignment_hash(obs: pd.DataFrame) -> str:
    columns = [key for key in ("is_seed", "is_ln_anchor", "seed_group", "figure3_group", "is_pseudotime_early") if key in obs]
    return hashlib.sha256(obs[columns].to_csv(index=True).encode()).hexdigest()


def matrix_fingerprint(matrix) -> dict[str, object]:
    sample = matrix[: min(matrix.shape[0], 1000), : min(matrix.shape[1], 1000)]
    if sparse.issparse(sample):
        sample = sample.toarray()
    sample = np.asarray(sample, dtype=float)
    return {
        "shape": list(matrix.shape),
        "sample_min": float(np.nanmin(sample)),
        "sample_max": float(np.nanmax(sample)),
        "sample_negative_fraction": float(np.mean(sample < 0)),
    }


def expression_object(adata: AnnData, args: argparse.Namespace) -> tuple[AnnData, dict[str, object]]:
    mode = args.resolved_score_mode
    if mode == "integrated-x":
        out = adata.copy()
        provenance = {
            "matrix_mode": mode,
            "matrix": "adata.X",
            "normalization_performed_by_this_script": False,
            "score_genes_use_raw": False,
            "normalization_scaling_state": "adata.X is used unchanged; the default workflow supplies per-sample normalized, log-transformed and scaled expression",
            "interpretation": "score_genes uses the supplied integrated expression matrix",
        }
        provenance.update(matrix_fingerprint(out.X))
        return out, provenance

    if bool(args.counts_layer) == bool(args.use_raw_counts):
        raise ValueError(
            "lognorm-unscaled mode requires exactly one of --counts-layer or --use-raw-counts"
        )
    if args.use_raw_counts:
        if adata.raw is None:
            raise ValueError("--use-raw-counts requested but adata.raw is absent")
        out = adata.raw.to_adata()
        out.obs = adata.obs.copy()
        out.obsm = {key: value.copy() for key, value in adata.obsm.items()}
        source = "adata.raw.X"
    else:
        if args.counts_layer not in adata.layers:
            raise KeyError(f"Counts layer not found: {args.counts_layer}")
        out = AnnData(
            X=adata.layers[args.counts_layer].copy(),
            obs=adata.obs.copy(),
            var=adata.var.copy(),
            obsm={key: value.copy() for key, value in adata.obsm.items()},
            uns=copy.deepcopy(adata.uns),
        )
        source = f"adata.layers[{args.counts_layer!r}]"
    values = out.X.data if sparse.issparse(out.X) else np.asarray(out.X).ravel()
    if (np.any(~np.isfinite(values)) or np.any(values < 0)
            or not np.allclose(values, np.rint(values), atol=1e-6)):
        raise ValueError("lognorm-unscaled requires verified finite, nonnegative integer raw counts")
    before = matrix_fingerprint(out.X)
    sc.pp.normalize_total(out, target_sum=10000)
    sc.pp.log1p(out)
    provenance = {
        "matrix_mode": mode,
        "matrix": source,
        "normalization_performed_by_this_script": True,
        "normalization": "scanpy.pp.normalize_total(target_sum=10000), then scanpy.pp.log1p",
        "score_genes_use_raw": False,
        "scaling_performed_by_this_script": False,
        "normalization_scaling_state": "raw counts normalized to 10000 then log1p; unscaled",
        "input_matrix_fingerprint": before,
    }
    provenance.update(matrix_fingerprint(out.X))
    return out, provenance


def plot_umap_by_case(adata: AnnData, case_key: str, outdir: Path, cmap: str,
                      point_size: float, dpi: int, score_mode: str) -> list[Path]:
    if "X_umap" not in adata.obsm:
        return []
    xy = np.asarray(adata.obsm["X_umap"])
    scores = adata.obs["meta_score"].astype(float)
    vmin, vmax = scores.quantile([0.01, 0.99])
    written = []
    for case_id, indexes in adata.obs.groupby(case_key, observed=True, sort=True).groups.items():
        positions = adata.obs_names.get_indexer(indexes)
        fig, ax = plt.subplots(figsize=(6.2, 4.9))
        points = ax.scatter(xy[positions, 0], xy[positions, 1], c=scores.loc[indexes],
                            cmap=cmap, s=point_size, linewidths=0, vmin=vmin, vmax=vmax,
                            rasterized=True)
        ax.set_title(f"{case_id}: SEED meta-score\n{score_mode}", fontsize=10)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.set_aspect("equal", adjustable="datalim")
        colorbar = plt.colorbar(points, ax=ax, fraction=0.035, pad=0.04)
        colorbar.set_label("meta_score (1–99% display)", fontsize=8)
        fig.tight_layout()
        path = outdir / f"{safe_name(case_id)}_metascore_{safe_name(score_mode)}_umap.png"
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        written.append(path)
    return written


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    signature_input = parser.add_mutually_exclusive_group(required=True)
    signature_input.add_argument("--signature", help="Signature TSV with a gene column")
    signature_input.add_argument("--use-paper-signature", action="store_true",
                                 help="Use the bundled 93-gene signature")
    parser.add_argument(
        "--score-mode", choices=SCORE_MODES, default=None,
        help=("Scoring profile. Defaults to integrated-x and scores integrated adata.X."),
    )
    parser.add_argument(
        "--matrix-mode", choices=SCORE_MODES, default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--counts-layer")
    parser.add_argument("--use-raw-counts", action="store_true")
    parser.add_argument("--output-h5ad")
    parser.add_argument("--output-scores")
    parser.add_argument("--output-summary")
    parser.add_argument("--output-provenance")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--case-key", default="public_case_id")
    parser.add_argument("--sample-key", default="public_sample_id")
    parser.add_argument("--cmap", default="viridis")
    parser.add_argument("--point-size", type=float, default=12)
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.score_mode and args.matrix_mode and args.score_mode != args.matrix_mode:
        raise ValueError("--score-mode and --matrix-mode disagree")
    args.score_mode_requested = args.score_mode or args.matrix_mode or "default"
    args.resolved_score_mode = args.score_mode or args.matrix_mode or "integrated-x"
    if args.resolved_score_mode == "integrated-x" and (args.counts_layer or args.use_raw_counts):
        raise ValueError(
            "integrated-x uses adata.X; do not supply --counts-layer or --use-raw-counts"
        )
    mode_token = safe_name(args.resolved_score_mode)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    defaults = {
        "output_h5ad": output_dir / f"metascore_{mode_token}.h5ad",
        "output_scores": output_dir / f"metascore_{mode_token}_scores.tsv",
        "output_summary": output_dir / f"metascore_{mode_token}_summary.tsv",
        "output_provenance": output_dir / f"metascore_{mode_token}_provenance.json",
    }
    for field, default in defaults.items():
        value = Path(getattr(args, field) or default)
        if mode_token not in value.stem:
            raise ValueError(
                f"{field.replace('_', '-')} filename must include resolved mode token {mode_token!r}: {value}"
            )
        setattr(args, field, str(value))
    input_path = Path(args.input).resolve()
    for field in defaults:
        target = Path(getattr(args, field)).resolve()
        if target == input_path:
            raise ValueError("A scoring output must not overwrite its input object")
        if target.exists():
            raise FileExistsError(f"Preserve an existing scoring output: {target}")
    adata_input = sc.read_h5ad(args.input)
    if not adata_input.obs_names.is_unique:
        raise ValueError("Input canonical spot IDs must be unique")
    for key in (args.case_key, args.sample_key):
        if key not in adata_input.obs:
            raise KeyError(f"Missing public metadata column: {key}")
    signature_source = "bundled" if args.use_paper_signature else "file"
    if args.use_paper_signature:
        args.signature = str(CANONICAL_SIGNATURE)
    genes, signature_table = read_signature(args.signature)
    input_umap_hash = array_hash(adata_input.obsm["X_umap"]) if "X_umap" in adata_input.obsm else None
    input_assignment_hash = assignment_hash(adata_input.obs)
    adata, provenance = expression_object(adata_input, args)
    present = [gene for gene in genes if gene in adata.var_names]
    missing = [gene for gene in genes if gene not in adata.var_names]
    if not present:
        raise ValueError("No signature genes are present in the selected expression matrix")
    score_parameters = {
        "ctrl_size": 50, "n_bins": 25, "score_name": "meta_score",
        "random_state": RANDOM_SEED, "use_raw": False,
    }
    if "ctrl_as_ref" in inspect.signature(sc.tl.score_genes).parameters:
        score_parameters["ctrl_as_ref"] = True
    sc.tl.score_genes(adata, gene_list=present, gene_pool=adata.var_names.tolist(), **score_parameters)
    output_umap_hash = array_hash(adata.obsm["X_umap"]) if "X_umap" in adata.obsm else None
    if output_umap_hash != input_umap_hash or assignment_hash(adata.obs) != input_assignment_hash:
        raise AssertionError("Scoring changed stored UMAP coordinates or SEED/early assignments")
    written = plot_umap_by_case(
        adata, args.case_key, output_dir, args.cmap, args.point_size, args.dpi,
        args.resolved_score_mode,
    )
    score_table = adata.obs[[args.case_key, args.sample_key, "meta_score"]].copy()
    score_table.insert(0, "canonical_spot_id", adata.obs_names.astype(str))
    summary = (
        score_table.groupby([args.case_key, args.sample_key], observed=False)["meta_score"]
        .agg(["count", "mean", "median", "std", "min", "max"])
        .reset_index()
    )
    provenance.update({
        "case_id": (
            str(adata.obs[args.case_key].astype(str).iloc[0])
            if adata.obs[args.case_key].astype(str).nunique() == 1 else None
        ),
        "case_ids": sorted(adata.obs[args.case_key].astype(str).unique().tolist()),
        "scoring_mode_requested": args.score_mode_requested,
        "scoring_mode_resolved": args.resolved_score_mode,
        "input_object": args.input,
        "signature_file": args.signature,
        "signature_file_sha256": hashlib.sha256(Path(args.signature).read_bytes()).hexdigest(),
        "signature_source": signature_source,
        "canonical_resource_read": bool(args.use_paper_signature),
        "signature_genes_used": present,
        "signature_genes_total": len(genes),
        "signature_genes_present": len(present),
        "signature_genes_missing": missing,
        "complete_gene_pool_size": int(adata.n_vars),
        "control_gene_pool": "all genes in selected scoring matrix",
        "score_genes_parameters": score_parameters,
        "scoring_method": "scanpy.tl.score_genes",
        "score_column": "meta_score",
        "spots_scored": int(adata.n_obs),
        "spot_population": "all spots in the supplied stored object; no SEED redefinition or spot selection",
        "input_object_sha256": file_sha256(args.input),
        "input_raw_present": adata_input.raw is not None,
        "input_layers": list(adata_input.layers),
        "stored_umap_sha256_before": input_umap_hash,
        "stored_umap_sha256_after": output_umap_hash,
        "stored_assignments_sha256": input_assignment_hash,
        "use_raw": False,
        "score_minimum": float(score_table["meta_score"].min()),
        "score_median": float(score_table["meta_score"].median()),
        "score_mean": float(score_table["meta_score"].mean()),
        "score_maximum": float(score_table["meta_score"].max()),
        "random_seed": RANDOM_SEED,
        "scanpy_version": sc.__version__,
        "python_version": platform.python_version(),
        "figures": [str(path) for path in written],
        "figure_display": {
            "coordinates": "stored X_umap; no recomputation",
            "cmap": args.cmap, "point_size": args.point_size, "dpi": args.dpi,
            "color_limits": "global score 1st and 99th percentiles; score values unchanged",
            "vmin": float(score_table["meta_score"].quantile(0.01)),
            "vmax": float(score_table["meta_score"].quantile(0.99)),
        },
        "discovery_provenance": signature_table.iloc[0].to_dict(),
    })
    for path in [args.output_h5ad, args.output_scores, args.output_summary, args.output_provenance]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output_h5ad)
    score_table.to_csv(args.output_scores, sep="\t", index=False)
    summary.to_csv(args.output_summary, sep="\t", index=False)
    gene_audit_path = output_dir / f"metascore_{mode_token}_signature_gene_coverage.tsv"
    pd.DataFrame({"signature_rank": np.arange(1, len(genes) + 1), "gene": genes,
                  "used_in_scoring": [gene in present for gene in genes]}).to_csv(gene_audit_path, sep="\t", index=False)
    provenance["output_hashes"] = {
        "h5ad": {"path": args.output_h5ad, "sha256": file_sha256(args.output_h5ad)},
        "scores": {"path": args.output_scores, "sha256": file_sha256(args.output_scores)},
        "summary": {"path": args.output_summary, "sha256": file_sha256(args.output_summary)},
        "gene_coverage": {"path": str(gene_audit_path), "sha256": file_sha256(gene_audit_path)},
        "figures": [{"path": str(path), "sha256": file_sha256(path)} for path in written],
    }
    Path(args.output_provenance).write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(f"[done] spots={adata.n_obs} signature_genes={len(present)}/{len(genes)} output={output_dir}")


if __name__ == "__main__":
    main()
