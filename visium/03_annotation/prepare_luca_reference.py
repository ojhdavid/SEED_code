#!/usr/bin/env python3
"""Convert LuCA core raw counts and cell metadata to CARD input formats."""
import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd

LUCA_CORE_URL = "https://datasets.cellxgene.cziscience.com/dd538ee7-f5e4-49e9-9f1e-2a1ea5246cf4.h5ad"
COLLECTION = "edb893ee-4066-4128-9aec-5eb2b03f8287"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        parser.error("Refusing to reuse nonempty output directory")
    obj = ad.read_h5ad(args.input_h5ad, backed="r")
    required = ["sample", "cell_type", "cell_type_major", "cell_type_tumor"]
    if obj.raw is None or "feature_name" not in obj.raw.var or not all(x in obj.obs for x in required):
        parser.error("Expected verified LuCA core raw counts and documented cell annotations")
    symbols = obj.raw.var["feature_name"].astype(str).to_numpy()
    ids = obj.raw.var_names.astype(str).to_numpy()
    cells = obj.obs_names.astype(str).to_numpy()
    if len(set(symbols)) != len(symbols) or len(set(cells)) != len(cells):
        parser.error("Gene symbols and cell IDs must be unique.")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    metadata = obj.obs[required].copy()
    metadata.insert(0, "cellID", cells)
    metadata = metadata.rename(columns={"sample": "sampleInfo"})
    if metadata.isna().any().any():
        parser.error("Missing annotation values")
    metadata.to_csv(args.output_dir / "CZI_core_reference_celltype_all.csv", index=False)
    pd.DataFrame({"gene_symbol": symbols, "feature_id": ids}).to_csv(args.output_dir / "reference_gene_mapping.tsv", sep="\t", index=False)
    target = args.output_dir / "CZI_core_reference_counts_10x.h5"
    with h5py.File(args.input_h5ad, "r") as source, h5py.File(target, "w") as output:
        raw = source["raw/X"]
        if raw.attrs["encoding-type"] != "csr_matrix":
            parser.error("Expected cell-by-gene CSR raw counts; refusing inferred matrix conversion")
        matrix = output.create_group("matrix")
        # Cell-by-gene CSR and its gene-by-cell CSC transpose share sparse arrays.
        n_nonzero = len(raw["data"])
        value_min, value_max = float("inf"), float("-inf")
        for name in ("data", "indices", "indptr"):
            input_array = raw[name]
            dest = matrix.create_dataset(name, shape=input_array.shape, dtype=input_array.dtype, compression="gzip", compression_opts=1)
            for begin in range(0, len(input_array), 2_000_000):
                values = input_array[begin:begin + 2_000_000]
                if name == "data":
                    if not np.all(np.isfinite(values)) or np.any(values < 0) or np.any(np.abs(values - np.rint(values)) > 1e-6):
                        parser.error("raw.X is not finite nonnegative integer counts")
                    value_min = min(value_min, float(values.min()))
                    value_max = max(value_max, float(values.max()))
                dest[begin:begin + len(values)] = values
        matrix.create_dataset("shape", data=np.array([len(symbols), len(cells)], dtype="int64"))
        string_type = h5py.string_dtype(encoding="utf-8")
        matrix.create_dataset("barcodes", data=cells, dtype=string_type)
        features = matrix.create_group("features")
        for name, values in {"name": symbols, "id": ids, "feature_type": np.repeat("Gene Expression", len(symbols))}.items():
            features.create_dataset(name, data=np.asarray(values, dtype=object), dtype=string_type)
    citation = str(obj.uns.get("citation", ""))
    obj.file.close()
    provenance = {"input": str(args.input_h5ad.resolve()), "input_sha256": sha256(args.input_h5ad),
                  "reference_asset": LUCA_CORE_URL, "collection": COLLECTION,
                  "publication_doi": "10.1016/j.ccell.2022.10.008", "embedded_citation": citation,
                  "matrix": "raw.X", "normalization": "none", "transpose": "cell-gene CSR to gene-cell CSC; sparse arrays unchanged",
                  "genes": len(symbols), "cells": len(cells), "nonzero": n_nonzero,
                  "count_min": value_min, "count_max": value_max, "gene_identifier": "raw.var.feature_name (symbols); feature IDs retained separately",
                  "cell_type_column": "cell_type_major", "sample_column": "sampleInfo from obs.sample",
                  "counts_10x_sha256": sha256(target)}
    (args.output_dir / "conversion_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps({k: provenance[k] for k in ("genes", "cells", "nonzero", "matrix", "count_min", "count_max")}))


if __name__ == "__main__":
    main()
