#!/usr/bin/env python3
"""Prepare LuCA tumor counts, BBKNN/ridge clusters and downstream inputs."""
import argparse
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def prepare_counts(obj, gene_map):
    if not obj.obs_names.is_unique or not obj.var_names.is_unique:
        raise ValueError('Input cell and Ensembl IDs must be unique')
    for column in ['sample', 'dataset', 'uicc_stage']:
        if column not in obj.obs:
            raise ValueError(f'Missing cell metadata: {column}')
    if 'count' not in obj.layers:
        raise ValueError('Expected raw UMI counts in layers["count"]')
    counts = obj.layers['count']
    values = counts.data if sparse.issparse(counts) else np.asarray(counts)
    if not np.isfinite(values).all() or (values < 0).any() or not np.allclose(values, np.rint(values), rtol=0, atol=1e-6):
        raise ValueError('count must contain finite, nonnegative integer values')
    symbols = pd.Series([gene_map.get(gene, gene) for gene in obj.var_names], index=obj.var_names)
    keep = symbols.notna().to_numpy()
    result = ad.AnnData(X=counts[:, keep].copy(), obs=obj.obs.copy(), var=obj.var.iloc[np.flatnonzero(keep)].copy())
    result.var['ensembl_id'] = result.var_names.astype(str)
    result.var['gene_symbol'] = symbols.iloc[np.flatnonzero(keep)].astype(str).to_numpy()
    result.var_names = result.var['gene_symbol'].to_numpy()
    if not result.var_names.is_unique:
        raise ValueError('Gene symbols must be unique; supply an unambiguous mapping')
    result.layers['count'] = result.X.copy()
    return result


def align_cluster_metadata(obj, path, cluster_col):
    metadata = pd.read_csv(path, index_col=0, dtype={cluster_col: str}, low_memory=False)
    if not metadata.index.is_unique or set(metadata.index) != set(obj.obs_names):
        raise ValueError('Cluster metadata must contain exactly the input cell IDs, once each')
    metadata = metadata.loc[obj.obs_names]
    if cluster_col not in metadata or metadata[cluster_col].isna().any():
        raise ValueError(f'Complete {cluster_col} labels are required')
    for column in ['sample', 'dataset', 'uicc_stage']:
        if column not in metadata or not np.array_equal(metadata[column].astype(object).fillna('').astype(str), obj.obs[column].astype(object).fillna('').astype(str)):
            raise ValueError(f'Cluster metadata disagrees with input {column}')
    return metadata[cluster_col].astype(str)


def integrate(obj, resolution):
    import bbknn
    obj = obj.copy()
    sc.pp.normalize_total(obj, target_sum=1e4)
    sc.pp.log1p(obj)
    obj.raw = obj
    sc.pp.highly_variable_genes(obj, flavor='seurat_v3', n_top_genes=2000,
                               layer='count', batch_key='dataset', subset=True)
    sc.pp.scale(obj, max_value=10)
    sc.tl.pca(obj)
    sc.pp.neighbors(obj, n_pcs=50)
    sc.tl.umap(obj)
    bbknn.bbknn(obj, batch_key='dataset')
    sc.tl.umap(obj)
    sc.tl.leiden(obj, resolution=0.5)
    bbknn.ridge_regression(obj, batch_key=['dataset'], confounder_key=['leiden'])
    sc.pp.pca(obj)
    bbknn.bbknn(obj, batch_key='dataset')
    sc.tl.umap(obj)
    sc.tl.leiden(obj, resolution=resolution, key_added=f'leiden_res_{resolution:g}')
    return obj


def write_counts_h5(obj, path):
    counts = sparse.csr_matrix(obj.layers['count'])
    counts.eliminate_zeros()
    string = h5py.string_dtype('utf-8')
    with h5py.File(path, 'w') as handle:
        matrix = handle.create_group('matrix')
        # Cells x genes CSR arrays represent genes x cells CSC after transposition.
        for key, value in [('data', counts.data), ('indices', counts.indices), ('indptr', counts.indptr)]:
            matrix.create_dataset(key, data=value, compression='gzip')
        matrix.create_dataset('shape', data=np.array([obj.n_vars, obj.n_obs], dtype=np.int64))
        matrix.create_dataset('barcodes', data=obj.obs_names.to_numpy(dtype=object), dtype=string)
        features = matrix.create_group('features')
        for key, value in [('id', obj.var['ensembl_id']), ('name', obj.var_names),
                           ('feature_type', ['Gene Expression'] * obj.n_vars)]:
            features.create_dataset(key, data=np.asarray(value, dtype=object), dtype=string)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True, help='Tumor-only LuCA H5AD with raw count layer')
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--gene-map', type=Path, help='Saved TSV with ensembl_id and symbol; otherwise query MyGene.info')
    parser.add_argument('--resolution', type=float, default=1.1)
    parser.add_argument('--cluster-metadata', type=Path, help='Cell-indexed CSV supplying downstream cluster labels')
    parser.add_argument('--prepare-only', action='store_true', help='Export counts and supplied labels without computing a new graph')
    args = parser.parse_args()
    if args.prepare_only and args.cluster_metadata is None:
        parser.error('--prepare-only requires --cluster-metadata')
    if not np.isfinite(args.resolution) or args.resolution <= 0:
        parser.error('--resolution must be positive')
    if args.output_dir.exists():
        parser.error('Choose a new output directory')
    obj = ad.read_h5ad(args.input)
    if args.gene_map:
        table = pd.read_csv(args.gene_map, sep='\t')
        if not {'ensembl_id', 'symbol'} <= set(table) or table['ensembl_id'].duplicated().any():
            parser.error('Gene map requires unique ensembl_id rows and a symbol column')
        if not set(obj.var_names) <= set(table['ensembl_id']):
            parser.error('Gene map must cover every input feature, including unmapped rows')
        mapping = dict(zip(table['ensembl_id'], table['symbol']))
    else:
        import mygene
        hits = pd.DataFrame(mygene.MyGeneInfo().querymany(obj.var_names.tolist(), scopes='ensembl.gene', fields='symbol', species='human'))
        mapping = dict(zip(hits['query'], hits['symbol']))
    counts = prepare_counts(obj, mapping)
    cluster_col = f'leiden_res_{args.resolution:g}'
    supplied = align_cluster_metadata(counts, args.cluster_metadata, cluster_col) if args.cluster_metadata else None
    args.output_dir.mkdir(parents=True)
    pd.DataFrame({'ensembl_id': obj.var_names, 'symbol': [mapping.get(gene, gene) for gene in obj.var_names]}).to_csv(args.output_dir/'gene_mapping.tsv', sep='\t', index=False)
    del obj
    integrated = None if args.prepare_only else integrate(counts, args.resolution)
    if integrated is not None:
        counts.obs = integrated.obs.copy()
        if supplied is not None:
            counts.obs[f'computed_{cluster_col}'] = counts.obs[cluster_col].copy()
    if supplied is not None:
        counts.obs[cluster_col] = pd.Categorical(supplied.loc[counts.obs_names])
    if integrated is not None:
        integrated.obs = counts.obs.copy()
        integrated.write_h5ad(args.output_dir/'LuCA_Ref_Tumor_BBKNN.h5ad', compression='gzip')
    counts.obs.index.name = None
    counts.write_h5ad(args.output_dir/'LuCA_Ref_Tumor_counts.h5ad', compression='gzip')
    counts.obs.to_csv(args.output_dir/'LuCA_Ref_Tumor.obs.csv')
    write_counts_h5(counts, args.output_dir/'LuCA_Ref_Tumor_counts.h5')
    summary = {'cells': counts.n_obs, 'genes': counts.n_vars, 'cluster_column': cluster_col,
               'cluster_source': 'supplied_metadata' if supplied is not None else 'computed',
               'graph_computed': integrated is not None,
               'cluster_sizes': counts.obs[cluster_col].value_counts().to_dict()}
    (args.output_dir/'input_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
