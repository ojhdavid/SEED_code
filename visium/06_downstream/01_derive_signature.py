#!/usr/bin/env python3
"""SEED group selection, raw-count assembly and marker-filtered signature derivation.

Called by 03_raw_to_signature.py; this module does not launch a separate workflow.
"""
import anndata as ad
import numpy as np
import pandas as pd

def restore_raw_counts(adata_selected, raw_data_dict):
    # Join raw counts on case barcodes and retain the shared gene order.
    selected_barcodes=adata_selected.obs_names
    sample_labels=adata_selected.obs['sample'].unique()
    adata_selected.obs['sample_sub']=adata_selected.obs['sample']
    raw_data_list=[];all_common_genes=set(adata_selected.var_names)
    for sample in sample_labels:
        if sample in raw_data_dict:all_common_genes &= set(raw_data_dict[sample].to_df().columns)
    common_genes=sorted(all_common_genes)
    for sample in sample_labels:
        if sample in raw_data_dict:
            raw_adata=raw_data_dict[sample]
            matching_barcodes=selected_barcodes.intersection(raw_adata.obs_names)
            if len(matching_barcodes)>0:
                raw_df_filtered=raw_adata.to_df().loc[matching_barcodes,common_genes]
                raw_df_filtered=raw_df_filtered[~raw_df_filtered.index.duplicated(keep='first')]
                raw_df_filtered=raw_df_filtered.loc[:,~raw_df_filtered.columns.duplicated()]
                raw_df_filtered=raw_df_filtered[common_genes]
                raw_data_list.append(raw_df_filtered)
    if not raw_data_list:raise ValueError('No matching raw-count rows found')
    raw_data_selected=pd.concat(raw_data_list,axis=0)
    restored=ad.AnnData(X=raw_data_selected.values,obs=adata_selected.obs.loc[raw_data_selected.index].copy(),var=adata_selected.var.loc[common_genes].copy())
    restored.obsm=adata_selected.obsm.copy()
    return restored

def process_alk_sample_and_extract_meta_markers(sample_name,integrated_adata,sample_metadata_dict,extract_raw_data_func,malignant_annotations,coord_csv_path,top_n_genes=40,stromal_marker_list=None,lymphocyte_marker_list=None,macrophage_marker_list=None,endothelial_marker_list=None,additional_remove_list=None):
 """Select distance-based groups and derive a marker-filtered DEG signature."""
 import anndata as ad,numpy as np,pandas as pd,scanpy as sc
 from scipy.spatial.distance import cdist
 sample_df=sample_metadata_dict[sample_name]
 raw={row['sample_sub']:extract_raw_data_func(row) for _,row in sample_df.iterrows()}
 a=integrated_adata[integrated_adata.obs['sample'].str.contains(sample_name)].copy()
 if not a.obs_names.is_unique:raise ValueError('Raw-count joining requires unique case barcodes')
 mal=a[a.obs.annotation.isin(malignant_annotations)].copy()
 anchors=mal.obs_names[mal.obs.annotation.eq('Malignant_cell_border_mixed')]
 coords=mal.obsm['X_umap'];target=mal.obs_names.isin(anchors)
 if not target.any():raise ValueError('No qualifying anchors for signature analysis')
 distances=cdist(coords[target],coords).min(axis=0)
 nearby=mal.obs_names[distances<1.5]
 nearby_mal=nearby[mal.obs.loc[nearby,'annotation'].eq('Malignant_cell')]
 seeds=pd.Index(np.concatenate([anchors,nearby_mal]))
 trajectory=pd.read_csv(coord_csv_path)
 malignant=a.obs_names[a.obs.annotation.eq('Malignant_cell')]
 early=pd.Index(set(malignant)&set(trajectory['x']))
 a.obs['group']='Cells';a.obs.loc[early,'group']='malignant_pseudotime_early';a.obs.loc[seeds,'group']='metastatic_seeds'
 restored=restore_raw_counts(a,raw)
 if not np.array_equal(restored.obs_names,a.obs_names):raise ValueError('Raw-count row order does not match the selected spots')
 # Check each joined sample against its raw-count matrix.
 for sample,r in raw.items():
  idx=np.where(restored.obs['sample'].eq(sample))[0]
  r=r[:,~r.var_names.duplicated(keep='first')].copy();expected=r[restored.obs_names[idx],restored.var_names].X
  if hasattr(expected,'toarray'):expected=expected.toarray()
  if not np.array_equal(restored.X[idx],expected):raise ValueError('Joined counts do not match the sample raw-count matrix')
 grouped=restored[a.obs.group.isin(['malignant_pseudotime_early','metastatic_seeds']).values].copy()
 counts=grouped.copy()
 if not np.isfinite(counts.X).all() or (counts.X<0).any() or not np.array_equal(counts.X,np.rint(counts.X)):raise ValueError('Restored counts not nonnegative integers')
 sc.pp.normalize_total(grouped,target_sum=1e4);sc.pp.log1p(grouped);grouped.obs['group']=grouped.obs.group.astype('category')
 if grouped.raw is not None or len(grouped.layers):raise ValueError('Unexpected alternate DEG matrix')
 sc.tl.rank_genes_groups(grouped,groupby='group',method='wilcoxon')
 ranked=sc.get.rank_genes_groups_df(grouped,group='metastatic_seeds')
 ranked=ranked.sort_values(by='scores',ascending=False).copy()
 ranked.insert(0,'DEG_rank',range(1,len(ranked)+1))
 top=ranked.iloc[:top_n_genes].copy()
 lists=[stromal_marker_list,lymphocyte_marker_list,macrophage_marker_list,endothelial_marker_list,additional_remove_list]
 names=['stromal_genes','lymphocyte_genes','macrophage_genes','endothelial_genes','additional_remove_genes']
 union=set(sum([x or [] for x in lists],[]));top['matched_marker_lists']=[';'.join(n for n,ls in zip(names,lists) if g in (ls or [])) for g in top.names]
 top['excluded']=top.names.isin(union);final=top.loc[~top.excluded].copy()
 return {'case':a,'counts':counts,'log_deg':grouped,'rank':ranked,'top':top,'signature':final,'anchors':anchors,'seeds':seeds,'early':early,'min_distances':dict(zip(mal.obs_names,distances))}
