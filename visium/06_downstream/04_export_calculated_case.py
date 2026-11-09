#!/usr/bin/env python3
"""Export case assignments for plotting and scoring."""
import argparse, hashlib, json
from pathlib import Path

def array_hash(x):
 import numpy as np
 return hashlib.sha256(np.ascontiguousarray(x).tobytes()).hexdigest()

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--analysis-dir',required=True);p.add_argument('--output',required=True)
 p.add_argument('--sample-map',default=str(Path(__file__).resolve().parents[1]/'resources/sample_id_mapping.tsv'))
 args=p.parse_args()
 import numpy as np,pandas as pd,scanpy as sc
 src=Path(args.analysis_dir);dest=Path(args.output);dest.mkdir(parents=True,exist_ok=False)
 a=sc.read_h5ad(src/'integrated_cohort.h5ad');m=pd.read_csv(src/'SEED_early_membership.tsv',sep='\t').set_index('canonical_spot_id')
 if not a.obs_names.is_unique or not m.index.is_unique:raise ValueError('Duplicate canonical IDs')
 if not m.index.isin(a.obs_names).all():raise ValueError('Calculated assignments are absent from the computed cohort')
 a=a[a.obs_names.isin(m.index)].copy();m=m.loc[a.obs_names]
 mapping=pd.read_csv(args.sample_map,sep='\t',dtype=str).set_index('public_sample_id')
 if not mapping.index.is_unique:raise ValueError('Sample map is not one-to-one')
 for col in ['source_case_id','public_case_id','cohort','tissue']:
  a.obs[col]=a.obs.public_sample_id.astype(str).map(mapping[col])
  if a.obs[col].isna().any():raise ValueError('Incomplete mapping: '+col)
 if a.obs.public_case_id.nunique()!=1:raise ValueError('Expected the one calculated discovery case')
 before={k:array_hash(a.obsm[k]) for k in ['X_pca','X_umap','spatial']}
 if not np.array_equal(a.obs.annotation.astype(str).to_numpy(),m.annotation.astype(str).to_numpy()):raise ValueError('Assignment annotation differs from computed object')
 a.obs['case']=a.obs.public_case_id.astype(str)
 for old,new in [('is_seed','is_seed'),('is_anchor','is_ln_anchor'),('is_early','is_pseudotime_early')]:
  if not set(m[old].dropna().unique()).issubset({True,False}):raise ValueError('Expected boolean assignment '+old)
  a.obs[new]=m[old].to_numpy(dtype=bool)
 a.obs['is_primary_seed']=(a.obs.is_seed & a.obs.tissue.eq('primary')).to_numpy()
 a.obs['is_seed_anchor']=a.obs.is_ln_anchor.to_numpy()
 a.obs['is_malignant_for_seed']=a.obs.annotation.isin(['Malignant_cell','Malignant_cell_border_mixed']).to_numpy()
 a.obs['seed_radius']=1.5
 a.obs['seed_min_anchor_distance']=m.min_anchor_distance.to_numpy()
 a.obs['group']=m.group.to_numpy()  # SEED membership takes precedence over the early group for DEG analysis.
 a.obs['seed_group']=np.where(a.obs.is_seed,'SEED','Other malignant spots')
 seed=a.obs.is_seed.to_numpy();early=a.obs.is_pseudotime_early.to_numpy();scope=a.obs.is_malignant_for_seed.to_numpy()|seed|early
 labels=np.full(a.n_obs,'not_shown',dtype=object);labels[scope]='Malignant cells';labels[seed]='SEED';labels[early]='Malignant pseudotime early'
 a.obs['figure3_group']=labels;a.obs['figure3_display_scope']=scope
 a.obs['pseudotime_early_group']=np.where(early,'Malignant pseudotime early','Other malignant spots')
 def ids(mask):return hashlib.sha256('\n'.join(sorted(a.obs_names[mask].astype(str))).encode()).hexdigest()
 a.uns['seed_display_provenance']={'seed_id_sha256':ids(seed),'early_id_sha256':ids(early),'scope_id_sha256':ids(scope),'x_umap_sha256':before['X_umap'],'display_order':'default -> SEED -> early','membership_source':'SEED_early_membership.tsv'}
 after={k:array_hash(a.obsm[k]) for k in before}
 if before!=after:raise AssertionError('Metadata export changed coordinates')
 path=dest/'case_for_plots_and_scoring.h5ad';a.write_h5ad(path)
 a.obs.to_csv(dest/'assignments.tsv',sep='\t',index=False)
 record={'analysis_directory':str(src),'membership_sha256':hashlib.sha256((src/'SEED_early_membership.tsv').read_bytes()).hexdigest(),'coordinates_before':before,'coordinates_after':after,'seed_recomputed':False,'embedding_recomputed':False,'primary_seed':int((seed&a.obs.tissue.eq('primary')).sum()),'ln_anchors':int((a.obs.is_ln_anchor&a.obs.tissue.eq('lymph_node')).sum()),'nearby_ln':int((seed&~a.obs.is_ln_anchor&a.obs.tissue.eq('lymph_node')).sum()),'total_seed':int(seed.sum()),'early':int(early.sum()),'early_overlap':int((seed&early).sum()),'spots':a.n_obs}
 (dest/'metadata_export.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({k:record[k] for k in ['spots','primary_seed','ln_anchors','total_seed','early']}))
if __name__=='__main__':main()
