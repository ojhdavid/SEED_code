#!/usr/bin/env python3
"""Raw Space Ranger counts -> case embeddings -> SEED groups -> gene signature.
"""
import argparse,hashlib,importlib,importlib.metadata as md,importlib.util,json,os,site,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent

def module(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def runtime():
 from threadpoolctl import threadpool_info
 mods={}
 for name in ['scanpy','numpy','scipy','pandas','anndata','sklearn','umap','numba','bbknn','louvain','leidenalg']:
  m=importlib.import_module(name);mods[name]={'version':getattr(m,'__version__',None),'file':m.__file__}
 return {'executable':sys.executable,'python':sys.version,'modules':mods,'threadpools':threadpool_info(),'user_site_enabled':site.ENABLE_USER_SITE,'environment':{k:os.environ.get(k) for k in ['PYTHONNOUSERSITE','PYTHONHASHSEED','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']}}

def build_sample_groups(coordinate_table, sample_map):
 import pandas as pd
 table=pd.read_csv(coordinate_table);ids=['ALK51','ALK59','ALK101','ALK102','ALK103']
 primary=pd.DataFrame({'sample':ids,'ALK':ids,'sample_sub':[s+'P' for s in ids]})
 merged=pd.merge(table.T,primary.T,left_index=True,right_index=True,how='outer').T.reset_index(drop=True)
 merged=merged.loc[merged.ALK.notna()]
 mapping=pd.read_csv(sample_map,sep='\t',dtype=str).set_index('source_sample_id')
 if not mapping.index.is_unique:raise ValueError('Sample map is not one-to-one')
 included=mapping['include_in_concat'].str.lower().map({'true':True,'false':False})
 if included.isna().any():raise ValueError('include_in_concat must contain only true/false')
 unknown=sorted(set(merged.sample_sub.astype(str))-set(mapping.index))
 if unknown:raise ValueError('Samples absent from mapping: '+', '.join(unknown))
 merged=merged.loc[merged.sample_sub.isin(mapping.index[included])]
 groups={alk:df.sort_values(by=['sample_sub']) for alk,df in merged.groupby('ALK')[merged.columns[1:]]}
 return groups,mapping

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for k in ['visium-root','coordinate-table','annotation-root','trajectory-csv','output']:p.add_argument('--'+k,required=True)
 p.add_argument('--sample-map',default=str(HERE.parent/'resources/sample_id_mapping.tsv'));p.add_argument('--signature-config',default=str(HERE/'resources/signature_config.json'));p.add_argument('--save-intermediates',action='store_true');args=p.parse_args()
 import numpy as np,pandas as pd,scanpy as sc,anndata as ad
 w=module('preprocessing_ops',HERE/'preprocessing_ops.py');h=module('signature_utils',HERE/'signature_utils.py');derivation=module('derive_signature',HERE/'01_derive_signature.py');c=h.load_signature_config(args.signature_config);out=Path(args.output);out.mkdir(parents=True,exist_ok=False);inter=out/'intermediates';inter.mkdir();t0=time.monotonic();before=runtime()
 np.random.seed(42)
 groups,mapping=build_sample_groups(args.coordinate_table,args.sample_map)
 all_samples={};all_raw={};stages=[];copyback=[];ctx={'sample':None};input_hashes={str(Path(q).resolve()):h.sha(q) for q in [args.coordinate_table,args.trajectory_csv,args.sample_map,args.signature_config]}
 exclusion_path=Path(args.signature_config).resolve().parent/c['exclusions_file'];input_hashes[str(exclusion_path)]=h.sha(exclusion_path)
 def qualified(x):
  x.obs['source_obs_name']=x.obs_names.astype(str);x.obs['source_sample_id']=x.obs['sample'].astype(str);x.obs['public_sample_id']=x.obs['source_sample_id'].map(mapping.public_sample_id)
  if x.obs['public_sample_id'].isna().any():raise ValueError('Incomplete sample mapping')
  x.obs['canonical_spot_id']=x.obs.public_sample_id.astype(str)+'::'+x.obs.source_obs_name;x.obs_names=x.obs.canonical_spot_id.to_numpy();x.obs['sample']=x.obs.public_sample_id.astype(str)
  if not x.obs_names.is_unique:raise ValueError('Nonunique canonical spot IDs')
 def state():
  r=np.random.get_state();return {'algorithm':r[0],'state':r[1].tolist(),'position':r[2],'has_gauss':r[3],'cached_gaussian':r[4]}
 orig_pca=sc.tl.pca
 def pca(x,*a,**kw):
  z=np.asarray(x.X);record={'sample':ctx['sample'],'PCA_arguments':kw,'input':h.arrhash(z),'spot_order':x.obs_names.tolist(),'gene_order':x.var_names.tolist(),'is_view':x.is_view,'C_contiguous':bool(z.flags.c_contiguous),'F_contiguous':bool(z.flags.f_contiguous),'strides':list(z.strides),'global_rng_before':state()}
  if args.save_intermediates and ctx['sample'].startswith('ALK103'):x.copy().write_h5ad(inter/(ctx['sample']+'_before_PCA.h5ad'))
  result=orig_pca(x,*a,**kw);record['PCA_output']=h.arrhash(x.obsm['X_pca']);record['global_rng_after']=state();stages.append(record)
  return result
 sc.tl.pca=pca
 for case,df in groups.items():
  rows=df[~df.sample_sub.astype(str).str.contains('T2',regex=False)];items={};print('case',case,flush=True)
  for _,row in rows.iterrows():
   raw_path=Path(args.visium_root)/str(row['sample'])/'Spatial_matrix';raw=sc.read_visium(raw_path);bounds=[row[k] for k in ['min_row','max_row','min_col','max_col']]
   if not any(pd.isna(v) for v in bounds):raw=w.crop(raw,*bounds)
   all_raw[str(row.sample_sub)]=ad.AnnData(X=raw.X.copy(),obs=raw.obs.copy(),var=raw.var.copy())
   for file in raw_path.rglob('*'):
    if file.is_file():input_hashes[str(file.resolve())]=h.sha(file)
  for _,row in rows.iterrows():
   ctx['sample']=str(row.sample_sub);items[ctx['sample']]=w.make_adata(row,Path(args.visium_root),Path(args.annotation_root),20000)[0];key=ctx['sample'][:-1] if ctx['sample'].endswith('P') else ctx['sample'];ann=Path(args.annotation_root)/('Annotation_professor_'+key+'.csv');input_hashes[str(ann.resolve())]=h.sha(ann)
  items=w.adata_intersection(items);case_obj=w.integration(items);bbknn_hash=h.arrhash(case_obj.obsp['connectivities'].data);sc.pp.neighbors(case_obj,n_neighbors=8,random_state=42);sc.tl.leiden(case_obj,resolution=0.3,random_state=42);sc.tl.umap(case_obj,random_state=42)
  if case=='ALK103' and args.save_intermediates:case_obj.copy().write_h5ad(inter/'ALK103_case_embedding.h5ad')
  for name,x in items.items():
   mask=case_obj.obs['sample'].astype(str).eq(name).to_numpy();x.obsm['X_umap']=np.asarray(case_obj.obsm['X_umap'])[mask];all_samples[name]=x;copyback.append({'sample':name,'UMAP_hash':h.arrhash(x.obsm['X_umap']),'case_BBKNN_data_hash':bbknn_hash})
 sc.tl.pca=orig_pca
 cohort=w.integration(w.adata_intersection(all_samples))
 for sample in cohort.obs['sample'].unique():
  if '_' not in str(sample):continue
  mask=cohort.obs['sample'].eq(sample);labels=cohort.obs.loc[mask,'annotation'];update=labels.str.contains('Malignant',case=False,na=False)&labels.ne('Malignant_cell_mixed_with_stromal');cohort.obs.loc[mask,'annotation']=labels.where(~update,'Malignant_cell_border_mixed')
 cohort.obs['annotation']=cohort.obs.annotation.replace({'Malignant_cell(Hot)':'Malignant_cell','Malignant_cell(Cold)':'Malignant_cell'})
 result=derivation.process_alk_sample_and_extract_meta_markers(sample_name=c['source_sample_name'],integrated_adata=cohort,sample_metadata_dict=groups,extract_raw_data_func=lambda row:all_raw[str(row.sample_sub)].copy(),malignant_annotations=c['malignant_annotations'],coord_csv_path=args.trajectory_csv,top_n_genes=c['top_n_genes'],stromal_marker_list=c['marker_lists']['stromal_genes'],lymphocyte_marker_list=c['marker_lists']['lymphocyte_genes'],macrophage_marker_list=c['marker_lists']['macrophage_genes'],endothelial_marker_list=c['marker_lists']['endothelial_genes'],additional_remove_list=c['marker_lists']['additional_remove_genes'])
 case=result['case'];bars=case.obs_names.astype(str);source=case.obs['sample'].astype(str);seed=bars.isin(result['seeds']);early=bars.isin(result['early']);anchor=bars.isin(result['anchors']);primary_mask=source.str.endswith('P').to_numpy();qids=source.map(mapping.public_sample_id).to_numpy()+'::'+bars.to_numpy();membership=pd.DataFrame({'canonical_spot_id':qids,'source_sample_id':source.to_numpy(),'source_barcode':bars,'annotation':case.obs.annotation.astype(str).to_numpy(),'is_seed':seed,'is_anchor':anchor,'is_early':early,'group':case.obs.group.to_numpy(),'min_anchor_distance':[result['min_distances'].get(b,np.nan) for b in bars]});membership.to_csv(out/'SEED_early_membership.tsv',sep='\t',index=False)
 counts=result['counts'];g=result['log_deg'];rawhash=h.arrhash(counts.X);loghash=h.arrhash(g.X)
 for obj in [case,counts,g,cohort]:qualified(obj)
 cohort.write_h5ad(out/'integrated_cohort.h5ad');counts.write_h5ad(out/'deg_counts.h5ad');g.write_h5ad(out/'deg_log_ranked.h5ad');counts.obs.to_csv(out/'DEG_group_spots.tsv',sep='\t',index=False);pd.DataFrame({'gene':counts.var_names}).to_csv(out/'DEG_gene_universe.tsv',sep='\t',index=False)
 for key,name in [('rank','DEG_sorted.tsv'),('top','top260.tsv'),('signature','derived_signature.tsv')]:result[key].rename(columns={'names':'gene'}).to_csv(out/name,sep='\t',index=False)
 result['top'].loc[result['top'].excluded].rename(columns={'names':'gene'}).to_csv(out/'excluded_markers.tsv',sep='\t',index=False)
 full=[]
 for group in g.obs.group.cat.categories:
  frame=sc.get.rank_genes_groups_df(g,group=group);frame.insert(0,'group',group);full.append(frame)
 pd.concat(full).to_csv(out/'DEG_all_groups.tsv',sep='\t',index=False);np.savez_compressed(out/'DEG_native_arrays.npz',**{k:g.uns['rank_genes_groups'][k] for k in ['names','scores','logfoldchanges','pvals','pvals_adj']})
 h.dump(out/'pca_summary.json',stages);h.dump(out/'embedding_transfer_summary.json',copyback);h.dump(out/'input_hashes.json',input_hashes)
 status={'scope':'cohort','embedding_computed':True,'case_order':list(groups),'sample_order':list(all_samples),'cohort_shape':list(cohort.shape),'primary_seed':int((seed&primary_mask).sum()),'ln_anchors':int((anchor&~primary_mask).sum()),'nearby_ln':int((seed&~primary_mask&~anchor).sum()),'total_seed':int(seed.sum()),'early':int(early.sum()),'early_overlap':int((seed&early).sum()),'groups':g.obs.group.value_counts().to_dict(),'DEG_genes':g.n_vars,'raw_counts':rawhash,'log_normalized':loghash,'top_n':len(result['top']),'derived_count':len(result['signature']),'rank_params':g.uns['rank_genes_groups']['params'],'runtime_before':before,'runtime_after':runtime(),'seconds':time.monotonic()-t0};h.dump(out/'run_summary.json',status);print(json.dumps({k:status[k] for k in ['primary_seed','ln_anchors','total_seed','early','DEG_genes','derived_count','seconds']}),flush=True)
if __name__=='__main__':main()
