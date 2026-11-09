#!/usr/bin/env python
"""Render primary SEED membership on the Space Ranger tissue image."""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd,scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Patch
from matplotlib.collections import PatchCollection

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--h5ad',required=True);p.add_argument('--visium-dir',required=True);p.add_argument('--output',required=True);p.add_argument('--label',required=True);p.add_argument('--sample',default='LUAD-SMC-103P');args=p.parse_args()
 out=Path(args.output);out.mkdir(parents=True,exist_ok=False);a=sc.read_h5ad(args.h5ad);raw=sc.read_visium(args.visium_dir)
 obs=a.obs.loc[a.obs['sample'].astype(str).eq(args.sample)&a.obs.is_malignant_for_seed.astype(bool)].copy()
 if not len(obs):raise ValueError('No eligible primary spots for --sample')
 bar=obs.source_obs_name.astype(str);idx=raw.obs_names.get_indexer(bar)
 if (idx<0).any() or bar.duplicated().any() or not raw.obs_names.is_unique:raise ValueError('Spatial barcode join has missing or duplicate entries')
 stored=a[obs.index].obsm['spatial'];xyraw=raw.obsm['spatial'][idx]
 if not np.array_equal(stored,xyraw):raise ValueError('Stored and Space Ranger spatial coordinates differ')
 lib=next(iter(raw.uns['spatial'].values()));img=lib['images']['hires'];sf=float(lib['scalefactors']['tissue_hires_scalef']);xy=xyraw*sf;rad=float(lib['scalefactors']['spot_diameter_fullres'])*sf/2*.65;seed=obs.is_seed.astype(bool).to_numpy();H,W=img.shape[:2]
 obs['fullres_x']=xyraw[:,0];obs['fullres_y']=xyraw[:,1];obs['hires_x']=xy[:,0];obs['hires_y']=xy[:,1];obs.to_csv(out/'spatial_join.tsv',sep='\t',index=False)
 for zoom in [False,True]:
  fig,ax=plt.subplots(figsize=(10,9));ax.imshow(img,origin='upper')
  for mask,color,scale,alpha in [(~seed,'#ffe082',1,.45),(seed,'#d7191c',3,1)]:
   ax.add_collection(PatchCollection([Circle(tuple(v),rad*scale) for v in xy[mask]],facecolor=color,edgecolor='white' if scale>1 else 'none',linewidth=.3,alpha=alpha))
  if zoom:ax.set_xlim(xy[:,0].min()-30,xy[:,0].max()+30);ax.set_ylim(xy[:,1].max()+30,xy[:,1].min()-30)
  else:ax.set_xlim(0,W);ax.set_ylim(H,0)
  ax.axis('off');ax.set_aspect('equal');ax.set_title(args.label+f'\nPrimary SEED (n={seed.sum()})',fontsize=12);ax.legend(handles=[Patch(facecolor='#d7191c',label=f'Primary SEED ({seed.sum()})'),Patch(facecolor='#ffe082',label='Other eligible malignant')],loc='upper center',bbox_to_anchor=(.5,-.02),frameon=False);fig.text(.5,.01,'SEED markers shown at 3× radius.',ha='center',fontsize=9)
  for ext in ['png','pdf']:fig.savefig(out/(('primary_zoom' if zoom else 'primary_full')+'.'+ext),dpi=200,bbox_inches='tight')
  plt.close(fig)
 (out/'spatial_summary.json').write_text(json.dumps(dict(missing=0,duplicates=0,stored_spatial_exact=True,primary_seed=int(seed.sum()),early_overlap=int((obs.is_seed.astype(bool)&obs.is_pseudotime_early.astype(bool)).sum()),scalefactor=sf,recomputed_embedding=False),indent=2))
if __name__=='__main__':main()
