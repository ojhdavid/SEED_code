#!/usr/bin/env python3
"""Run SMC-103 signature derivation, spatial/UMAP plotting and scoring."""
import argparse,json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',required=True);a=p.parse_args()
 config_path=Path(a.config).resolve();cfg=json.loads(config_path.read_text())
 expected={'visium_root','coordinate_table','annotation_root','trajectory_csv','output'}
 if set(cfg)!=expected:raise ValueError('Config keys must be exactly: '+', '.join(sorted(expected)))
 def path(v):
  q=Path(v).expanduser();return q.resolve() if q.is_absolute() else (config_path.parent/q).resolve()
 cfg={k:path(v) for k,v in cfg.items()};out=cfg['output']
 if out.exists():raise FileExistsError('Refusing existing output: '+str(out))
 for k,v in cfg.items():
  if k!='output' and not v.exists():raise FileNotFoundError(str(v))
 out.mkdir(parents=True);logs=out/'logs';logs.mkdir()
 py=os.environ.get('PYTHON',sys.executable)
 env=os.environ.copy();env.update(PYTHON=py,PYTHONNOUSERSITE='1',PYTHONDONTWRITEBYTECODE='1',PYTHONHASHSEED='0',OMP_NUM_THREADS='64',OPENBLAS_NUM_THREADS='64',MKL_NUM_THREADS='64')
 env.setdefault('MPLCONFIGDIR',str(out/'cache/matplotlib'));env.setdefault('NUMBA_CACHE_DIR',str(out/'cache/numba'))
 commands=[]
 def run(name,argv):
  argv=list(map(str,argv));t=time.monotonic()
  with (logs/(name+'.log')).open('w')as handle:r=subprocess.run(argv,cwd=ROOT,env=env,stdout=handle,stderr=subprocess.STDOUT)
  commands.append({'step':name,'argv':argv,'returncode':r.returncode,'seconds':time.monotonic()-t});(out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
  if r.returncode:raise RuntimeError(name+' failed; see '+str(logs/(name+'.log')))
 analysis=out/'analysis';export=out/'case';cohort=cfg['visium_root']/'ALK103/Spatial_matrix'
 run('raw_to_signature',['bash',ROOT/'visium/06_downstream/run_raw_signature.sh','--visium-root',cfg['visium_root'],'--coordinate-table',cfg['coordinate_table'],'--annotation-root',cfg['annotation_root'],'--trajectory-csv',cfg['trajectory_csv'],'--output',analysis,'--save-intermediates'])
 run('export_case',[py,'-s',ROOT/'visium/06_downstream/04_export_calculated_case.py','--analysis-dir',analysis,'--output',export])
 obj=export/'case_for_plots_and_scoring.h5ad'
 run('primary_spatial',[py,'-s',ROOT/'visium/05_seed_definition/05_plot_primary_spatial.py','--h5ad',obj,'--visium-dir',cohort,'--output',out/'primary_spatial','--label','SMC-103 SEED'])
 run('seed_umap',[py,'-s',ROOT/'visium/05_seed_definition/03_plot_case_seed_umap.py','--input',obj,'--case-id','LUAD-SMC-103','--output-dir',out/'seed_umap','--formats','png,pdf'])
 run('early_umap',[py,'-s',ROOT/'visium/05_seed_definition/02_plot_seed_umap.py','--input',obj,'--output-dir',out/'early_umap','--formats','png,pdf'])
 run('derived_signature_scoring',[py,'-s',ROOT/'visium/06_downstream/02_plot_metascore.py','--input',obj,'--signature',analysis/'derived_signature.tsv','--score-mode','integrated-x','--output-dir',out/'scoring'])
 print(json.dumps({'output':str(out),'steps_completed':len(commands),'signature':str(analysis/'derived_signature.tsv'),'primary_spatial':str(out/'primary_spatial/primary_full.png')},indent=2))
if __name__=='__main__':main()
