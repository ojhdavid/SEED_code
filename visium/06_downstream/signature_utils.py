"""Run metadata and signature configuration helpers."""
import csv
import hashlib
import json
from pathlib import Path

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for x in iter(lambda:f.read(8*1024*1024),b''):h.update(x)
 return h.hexdigest()

def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,default=str)+'\n')

def arrhash(x):
 import numpy as np
 from scipy import sparse
 a=x.toarray() if sparse.issparse(x) else np.asarray(x)
 return {'shape':list(a.shape),'dtype':str(a.dtype),'sha256':hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()}

def load_signature_config(path):
    """Read the discovery settings and ordered marker categories from one TSV."""
    path = Path(path).resolve()
    config = json.loads(path.read_text())
    required = ["stromal_genes", "lymphocyte_genes", "macrophage_genes",
                "endothelial_genes", "additional_remove_genes"]
    if config["top_n_genes"] != 260:
        raise ValueError("The discovery configuration requires top_n_genes=260")
    exclusions = path.parent / config["exclusions_file"]
    marker_lists = {category: [] for category in required}
    with exclusions.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["category", "gene"]:
            raise ValueError("Exclusion TSV must contain category and gene columns")
        for row in reader:
            category, gene = row["category"], row["gene"]
            if category not in marker_lists or not gene or gene.strip() != gene:
                raise ValueError("Invalid exclusion category or gene")
            marker_lists[category].append(gene)
    if any(not genes for genes in marker_lists.values()):
        raise ValueError("Each exclusion category must contain genes")
    config["marker_lists"] = marker_lists
    return config
