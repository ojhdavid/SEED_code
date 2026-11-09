import os
import pickle
import subprocess
import pandas as pd
import scanpy as sc
from pyscenic.utils import load_motifs
from pyscenic.cli.utils import load_signatures
import argparse

# pySCENIC workflow for LuCA reference (tumor only, all stages)
# STEP 0: Build AnnData from UMI matrix and metadata
# STEP 1: Network inference (GRNBoost2)
# STEP 2-3: Regulon prediction (cisTarget)
# STEP 4: Cellular enrichment (AUCell)

sc.settings.verbosity = 3
sc.logging.print_header()

def run_command(cmd):
    #Run pySCENIC CLI command and stop if it fails
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

def pySCENIC_run(dataset_id, resources_dir, auxiliaries_dir, results_dir, grn_workers, ctx_workers, aucell_workers, seed):
    if not os.path.exists(results_dir):
        os.makedirs(results_dir, exist_ok=True)
        print("Created directory:", results_dir)
    else:
        print("Using directory:", results_dir)

    #Auxiliary data sets
    #TF list: https://github.com/aertslab/pySCENIC/tree/master/resources
    #Ranking databases and motif annotations: https://resources.aertslab.org/cistarget/
    human_tfs_fname = os.path.join(auxiliaries_dir, "allTFs_hg38.txt")
    ranking_dbs_fnames = [
        os.path.join(auxiliaries_dir, "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"),
        os.path.join(auxiliaries_dir, "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
    ]
    motif_annotations_fname = os.path.join(auxiliaries_dir, "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")

    #Input and output files
    metadata_fname = os.path.join(resources_dir, dataset_id + "_metadata.csv")
    umi_mtx_fname = os.path.join(resources_dir, dataset_id + "_UMI.csv")
    adjacencies_fname = os.path.join(results_dir, dataset_id + ".adjacencies.tsv")
    motifs_fname = os.path.join(results_dir, dataset_id + ".motifs.csv")
    regulons_dat_fname = os.path.join(results_dir, dataset_id + ".regulons.dat")
    aucell_mtx_fname = os.path.join(results_dir, dataset_id + ".auc.csv")
    anndata_fname = os.path.join(results_dir, dataset_id + ".h5ad")

    #STEP 0: Build AnnData from UMI matrix and metadata
    umi = pd.read_csv(umi_mtx_fname, index_col=0)
    metadata = pd.read_csv(metadata_fname, low_memory=False)
    obs = metadata[["cell_id", "sample", "dataset", "uicc_stage"]].set_index("cell_id").sort_index()
    adata = sc.AnnData(X=umi.sort_index())
    adata.obs = obs
    print(adata)
    adata.write_h5ad(anndata_fname)

    #STEP 1: Network inference based on GRNBoost2
    grn_cmd = ["pyscenic", "grn", umi_mtx_fname, human_tfs_fname,
               "-o", adjacencies_fname,
               "--num_workers", str(grn_workers)]
    if seed is not None:
        grn_cmd += ["--seed", str(seed)]
    run_command(grn_cmd)

    #STEP 2-3: Regulon prediction (cisTarget)
    run_command(["pyscenic", "ctx", adjacencies_fname] + ranking_dbs_fnames + [
        "--annotations_fname", motif_annotations_fname,
        "--expression_mtx_fname", umi_mtx_fname,
        "--output", motifs_fname,
        "--mask_dropouts",
        "--num_workers", str(ctx_workers)])

    df_motifs = load_motifs(motifs_fname)
    print(df_motifs.head())

    #Save regulons
    regulons = load_signatures(motifs_fname)
    print("Number of regulons:", len(regulons))
    with open(regulons_dat_fname, "wb") as f:
        pickle.dump(regulons, f)

    #STEP 4: Cellular enrichment (AUCell)
    run_command(["pyscenic", "aucell", umi_mtx_fname, motifs_fname,
                 "--output", aucell_mtx_fname,
                 "--num_workers", str(aucell_workers)])

    auc_mtx = pd.read_csv(aucell_mtx_fname, index_col=0)
    print("AUCell matrix:", auc_mtx.shape)
    print("Saved:", aucell_mtx_fname)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='pySCENIC for LuCA reference')
    parser.add_argument("--dataset_id", type=str, required=False, help='Dataset ID (prefix of input files)', default="LuCA_Ref_Tumor")
    parser.add_argument("--resources_dir", type=str, required=False, help='Directory of UMI matrix and metadata', default="../Data/pySCENIC/RESOURCE")
    parser.add_argument("--auxiliaries_dir", type=str, required=False, help='Directory of TF list, ranking databases and motif annotations', default="../Data/pySCENIC/AUXILLIARIES")
    parser.add_argument("--results_dir", type=str, required=False, help='Output directory', default="../Data/pySCENIC/RESULTS")
    parser.add_argument("--grn_workers", type=int, required=False, help='Number of workers for GRN step', default=12)
    parser.add_argument("--ctx_workers", type=int, required=False, help='Number of workers for cisTarget step', default=8)
    parser.add_argument("--aucell_workers", type=int, required=False, help='Number of workers for AUCell step', default=16)
    parser.add_argument("--seed", type=int, required=False, help='Random seed for GRNBoost2', default=None)
    args = parser.parse_args()

    pySCENIC_run(
        dataset_id = args.dataset_id,
        resources_dir = args.resources_dir,
        auxiliaries_dir = args.auxiliaries_dir,
        results_dir = args.results_dir,
        grn_workers = args.grn_workers,
        ctx_workers = args.ctx_workers,
        aucell_workers = args.aucell_workers,
        seed = args.seed
    )
