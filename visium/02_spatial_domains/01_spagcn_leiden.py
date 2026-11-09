import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import os
import random
import numpy as np
import argparse
from pathlib import Path

sc.set_figure_params(facecolor="white", figsize =(8,8))
sc.settings.verbosity =3

def Clustering_result(sample,method,nFeature_Spatial_min, nCount_Spatial_min,percent_mito_max, resolution,
                      sample_map, visium_root, coordinate_root):
    mapping = pd.read_csv(sample_map, sep="\t", dtype=str).fillna("")
    rows = mapping.loc[mapping["public_sample_id"].eq(sample)]
    if len(rows) != 1:
        raise ValueError(f"Expected one mapping row for public sample {sample}: {len(rows)}")
    row = rows.iloc[0]
    input_file_path = str(Path(visium_root) / row["space_ranger_sample"] / "Spatial_matrix")
    print(input_file_path)
    download_path = "h5ad_result/"
    if not os.path.exists(download_path):
        os.makedirs(download_path, exist_ok=True)
        print("Created directory:", download_path)
    else:
        print("Using directory:", download_path)
    adata = sc.read_visium(path = input_file_path, count_file = "filtered_feature_bc_matrix.h5")
    if row["tissue"] != "primary":
        table = pd.read_csv(Path(coordinate_root) / row["coordinate_table"])
        core = table.loc[table["sample_sub"].astype(str).eq(row["coordinate_sample_id"])]
        if len(core) != 1 or str(core.iloc[0]["sample"]) != row["space_ranger_sample"]:
            raise ValueError(f"Coordinate core mapping mismatch for {sample}")
        bounds = core.iloc[0]
        xy = adata.obsm["spatial"]
        keep = ((xy[:, 1] >= bounds["min_row"]) & (xy[:, 1] <= bounds["max_row"]) &
                (xy[:, 0] >= bounds["min_col"]) & (xy[:, 0] <= bounds["max_col"]))
        adata = adata[keep].copy()
        if adata.n_obs == 0:
            raise ValueError(f"No spots after mapped core crop for {sample}")
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    sc.pp.filter_cells(adata, min_genes=nFeature_Spatial_min)
    sc.pp.filter_genes(adata, min_cells=nCount_Spatial_min)
    adata = adata[adata.obs.total_counts < 50000, :]
    adata = adata[adata.obs.pct_counts_mt < percent_mito_max, :]
    if method == "SpaGCN" :
        import SpaGCN as spg
        import torch
        try:
            spatial = pd.read_csv(input_file_path + "/spatial/tissue_positions_list.csv",sep=",",header=None,na_filter=False,index_col=0)
            adata.obs["x1"]=spatial[1]
            adata.obs["x2"]=spatial[2]
            adata.obs["x3"]=spatial[3]
            adata.obs["x4"]=spatial[4]
            adata.obs["x5"]=spatial[5]
        except:
            spatial = pd.read_csv(input_file_path + "/spatial/tissue_positions.csv",sep=",",na_filter=False,index_col=0)
            adata.obs["x1"]=spatial["in_tissue"]
            adata.obs["x2"]=spatial["array_row"]
            adata.obs["x3"]=spatial["array_col"]
            adata.obs["x4"]=spatial["pxl_row_in_fullres"]
            adata.obs["x5"]=spatial["pxl_col_in_fullres"]

        adata.obs["x_array"]=adata.obs["x2"]
        adata.obs["y_array"]=adata.obs["x3"]
        adata.obs["x_pixel"]=adata.obs["x4"]
        adata.obs["y_pixel"]=adata.obs["x5"]
        adata=adata[adata.obs["x1"]==1]
        adata.var_names=[i.upper() for i in list(adata.var_names)]
        adata.var["genename"]=adata.var.index.astype("str")
        x_array=adata.obs["x_array"].tolist()
        y_array=adata.obs["y_array"].tolist()
        x_pixel=adata.obs["x_pixel"].tolist()
        y_pixel=adata.obs["y_pixel"].tolist()
        # Calculate the spatial adjacency matrix.
        s=1
        b=49
        adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, histology=False)
        adata.var_names_make_unique()
        spg.prefilter_genes(adata,min_cells=3) # avoiding all genes are zeros
        spg.prefilter_specialgenes(adata)
        #Normalize and take log for UMI
        sc.pp.normalize_per_cell(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=2000)
        adata = adata[:, adata.var.highly_variable]
        sc.pp.scale(adata)
        p=0.5
        #Find the l value given p
        l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)
        #Set seed
        r_seed=t_seed=n_seed=100
        clf=spg.SpaGCN()
        clf.set_l(l)
        #Set seed
        random.seed(r_seed)
        torch.manual_seed(t_seed)
        np.random.seed(n_seed)
        #Run
        clf.train(adata,adj,init_spa=True,init="louvain",res=resolution, tol=5e-3, lr=0.05, max_epochs=200)
        y_pred, prob=clf.predict()
        adata.obs["pred"]= y_pred
        adata.obs["pred"]=adata.obs["pred"].astype('category')
        #Do cluster refinement(optional)
        #shape="hexagon" for Visium data, "square" for ST data.
        adj_2d=spg.calculate_adj_matrix(x=x_array,y=y_array, histology=False)
        refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d, shape="hexagon")
        adata.obs["refined_pred"]=refined_pred
        adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')
        plot_color=[
        "#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C","#15821E","#3A84E6","#997273",
        "#787878","#DB4C6C","#9E7A7A","#554236","#AF5F3C","#93796C","#F9BD3F","#DAB370","#877F6C","#268785"
        ]
        #Plot refined spatial domains
        domains="refined_pred"
        num_celltype=len(adata.obs[domains].unique())
        adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
        sc.tl.pca(adata, svd_solver='arpack')
        sc.pp.neighbors(adata, n_neighbors=10)
    elif method == "Leiden" :
        #Normalize and take log for UMI
        sc.pp.normalize_per_cell(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=2000)
        sc.pl.highly_variable_genes(adata)
        adata = adata[:, adata.var.highly_variable]
        sc.pp.scale(adata)
        sc.tl.pca(adata, svd_solver='arpack')
        sc.pp.neighbors(adata, n_neighbors=9)
        sc.tl.leiden(adata, resolution=resolution)
        sc.pl.spatial(adata, color=["leiden"], size=1.5)
    safe_sample = "".join(ch if ch.isalnum() or ch in "._" else "_" for ch in sample)
    adata.obs["public_sample_id"] = sample
    adata.write_h5ad(filename= download_path + "/" + method + "_" + safe_sample + ".h5ad")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='SpaGCN or Leiden clustering of a mapped public Visium sample')
    parser.add_argument("--sample", type=str, required=True, help='Specify a Visium sample number')
    parser.add_argument("--method", choices=["SpaGCN", "Leiden"], help='Clustering method', default="SpaGCN")
    parser.add_argument("--nFeature_Spatial_min", type=int, required=False, help='Minimum number of features', default=200)
    parser.add_argument("--nCount_Spatial_min", type=int, required=False, help='Minimum number of counts', default=3)
    parser.add_argument("--percent_mito_max", type=int, required=False, help='Maximum mitochondrial count percentage', default=5)
    parser.add_argument("--resolution", type=float, required=False, help='Resolution for clustering', default=0.3)
    parser.add_argument("--sample-map", default="visium/resources/sample_id_mapping.tsv")
    parser.add_argument("--visium-root", default="data/visium/Visium")
    parser.add_argument("--coordinate-root", default="visium/resources/coordinates")
    args = parser.parse_args()

    Clustering_result(
        sample = args.sample,
        method = args.method,
        nFeature_Spatial_min = args.nFeature_Spatial_min,
        nCount_Spatial_min = args.nCount_Spatial_min,
        percent_mito_max = args.percent_mito_max,
        resolution = args.resolution,
        sample_map = args.sample_map,
        visium_root = args.visium_root,
        coordinate_root = args.coordinate_root,
    )
