import scanpy as sc
import pandas as pd
import numpy as np
import os
import argparse

# LuCA reference (tumor only, all stages) preprocessing for pySCENIC
# Output: raw UMI count matrix (cells x genes, gene symbols) and cell metadata as CSV

sc.settings.verbosity = 3
sc.logging.print_header()

def pySCENIC_preprocessing(input_file, output_dir, prefix, count_layer, gene_id_type="ensembl"):
    print(input_file)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print("Created directory:", output_dir)
    else:
        print("Using directory:", output_dir)

    adata = sc.read(input_file)
    print(adata)

    #Use raw UMI counts as the main matrix
    adata.X = adata.layers[count_layer].copy()

    if gene_id_type == "ensembl":
        #Map Ensembl gene IDs to gene symbols (MyGene.info)
        import mygene
        mg = mygene.MyGeneInfo()
        ensembl_ids = adata.var_names.tolist()
        out = mg.querymany(ensembl_ids, scopes="ensembl.gene", fields="symbol", species="human")
        out_df = pd.DataFrame(out)
        mapping = dict(zip(out_df["query"], out_df["symbol"]))
        adata.var["gene_symbol"] = [mapping.get(gene, gene) for gene in adata.var_names]
        adata.var_names = adata.var["gene_symbol"]

        #Remove genes without a gene symbol (NaN after mapping)
        print("Duplicated gene symbols before filtering:", adata.var_names.duplicated().sum())
        mask = ~adata.var["gene_symbol"].isna()
        adata = adata[:, mask].copy()
        print("Duplicated gene symbols after filtering:", adata.var_names.duplicated().sum())
        adata.var_names = adata.var_names.astype(str)
        adata.var_names_make_unique()
        print(adata)

    elif not adata.var_names.is_unique:
        raise ValueError("Prepared gene symbols must be unique")

    #Build cells x genes UMI count table
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    count_df = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)
    count_df.index.name = "cell_id"
    adata.obs.index.name = "cell_id"

    #Save UMI matrix and metadata
    count_df.to_csv(os.path.join(output_dir, prefix + "_UMI.csv"))
    adata.obs.to_csv(os.path.join(output_dir, prefix + "_metadata.csv"))
    print("Saved:", os.path.join(output_dir, prefix + "_UMI.csv"))
    print("Saved:", os.path.join(output_dir, prefix + "_metadata.csv"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Preprocessing of LuCA reference for pySCENIC')
    parser.add_argument("--input_file", type=str, required=False, help='Input h5ad file', default="../Data/LuCA_ref_only_tumor_all_stage.h5ad")
    parser.add_argument("--output_dir", type=str, required=False, help='Output directory', default="../Data/pySCENIC/RESOURCE/")
    parser.add_argument("--prefix", type=str, required=False, help='Prefix of output files', default="LuCA_Ref_Tumor")
    parser.add_argument("--count_layer", type=str, required=False, help='Layer containing raw UMI counts', default="count")
    parser.add_argument("--gene_id_type", choices=["ensembl", "symbol"], default="ensembl")
    args = parser.parse_args()

    pySCENIC_preprocessing(
        input_file = args.input_file,
        output_dir = args.output_dir,
        prefix = args.prefix,
        count_layer = args.count_layer,
        gene_id_type = args.gene_id_type
    )
