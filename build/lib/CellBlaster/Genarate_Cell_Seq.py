import os
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import argparse
import pandas as pd
import numpy as np
import scanpy as sc
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import DBSCAN
from sklearn.neighbors import LocalOutlierFactor
import warnings

# Suppress unnecessary warnings for cleaner output
warnings.filterwarnings('ignore')

def parse_args():
    """
    Handles user input for file paths, identifiers, and filtering criteria.
        
    Returns:
        argparse.Namespace: An object containing the following attributes:
            - input (str): Path to the single-cell h5ad file.
            - symbol (str): Prefix used for naming all generated output files.
            - outdir (str): Destination directory for the results.
            - filter_keywords (list): Keywords to identify and remove specific genes.
                                    Genes whose names contain one or more strings you provide will be excluded.
    """

    parser = argparse.ArgumentParser(description="Single-cell analysis pipeline: Outlier removal, DEG identification, and Percentile Symbolic Encoding.")
    # Required inputs
    parser.add_argument("-i", "--input", required=True, help="Full path to the input .h5ad file (e.g., /path/to/data.h5ad)")
    parser.add_argument("-s", "--symbol", required=True, help="Output symbol/prefix (e.g., SRP309176)")
    parser.add_argument("-o", "--outdir", required=True, help="Directory to save the output files")
    parser.add_argument("-f", "--filter_keywords", nargs='+', default=['LNC'],
                         help="One or more keywords to filter out genes (e.g., AthLNC mt)")
    return parser.parse_args()


def main():
    args = parse_args()
    # Ensure output directory exists
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
        print(f"01.Read data and creat directory: {args.outdir}")

    # Define output file paths based on requirements
    all_txt_path = os.path.join(args.outdir, f"{args.symbol}.all.txt")
    celltype_txt_path = os.path.join(args.outdir, f"{args.symbol}.Celltype.txt")
    degs_csv_path = os.path.join(args.outdir, f"{args.symbol}.topDEGs.csv")


    # Step 1: Data Loading & Initial Preprocessing
    adata = sc.read_h5ad(args.input) 

    # Filter Long Non-Coding RNAs (starting with 'LNC')
    adata.var["name"] = adata.var.index.tolist()
    filter_pattern = '|'.join(args.filter_keywords)
    print(f"Filtering genes containing: {filter_pattern}")
    filtered_genes = adata.var["name"].str.contains(filter_pattern, case=False, na=False)
    adata = adata[:, ~filtered_genes]

    # Filter low-expression genes (min 5 cells)
    sc.pp.filter_genes(adata, min_cells=5)
    adata.uns['raw'] = adata.X.copy()

    # Store raw counts and perform normalization
    sc.pp.normalize_total(adata,target_sum=10000)
    sc.pp.log1p(adata)


    # Step 2: Outlier Detection (DBSCAN + LOF)
    print("02.Performing Outlier Detection...")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata.raw = adata.copy()
    adata = adata[:, adata.var['highly_variable']]
    sc.pp.scale(adata)
    sc.tl.pca(adata, svd_solver='arpack')

    # Loop through each cell type to find local outliers
    eps = 10 
    min_samples = 10
    outliers = []
    lof_n_neighbors = 20
    lof_contamination = 0.05
    for cell_type in adata.obs['Celltype'].unique():
        ct_mask = adata.obs['Celltype'] == cell_type
        ct_data = adata[ct_mask]
        X_scaled = StandardScaler().fit_transform(ct_data.obsm['X_pca'])
        # Method 1: DBSCAN
        db_labels  = DBSCAN(eps=eps, min_samples=min_samples).fit_predict(X_scaled)
        db_outliers = db_labels  == -1
        # Method 2: Local Outlier Factor
        lof_labels = LocalOutlierFactor(n_neighbors=lof_n_neighbors, contamination=lof_contamination).fit_predict(X_scaled)
        lof_outliers = lof_labels == -1
        # Combine (Outlier if either method flags it)
        combined_outliers = db_outliers | lof_outliers
        ct_indices = np.where(ct_mask)[0]
        outliers.extend(ct_indices[combined_outliers])

    # Mark and filter outliers
    outlier_mask = np.zeros(adata.n_obs, dtype=bool)
    outlier_mask[outliers] = True
    adata.obs['outlier'] = 'Inlier'
    adata.obs.loc[outlier_mask, 'outlier'] = 'Outlier'

    # Keep only Inliers and revert to raw counts for DEG analysis
    adata = adata.raw.to_adata()
    adata.X = adata.uns['raw'].copy()
    adata = adata[ adata.obs['outlier'] == 'Inlier'].copy()
    adata.raw = adata.copy()
    # Re-normalize cleaned data
    sc.pp.normalize_total(adata,target_sum=10000)
    sc.pp.log1p(adata)


    # Step 3: DEG Analysis (Wilcoxon Rank-Sum Test)
    print("03.DEG Analysis (Wilcoxon Rank-Sum Test)...")
    sc.tl.rank_genes_groups(adata, groupby='Celltype',  method='wilcoxon', key='rank_genes_groups')
    degs = sc.get.rank_genes_groups_df(adata, group = None, pval_cutoff=0.05, log2fc_min=0.5) 
    # Select top 150 DEGs per group
    top200_degs = degs.groupby('group').head(150) 
    top200_degs.to_csv(degs_csv_path,sep='\t')

    # Filter adata to include only these unique top DEGs
    gene_3 = top200_degs['names'].unique().tolist()
    adata = adata.raw.to_adata()
    adata1 = adata[:, gene_3].copy()
    sc.pp.normalize_total(adata1,target_sum=10000)
    sc.pp.log1p(adata1)
    sc.pp.scale(adata1)

    # Step 4: Symbolic Percentile Encoding & Saving Results
    print("04.Generating symbolic outputs...")
    # Convert expression to rank-based percentiles
    df = pd.DataFrame(adata1.X,index = adata1.obs.index,columns = adata1.var.index)
    df = df.rank(pct=True, axis=1) * 100
    df = df.round().astype(int)
    # Define encoding symbols (A-Y, a-y)
    symbols = [chr(ord('A') + i) for i in range(25)] + [chr(ord('a') + i) for i in range(25)]
    bins = np.linspace(0, 100, 51) 
    # Bin percentiles into 50 categories and map to symbols
    df_symbolic = df.apply(lambda col: pd.cut(col, bins=bins, labels=symbols, include_lowest=True))
    # Save Celltype mapping
    df_symbolic.to_csv(all_txt_path, sep='\t')
    Celltype = adata1.obs[['Celltype']]
    Celltype['Cell'] = df.index
    Celltype.to_csv(celltype_txt_path, sep='\t', index=False)


if __name__ == "__main__":
    main()