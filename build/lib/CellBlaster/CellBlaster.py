import os
import sys
import argparse
import requests
import pandas as pd
import numpy as np
import seaborn as sns 
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from glob import glob
from tqdm import tqdm
import subprocess

def parse_args():
    parser = argparse.ArgumentParser(description="Process and merge multiple species sequence datasets.")
    # -t / --dabase_type: Determines the working directory for reference data.
    #Choices are strictly limited to 'Dicot' (双子叶植物) or 'Monocot' (单子叶植物).
    parser.add_argument("-t", "--dabase_type", required=True, choices=['Dicot', 'Monocot'], 
                        help="Query database of Dicot or Monocot.")
    
    # -s / --symbols: A list of identifiers for datasets already processed and stored in the CellBlaster.
    # Accepts multiple values separated by spaces. User can choose the database by yourself.
    # Example: -s SRP406470 SRP390780 (The script will look for .all.txt and .Celltype.txt for these).
    parser.add_argument("-s", "--symbols", nargs='+', required=True, 
                        help="One or more dataset symbols (e.g., SRP406470 SRP390780).")
    
    # -o / --output_path: The directory where final results (merged files and OG matrix) will be saved.
    # Defaults to the current working directory.    
    parser.add_argument("-o", "--output_path", default="./", help="Output files.")

    # -q / --query: The full file path to the user's specific single-cell data in .h5ad format.
    # This is the dataset that will be compared against the database.
    parser.add_argument("-q", "--query", required=True, 
                        help="Full path to the input .h5ad file (e.g., /path/to/data.h5ad)")
    
    # -qs / --query_symbol: A unique name/prefix for the query dataset.
    # Used for naming query-specific output files (e.g., <query_symbol>_topDEGs.csv, <query_symbol>.all.txt).
    parser.add_argument("-qs", "--query_symbol", required=True, 
                        help="Query symbol (e.g., SRP309176)")
    
    # -f / --filter_keywords: Keywords used to exclude unwanted genes (e.g., non-coding, organelle genes).
    # Any gene names containing these strings (case-insensitive) will be dropped during preprocessing.
    # Default is 'LNC' (Long Non-Coding RNA).
    parser.add_argument("-f", "--filter_keywords", nargs='+', default=['LNC'],
                        help="One or more keywords to filter out genes (e.g., AthLNC mt)")
    return parser.parse_args()

#-------------------------------Download databasen data-------------------------
def download_dataset(url, dest_dir, filename):
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)
    target_path = os.path.join(dest_dir, filename)
    if os.path.exists(target_path):
        print(f"File {filename} already exists. Skipping...")
        return target_path
    print(f"Downloading {filename} from {url}...")
    try:
        with requests.get(url, stream=True, timeout=60) as r:
            r.raise_for_status()
            with open(target_path, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        print(f"Successfully downloaded: {filename}")
        return target_path
    except Exception as e:
        print(f"Failed to download {filename}. Error: {e}")
        return None
    
args = parse_args()
data_dir = args.output_path
data_dir = data_dir if data_dir.endswith('/') else data_dir + '/'
data_dir = data_dir  + "01.DataBase"
data_dir = data_dir if data_dir.endswith('/') else data_dir + '/'
os.makedirs(data_dir, exist_ok=True)
os.chdir(data_dir)

for symbol in args.symbols:
    print(f"Downloading {symbol} files...")
    all_url= f"https://zenodo.org/records/19347241/files/{symbol}.all.txt?download=1"
    download_dataset(all_url, data_dir, f"{symbol}.all.txt")
    celltype_url= f"https://zenodo.org/records/19347241/files/{symbol}.Celltype.txt?download=1"
    download_dataset(celltype_url, data_dir, f"{symbol}.Celltype.txt")
    deg_url= f"https://zenodo.org/records/19347241/files/{symbol}.topDEGs.csv?download=1"
    download_dataset(deg_url, data_dir, f"{symbol}.topDEGs.csv")


base_url = "https://zenodo.org/records/19347241/files/"
download_dataset(f"{base_url}Celltype.txt?download=1", data_dir, "Celltype.txt")
if args.dabase_type == "Dicot":
    Ortho_url = f"{base_url}Dicot_Orthogroups.txt?download=1"
    download_dataset(Ortho_url, data_dir, "Dicot_Orthogroups.txt")
    Orthogroups = f"{data_dir}Dicot_Orthogroups.txt"
elif args.dabase_type == "Monocot":
    Ortho_url = f"{base_url}Monocot_Orthogroups.txt?download=1"
    download_dataset(Ortho_url, data_dir, "Monocot_Orthogroups.txt")
    Orthogroups = f"{data_dir}Monocot_Orthogroups.txt"
else:
    print("Error: Please enter 'Dicot' or 'Monocot'.")
    sys.exit(1)

print(f"Starting to process {len(args.symbols)} datasets from {args.dabase_type}...")


#------------------------Read the comparison results of the respective homologous genes------------------------
char_to_int = {char: idx for idx, char in enumerate(list("ABCDEFGHIJKLMNOPQRSTUVWXYabcdefghijklmnopqrstuvwxy")[:50])}
records = []
with open(Orthogroups, 'r') as f:
    for line in f:
        if ':' in line:
            og, genes_str = line.strip().split(':', 1)
            genes = genes_str.strip().split()
            for gene in genes:
                records.append((og, gene))

OG = pd.DataFrame(records, columns=['OG', 'Gene'])
gene_to_og = dict(zip(OG['Gene'], OG['OG']))

def map_and_group_by_og(df, gene_to_og):
    symbols = [chr(ord('A') + i) for i in range(25)] + [chr(ord('a') + i) for i in range(25)]
    symbol_to_num = {s: i for i, s in enumerate(symbols)}
    num_to_symbol = {i: s for i, s in enumerate(symbols)}
    df = df.loc[:, df.columns.isin(gene_to_og)]
    df.columns = [gene_to_og[gene] for gene in df.columns]
    df_num = df.map(symbol_to_num.get)
    df_median = df_num.T.groupby(level=0).median().T
    df_symbol = df_median.round().astype(int).map(num_to_symbol.get)
    return df_symbol

def encode_df(df, colnames):
    df = df[colnames].apply(lambda col: col.map(lambda x: char_to_int.get(x, -1)))
    return df.to_numpy()

def cal_similarity(query, dataset):
    common_col = sorted(set(query.columns) & set(dataset.columns))
    print(f"Number of common OGs: {len(common_col)}")
    query_encoded = encode_df(query, common_col)
    dataset_encoded = encode_df(dataset, common_col)
    similarity_matrix = np.zeros((query_encoded.shape[0], dataset_encoded.shape[0]))
    for i in tqdm(range(query_encoded.shape[0]), desc="Computing similarity"):
        diff = np.abs(dataset_encoded - query_encoded[i])
        similarity_matrix[i] = diff.sum(axis=1)
    similarity = pd.DataFrame(similarity_matrix, index=query.index, columns=dataset.index)
    similarity = similarity / len(common_col)
    return similarity

def draw_cell_count(data,filename):
    custom_cmap = LinearSegmentedColormap.from_list("custom_blue_gradient",["#ffffff", "#5b9ff1","#9c1019"])
    plt.figure(figsize=(12,12)) 
    ax = sns.heatmap(data, cmap=custom_cmap, annot=True, fmt="d",linewidths=2,
                linecolor='#c9d6df', square=True,cbar_kws={'shrink': 0.8}, 
                annot_kws={"size": 12,"color": "black"} )
    plt.title("Top1 Cell count Between Cell Types", fontsize=16, weight='bold')
    plt.xlabel("Reference", fontsize=18, weight='bold')
    plt.ylabel("Query", fontsize=18, weight='bold')
    plt.xticks(fontsize=14, rotation=45, ha='right', weight='bold')
    plt.yticks(fontsize=14, rotation=0, weight='bold')
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label("Cell count percent", fontsize=18)
    plt.tight_layout()
    for y_index, row in enumerate(data.values):
        max_col_index = row.argmax()
        rect = Rectangle((max_col_index, y_index), 1, 1, fill=False, edgecolor="#fd1158", linewidth=3)
        ax.add_patch(rect)
    plt.savefig(filename)

def draw_cell_percent(data,filename):
    custom_cmap = LinearSegmentedColormap.from_list("custom_blue_gradient",["#ffffff", "#5b9ff1","#9c1019"])
    plt.figure(figsize=(12,12)) 
    ax = sns.heatmap(data, cmap=custom_cmap, annot=True, fmt=".1f",linewidths=2,
                linecolor='#c9d6df', square=True,cbar_kws={'shrink': 0.8}, 
                annot_kws={"size": 12,"color": "black"} )
    plt.title("Prediction of Cell Types (%)", fontsize=16, weight='bold')
    plt.xlabel("Reference", fontsize=18, weight='bold')
    plt.ylabel("Query", fontsize=18, weight='bold')
    plt.xticks(fontsize=14, rotation=45, ha='right', weight='bold')
    plt.yticks(fontsize=14, rotation=0, weight='bold')
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label("Cell count percent (%)", fontsize=18)
    plt.tight_layout()
    for y_index, row in enumerate(data.values):
        max_col_index = row.argmax()
        rect = Rectangle((max_col_index, y_index), 1, 1, fill=False, edgecolor="#fd1158", linewidth=3)
        ax.add_patch(rect)
    plt.savefig(filename)

def draw_cell_percent_clustering(data, filename):
    custom_cmap = LinearSegmentedColormap.from_list("custom_blue_gradient", ["#ffffff", "#5b9ff1", "#9c1019"])
    n_rows, n_cols = data.shape
    cell_size = 0.5 
    figsize = (cell_size * n_cols + 4, cell_size * n_rows + 4)  
    g = sns.clustermap(
        data,
        cmap=custom_cmap,
        annot=True,
        fmt=".1f",
        linewidths=2,
        linecolor="#c9d6df",
        cbar_kws={"shrink": 0.7},
        annot_kws={"size": 12, "color": "black"},
        figsize=figsize,
        dendrogram_ratio=(0.1, 0.1),
        cbar_pos=(1.1, 0.3, 0.03, 0.4), 
        tree_kws={'linewidths': 2, 'colors': '#524748'}  
    )
    g.ax_heatmap.set_title("Prediction of Cell Types (%)",fontsize=18,weight="bold",pad=60)
    g.ax_heatmap.set_xlabel("Reference", fontsize=18, weight="bold", labelpad=10)
    g.ax_heatmap.set_ylabel("Query", fontsize=18, weight="bold", labelpad=10)
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=13)
    plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=13)
    cbar = g.ax_cbar
    cbar.set_ylabel("Cell count percent (%)", fontsize=16, labelpad=15)
    cbar.tick_params(labelsize=14)
    reordered_data = data.iloc[g.dendrogram_row.reordered_ind, g.dendrogram_col.reordered_ind]
    for y_index, row in enumerate(reordered_data.values):
        max_col_index = row.argmax()
        rect = Rectangle((max_col_index, y_index), 1, 1, fill=False, edgecolor="#fd1158", linewidth=3)
        g.ax_heatmap.add_patch(rect)
    plt.tight_layout()
    plt.savefig(filename, bbox_inches="tight", dpi=300)
    plt.close()

def  draw_ctype_bar(df,filename):
    n_rows, n_cols = df.shape
    custom_colors = [
        "#469cd8", "#ff7f0e", "#a7ff83", "#FF00FF", "#fa4659", "#9467bd",
        "#8c564b", "#ff847b", "#33a02c", "#fc5c9c", "#f0eca4", "#8A2BE2" ,
        "#a0ecf6", "#7f7f7f", "#00FFFF", "#46ecd8", "#eea291", "#f5c7f7",
        "#069df5", "#e7b37f", "#2f9296", "#de95ba", "#7FFF00","#fffb00"]
    colors = custom_colors[:n_cols]
    color_dict = dict(zip(df.columns, colors))
    fig, ax = plt.subplots(figsize=(max(8, n_rows * 0.8), 6))
    bar_width = 0.7
    indices = np.arange(n_rows)
    for i, idx in enumerate(df.index):
        row = df.loc[idx]
        sorted_items = row.sort_values(ascending=True)
        cum_height = 0
        for celltype, val in sorted_items.items():
            ax.bar(i, val, bottom=cum_height, width=bar_width, color=color_dict[celltype], edgecolor='black', linewidth=0.4)
            cum_height += val
    for spine in ax.spines.values():
        spine.set_linewidth(1.8)
    ax.tick_params(axis='both', width=1.6, length=6, labelsize=12)
    ax.set_xticks(indices)
    ax.set_xticklabels(df.index, rotation=45, ha='right', fontsize=12)
    ax.set_ylabel("Proportion (%)", fontsize=14)
    ax.set_xlabel("Query Cell Type", fontsize=14)
    ax.set_title("Celltype Assignment (%)", fontsize=16, fontweight='bold')
    handles = [plt.Rectangle((0,0),1,1, color=color_dict[ct]) for ct in df.columns]
    ax.legend(handles, df.columns, title="Predicted Cell Type", bbox_to_anchor=(1.02,1), loc='upper left', fontsize=12, title_fontsize=14)
    sns.despine()
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight') 


#---------------------Loop through each data set in the retrieval database that you have set up.------------------
datasets = []
Celltype_data = []
all_DEG = []
for i, symbol in enumerate(args.symbols):
    prefix = f"d{i+1}"
    all_file = f"{data_dir}{symbol}.all.txt"
    celltype_file = f"{data_dir}{symbol}.Celltype.txt"
    deg_file = f"{data_dir}{symbol}.topDEGs.csv"
    if not os.path.exists(all_file) or not os.path.exists(celltype_file):
        print(f"Warning: Files for {symbol} not found. Skipping...")
        continue
    print(f"Processing {symbol} with prefix {prefix}...")
    ds = pd.read_csv(all_file, sep='\t', index_col=0)
    ct = pd.read_csv(celltype_file, sep='\t')
    deg = pd.read_csv(deg_file, sep='\t',index_col=0)
    ds.index = [f"{prefix}_{str(idx)}" for idx in ds.index]
    try:
        ds = map_and_group_by_og(ds, gene_to_og)
    except NameError:
        print("Error: 'map_and_group_by_og' or 'gene_to_og' not defined. Please ensure they are defined in your script.")
        sys.exit(1)
    ct['Cell'] = [f"{prefix}_{str(idx)}" for idx in ct['Cell']]
    datasets.append(ds)
    Celltype_data.append(ct)
    all_DEG.append(deg)

Celltype_data = pd.concat(Celltype_data)
print(f"Merging complete! Total cells: {len(datasets)}")


# --- ----Call the Generate_Cell_Seq script to generate the dataset of the query data specified by the user--------
input_f = args.query
symbol = args.query_symbol
out_dir = args.output_path
out_dir = out_dir if out_dir.endswith('/') else out_dir + '/'
out_data = out_dir +"02.QueryData/"
filter_list = args.filter_keywords
current_script_dir = os.path.dirname(os.path.abspath(__file__))
code_path = os.path.join(current_script_dir, "Genarate_Cell_Seq.py")
os.makedirs(out_data, exist_ok=True)
os.chdir(out_data)
cmd = [
    "python", code_path,
    "-i", input_f,
    "-s", symbol,
    "-o", out_data,
    "-f"
] + filter_list
try:
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    print("Genarate Cell Seq Output:\n", result.stdout)
except subprocess.CalledProcessError as e:
    print("Error occurred while running the script:")
    print(e.stderr)

# -------------------Modify output path----------------------------------
blast_out = out_dir + "03.Blast_Result/"  
os.makedirs(blast_out, exist_ok=True)
os.chdir(blast_out)

blast_out_new = blast_out + symbol + "_blast_1"
os.makedirs(blast_out_new, exist_ok=True)
os.chdir(blast_out_new)

all_DEG = pd.concat(all_DEG)
all_DEG['names'] = all_DEG['names'].map(gene_to_og)
df_query = pd.read_csv(out_data + symbol+ ".topDEGs.csv", sep='\t',index_col=0)
df_query['names'] = df_query['names'].map(gene_to_og)
OMG_all = pd.DataFrame(0, index=df_query['group'].unique(), columns=all_DEG['group'].unique())
for g1 in df_query['group'].unique():
    ogs1 = {x for x in df_query[df_query['group'] == g1]['names'] if pd.notna(x)}  # 过滤 NaN
    for g2 in all_DEG['group'].unique():
        ogs2 = {x for x in all_DEG[all_DEG['group'] == g2]['names'] if pd.notna(x)}  # 过滤 NaN
        shared = ogs1 & ogs2
        OMG_all.loc[g1, g2] = len(shared)

OMG_all = OMG_all.sort_index().sort_index(axis=1)
OMG_all.to_csv("OMG_all_"+symbol+".csv")


###########################################Read the query data#######################################
query = pd.read_csv(out_data + symbol+ ".all.txt", sep='\t',index_col=0)
Celltype_query = pd.read_csv(out_data + symbol+ ".Celltype.txt", sep='\t')
query = map_and_group_by_og(query, gene_to_og)

################################Calculate the correlation for each dataset separately and select the cell with the minimum loss.###########################################
datasets = [ds[sorted(set(query.columns) & set(ds.columns))] for ds in datasets]
similarity_df = pd.concat([cal_similarity(query, ds) for ds in datasets])
similarity_df = similarity_df.fillna(similarity_df.max().max())
result_df = similarity_df.idxmin(axis=1).reset_index(name='TopColumn')
result_df.columns = ['query_cell', 'dataset_cell']
result_df["query_celltype"] = result_df["query_cell"].map(dict(zip(Celltype_query["Cell"], Celltype_query["Celltype"])))
result_df["dataset_celltype"] = result_df["dataset_cell"].map(dict(zip(Celltype_data["Cell"], Celltype_data["Celltype"])))
result_df.to_csv('01.Cell_match_'+symbol+'.csv')

##############################The number of each cell type detected in the top1 cells#######################################
pivot_df = result_df.pivot_table(index='query_celltype',columns='dataset_celltype',aggfunc='size',fill_value=0)
filename = '01.Ctype_count_'+symbol+'.png'
draw_cell_count(pivot_df,filename)
pivot_normalized = pivot_df.div(pivot_df.sum(axis=1), axis=0)*100
pivot_normalized = pivot_normalized.round(2)
filename = '02.Ctype_percent_'+symbol+'.png'
draw_cell_percent(pivot_normalized,filename)

##################################Calculate the quantity of OGs.########################################################
OMG = OMG_all[OMG_all.columns.intersection(pivot_df.columns)]
filename = '03.OMG_count_'+symbol+'.png'
draw_cell_count(OMG,filename)
###百分比
multiplied_df = OMG * pivot_normalized
multiplied_df = multiplied_df.div(multiplied_df.sum(axis=1), axis=0)*100
filename = '04.OMG_similar_'+symbol+'.png'
draw_cell_percent(multiplied_df,filename)

##################################Start multiple rounds of correction#######################################################
qualified_rows = multiplied_df[(multiplied_df.max(axis=1) - multiplied_df.apply(lambda x: x.nlargest(2).iloc[1], axis=1)) >= 30].index.tolist()
print("Qualified celltype: ",qualified_rows)  #预测正确的细胞类型
unqualified_rows = multiplied_df[~multiplied_df.index.isin(qualified_rows)].index.tolist()
df_quality = pd.DataFrame({'Celltype': qualified_rows + unqualified_rows,'Passed_QC': ['Yes'] * len(qualified_rows) + ['No'] * len(unqualified_rows)})
df_quality.to_csv('02.QC_Ctype_'+symbol+'.csv')
OMG_similar = multiplied_df.loc[multiplied_df.index.isin(qualified_rows)]
OMG_similar.to_csv('03.OMG_similar_'+symbol+'.csv')

i = 2
while  unqualified_rows:
    print(i," round of blast celltype")
    ##########################################Set the working path###########################################
    out_new = blast_out +"/"+ symbol + "_blast_" +str(i)
    os.makedirs(out_new, exist_ok=True)
    os.chdir(out_new)
    #########################Extract the cell types that need to be subjected to another blast process.##################################
    Celltype_query = Celltype_query[~Celltype_query['Celltype'].isin(qualified_rows)]
    query = query.loc[query.index.isin(Celltype_query['Cell'])]
    Celltype_data = Celltype_data[~Celltype_data['Celltype'].isin(qualified_rows)]
    datasets = [ds.loc[ds.index.isin(Celltype_data['Cell'])] for ds in datasets]
    ###########################Calculate the correlation with each dataset separately and obtain the most similar cells.#########################################
    similarity_df = pd.concat([cal_similarity(query, ds) for ds in datasets])
    similarity_df = similarity_df.fillna(similarity_df.max().max())
    result_df = similarity_df.idxmin(axis=1).reset_index(name='TopColumn')
    result_df.columns = ['query_cell', 'dataset_cell']
    result_df["query_celltype"] = result_df["query_cell"].map(dict(zip(Celltype_query["Cell"], Celltype_query["Celltype"])))
    result_df["dataset_celltype"] = result_df["dataset_cell"].map(dict(zip(Celltype_data["Cell"], Celltype_data["Celltype"])))
    result_df.to_csv('01.Cell_match_'+symbol+'.csv')
    ###############################The number of each cell type detected in the top1 cells#########################################
    pivot_df = result_df.pivot_table(index='query_celltype',columns='dataset_celltype',aggfunc='size',fill_value=0)
    filename = '01.Ctype_count_'+symbol+'.png'
    draw_cell_count(pivot_df,filename)
    pivot_normalized = pivot_df.div(pivot_df.sum(axis=1), axis=0)*100
    pivot_normalized = pivot_normalized.round(2)
    filename = '02.Ctype_percent_'+symbol+'.png'
    draw_cell_percent(pivot_normalized,filename)
    ###################################Calculate the quantity of OGs.#########################################################
    OMG = OMG_all[OMG_all.columns.intersection(pivot_df.columns)]
    OMG_temp = OMG.loc[:, ~OMG.columns.isin(qualified_rows)]
    OMG_temp = OMG_temp.loc[~OMG_temp.index.isin(qualified_rows)]
    filename = '03.OMG_count_'+symbol+'.png'
    draw_cell_count(OMG_temp,filename)
    multiplied_df = OMG_temp * pivot_normalized
    multiplied_df = multiplied_df.div(multiplied_df.sum(axis=1), axis=0)*100
    filename = '04.OMG_similar_'+symbol+'.png'
    draw_cell_percent(multiplied_df,filename)
    ########################################Determine whether further correction is necessary#############################################
    diff = multiplied_df.max(axis=1) - multiplied_df.apply(lambda x: x.nlargest(2).iloc[1], axis=1)
    qualified_new = diff[diff >= 30].index.tolist()
    if not qualified_new:
        qualified_new = [diff.idxmax()]
    print("Qualified celltype: ",qualified_new)
    qualified_rows = qualified_new + qualified_rows
    unqualified_rows = multiplied_df[~multiplied_df.index.isin(qualified_rows)].index.tolist()
    df_quality = pd.DataFrame({'Celltype': qualified_rows + unqualified_rows,'Passed_QC': ['Yes'] * len(qualified_rows) + ['No'] * len(unqualified_rows)})
    df_quality.to_csv('02.QC_Ctype_'+symbol+'.csv')
    #######################################Store the final result###################################################
    OMG_similar = multiplied_df.loc[multiplied_df.index.isin(qualified_new)]
    OMG_similar.to_csv('03.OMG_similar_'+symbol+'.csv')
    i = i+1

#######################################Visualization of the final result##########################
os.chdir(blast_out)
file_list = glob(os.path.join(out_new, '**', '03.OMG_similar_*.csv'), recursive=True)
combined_df = pd.concat([pd.read_csv(f, index_col=0) for f in file_list], axis=0).fillna(0) # 合并所有文件并填充空值为 0
final = "Final_Visiual"
os.makedirs(final, exist_ok=True)
os.chdir(blast_out +"/"+final)
combined_df.to_csv('01.OMG_similar_all.csv')
filename = '01.OMG_similar.png'
draw_cell_percent(combined_df,filename)
filename = '02.Ctype_bar.png'
draw_ctype_bar(combined_df,filename)

max_colnames = combined_df.idxmax(axis=1)
result = pd.DataFrame({'Query_CellType': combined_df.index,'Assigned_CellType': max_colnames.values})
ct_matrix = pd.read_csv( f"{data_dir}Celltype.txt",sep='\t')
def normalize(s):
    return s.strip()

ct_pairs = {frozenset(map(normalize, pair)) for pair in ct_matrix[['ct1', 'ct2']].values}
result['Matched'] = result.apply(lambda r: normalize(r['Query_CellType']) == normalize(r['Assigned_CellType']) or
frozenset([normalize(r['Query_CellType']), normalize(r['Assigned_CellType'])]) in ct_pairs,axis=1)
result['accuracy'] = result['Matched'].mean()
result.to_csv('02.Accuracy.csv')

combined_df = pd.read_csv('01.OMG_similar_all.csv', index_col=0)
filename = '01.OMG_similar_cluster.png'
draw_cell_percent_clustering(combined_df,filename)


