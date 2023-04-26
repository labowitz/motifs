import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler


# Import the atlas dataset of counts
counts = pd.read_csv('../../data/processed_data/' + "integrated_counts.csv", index_col = [0]).T
# Import metadata csv
meta = pd.read_csv('../../data/processed_data/' + "integrated_meta_data.csv", index_col = [0])
# Import reference gene list
all_pathways = pd.read_csv('../../data/raw_data/pathbank/' + "pathway_df.csv", index_col=0)

# Returns a clipped matrix given a list of max values (values to clip at) for each gene
def quantile_saturation(x = [], integrated_gene_quantiles = []):
    
    col_order = list(x.columns)
    
    x_clipped = pd.DataFrame()
    # Saturate the values based on 99% quantile (more robust than using the abs max)
    for i in x.index:
        x_clipped =  x_clipped.append(np.clip(x.loc[i], a_min=0, a_max = integrated_gene_quantiles))
        
    # Pandas automatically makes the new df columns alphabetically ordered
    # re-order columns how they were originally entered by the user
    
    x_clipped = x_clipped[col_order]
    
    return x_clipped


## Normalize the dataset

# Normalize the atlas counts to have 1e4 per cluster
integrated_counts = np.log1p(counts.div(counts.sum(axis=1),axis=0)*1e4)
integrated_counts["cell_id"] = meta["cell_id"]
integrated_counts["Cell_class"] = meta["Cell_class"]

# Pathway names
tgfb = all_pathways[all_pathways['pathway']=='Bmp_Tgfb']['gene'].values
notch = all_pathways[all_pathways['pathway']=='Notch']['gene'].values
eph_ephrin = np.concatenate((all_pathways[all_pathways['pathway']=='Eph_r']['gene'].values, all_pathways[all_pathways['pathway']=='Eph_l']['gene'].values))
wntr = all_pathways[all_pathways['pathway']=='Wnt']['gene'].values

def minmax_scaled_quantile(
    integrated_counts=integrated_counts, # the normalized, scaled integrated atlas
    pathway_genes=[], # pathway genes
    pathway_name="", # pathway name
    min_expr=0.2, # min. exp. cut-off
    min_genes_on=2, # min. genes exp.
):
    x = integrated_counts[pathway_genes]

    gene_quantiles = np.quantile(x, q=0.99, axis=0)

    x_clipped = quantile_saturation(x, gene_quantiles)
    x_clipped.describe() 
    
    scaler= MinMaxScaler()

    scaler.fit(x_clipped)

    df_clip = pd.DataFrame(scaler.transform(x_clipped), columns = pathway_genes)
    df_clip['cell_id'] = meta['cell_id'].values

    # Only keep the clusters in the integrated atlas that actually express the pathway
    df_clip["ON"] = (
        np.sum(df_clip[pathway_genes].values > min_expr, axis=1) > min_genes_on
    )
    
    return df_clip