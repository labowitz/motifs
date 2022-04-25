import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler


# Import the atlas dataset of counts
counts = pd.read_csv('../../data/processed_data/' + "integrated_counts.csv", index_col = [0]).T
# Import metadata csv
meta = pd.read_csv('../../data/processed_data/' + "integrated_meta_data.csv", index_col = [0])
# Import reference gene list
all_pathways = pd.read_csv('../../data/raw_data/' + "allPathways_listGenes_dec2021.tsv", delimiter="\t")

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

# Pathway names
tgfb = all_pathways[all_pathways['pathway']=='Bmp_Tgfb']['gene'].values
notch = all_pathways[all_pathways['pathway']=='Notch']['gene'].values
eph_ephrin = np.concatenate((all_pathways[all_pathways['pathway']=='Eph_r']['gene'].values, all_pathways[all_pathways['pathway']=='Eph_l']['gene'].values))
wntr = all_pathways[all_pathways['pathway']=='Wnt']['gene'].values

