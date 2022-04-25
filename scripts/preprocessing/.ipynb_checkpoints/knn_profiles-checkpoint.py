import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

from sklearn.preprocessing import MinMaxScaler
from scipy.cluster.hierarchy import linkage
from sklearn.neighbors import NearestNeighbors

from python_imports import *

# Returns a clipped matrix given a list of max values (values to clip at) for each gene
def quantile_saturation(x = [], gene_quantiles = []):
    
    col_order = list(x.columns)
    
    x_clipped = pd.DataFrame()
    # Saturate the values based on 99% quantile (more robust than using the abs max)
    for i in x.index:
        x_clipped =  x_clipped.append(np.clip(x.loc[i], a_min=0, a_max = gene_quantiles))
        
    # Pandas automatically makes the new df columns alphabetically ordered
    # re-order columns how they were originally entered by the user
    
    x_clipped = x_clipped[col_order]
    
    return x_clipped

# based on pre-computed quantiles for each gene, we clip the data and apply min.max scaler 
# the range here is kept as in the main dataset 
def clip_and_scale(new_sample = [], gene_quantiles = [], scaler = []):
    # clip by quantiles 
    new_clip_sample = np.clip(new_sample, a_min=0, a_max = gene_quantiles)
    # min.max scaler -- same as for the integrated dataset 
    new_sample = new_clip_sample.values 
    new_sample = new_sample.reshape(1,-1)

    test_sample =scaler.transform(new_sample)
    return test_sample

def knn_mapping(
    raw_data, # the raw count data of dataset to do knn mapping to
    prof_anno=pd.DataFrame(), # the pathway class labels
    integrated_counts=integrated_counts, # the normalized, scaled integrated atlas
    pathway_genes=[], # pathway genes
    pathway_name="", # pathway name
    min_expr=0.2, # min. exp. cut-off
    min_genes_on=2, # min. genes exp.
    n_nearest=3, # no. of top __ profiles we want to map to each cluster using KNN
):
    
    # 3. quantiles from main dataset 
    x = integrated_counts[pathway_genes]
    gene_quantiles = np.quantile(x, q=0.99, axis=0)

    scaler= MinMaxScaler() 

    x_clipped = quantile_saturation(x, gene_quantiles)
    x_clipped.describe() 

    scaler.fit(x_clipped)

    df_clip = pd.DataFrame(scaler.transform(x_clipped), columns = pathway_genes)
    df_clip['cell_id'] = meta['cell_id'].values

    # Only keep the clusters in the integrated atlas that actually express the pathway
    df_clip["ON"] = (
        np.sum(df_clip[pathway_genes].values > min_expr, axis=1) > min_genes_on
    )

    # Isolate the gene counts
    pathway_mat = df_clip[df_clip["ON"]]
    pathway_mat.index = pathway_mat["cell_id"]

    # Add the profile annotations from R
    prof_anno = prof_anno
    pathway_mat.loc[:,"class_label"] = prof_anno.loc[:,"class_label"]

    # Average out the counts for each class
    pathway_mat = pathway_mat.groupby("class_label")[pathway_genes].mean()
    pathway_mat["class_label"] = pathway_mat.index

    pathway_mat["ON"] = "ON"

    # Initialize the nearest neighbors
    neigh = NearestNeighbors(n_neighbors=5, metric="euclidean", p=2, radius=0.4)
    neigh.fit(pathway_mat[pathway_genes])  # fit to all profiles
    
    data_counts=raw_data[pathway_genes]
    
    # DataFrames to save our knn mappings to
    df_save = pd.DataFrame(
        columns=["leiden"] + list(pathway_genes) + ["class_label", "ON"]
    )
    df_norm = pd.DataFrame(columns=list(pathway_genes))
    
    for target_cell_type in data_counts.index:
        # 2. which profile to map
        new_sample = data_counts[data_counts.index == target_cell_type][pathway_genes].mean()

        # 3.1 quantile-based saturation
        test_sample = clip_and_scale(new_sample, gene_quantiles, scaler)

        df_norm.loc[target_cell_type] = test_sample[0]
        # 4. find the top 5 matches to the target profile
        # we are using the min.max scaled data.
        nearest_n = neigh.kneighbors(test_sample, n_nearest, return_distance=False)

        nearest_one = pathway_mat.iloc[nearest_n[0]]

        nearest_one["leiden"] = target_cell_type

        df_save = pd.concat([df_save, nearest_one])

    return df_norm, df_save