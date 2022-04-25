#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib import cm, colors
import scipy.cluster.hierarchy as sch
import seaborn as sns
import datetime
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.cluster import KMeans
from itertools import combinations
from scipy.spatial.distance import cosine, squareform, euclidean
import math

# Here are the lists of genes for our main pathways of interest.
datadir = '../../data/raw_data/'
all_pathways = pd.read_csv(datadir + "allPathways_listGenes_dec2021.tsv", delimiter="\t")

bmpr = all_pathways[all_pathways['pathway']=='Bmp_Tgfb']['gene'].values
notch = all_pathways[all_pathways['pathway']=='Notch']['gene'].values
eph_ephrin = np.concatenate((all_pathways[all_pathways['pathway']=='Eph_r']['gene'].values, all_pathways[all_pathways['pathway']=='Eph_l']['gene'].values))
wntr = all_pathways[all_pathways['pathway']=='Wnt']['gene'].values

def get_genes(adata, genes):
    '''This function gets genes of interest that have not been filtered out.
    
    Input: 
    
    adata: AnnData object whose genes have been filtered
    genes: list of gene names
    
    Output:
    
    list of genes of interest that are actually in the filtered dataset
    '''
    list_genes = []
    for i in genes:
        if True in adata.var_names.str.endswith(i):
            list_genes.append(i)
    return list_genes

def vis_pre_processing(adata, genes_range = (0,10000), counts_range = (0, 400000),
                      genes_threshold = 2000, counts_threshold=2000, title=""):
    '''A histogram of genes/cell and counts/cell, a boxplot of 15 highest
    expressed genes, a scatterplot of genes against counts. This is to visualize 
    the data before we filter genes and cells.
    
    Inputs:
   
    adata: AnnData object
    genes_range: the range of values we want to display on our histogram for genes/cell
    counts_range: the range of values we want to display on our histogram for counts/cell
    genes_threshold: where to show a red line (threshold) for number of genes, for easy visualization
    counts_threshold: where to show a red line (threshold) for counts, for easy visualization
    title: title of figure, if you want to save
    
    Outputs:
    
    fig: figure with genes/cell and counts/cell, boxplot of highest expressed genes, and scatterplot of
    genes vs. counts
    '''
    fig, ax = plt.subplots(2, 2, figsize = (9, 9))
    ax[0,0].hist(adata.obs['n_genes_per_cell'][:], bins = 100, range = genes_range)
    ax[0,0].axvline(x=genes_threshold, color='r', linestyle='dashed', linewidth=2)
    ax[0,0].grid(False)
    ax[0,0].set_title('Histogram of Number of Genes per Cell')
    ax[0,0].set_xlabel('Number of Genes')
    ax[0,0].set_ylabel('Frequency (# of Cells)')
    
    ax[0,1].hist(adata.obs['n_total_counts_per_cell'][:], bins = 100, range=counts_range)
    ax[0,1].axvline(x=counts_threshold, color='r', linestyle='dashed', linewidth=2)
    ax[0,1].grid(False)
    ax[0,1].set_title('Histogram of Counts per Cell')
    ax[0,1].set_xlabel('Log Counts per Cell')
    ax[0,1].set_ylabel('Frequency')
    ax[0,1].set_xscale('log')
    
    sc.pl.highest_expr_genes(adata, n_top=15, ax=ax[1,0], show=False)
    ax[1,0].set_title('15 Highest Expressed Genes')
    sc.pl.scatter(adata, x='n_total_counts_per_cell', y='n_genes', size = 20, 
                  title = 'Counts vs. Genes', ax=ax[1,1], show=False)
    ax[1,1].set_xlabel('Counts')
    ax[1,1].set_ylabel('Genes')
    ax[1,1].grid(False)
    fig.suptitle(title)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    return fig

def filter_data(adata, min_counts=2000, min_genes=2000, min_cells=3):
    '''This function filters out cells with a minimum value of genes/cell and counts/cell. It also filters
    out genes that are expressed in at least a minimum number of cells.
    
    Inputs:
    
    adata: AnnData object
    min_counts: minimum counts/cell filter
    min_genes: minimum genes/cell filter
    min_cells: minimum number of cells gene should be expressed in 
    
    Output:
    
    adata: filtered AnnData object
    '''
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    return adata

def vis_post_processing(adata, genes_range = (0,10000), counts_range = (0, 400000),title="",
                       genes_threshold = 2000, counts_threshold=2000):
    '''Histograms of genes and total counts, and finally a scatter plot
    of genes against counts.
    
    Inputs:
   
    adata: AnnData object
    genes_range: the range of values we want to display on our histogram for genes/cell
    counts_range: the range of values we want to display on our histogram for counts/cell
    genes_threshold: where to show a red line (threshold) for number of genes, for easy visualization
    counts_threshold: where to show a red line (threshold) for counts, for easy visualization
    title: title of figure, if you want to save
    
    Outputs:
    
    fig: figure with genes/cell and counts/cell, boxplot of highest expressed genes, and scatterplot of
    genes vs. counts.
    '''
    
    fig,ax = plt.subplots(2, 2, figsize = (8,8))
    ax[0,0].hist(adata.obs['n_genes_per_cell'][:], bins = 100, range = genes_range)
    ax[0,0].axvline(x=genes_threshold, color='r', linestyle='dashed', linewidth=2)
    ax[0,0].grid(False)
    ax[0,0].set_title('Histogram of Number of Genes per Cell')
    ax[0,0].set_xlabel('Number of Genes')
    ax[0,0].set_ylabel('Frequency (# of Cells)')
    ax[0,1].hist(adata.obs['n_total_counts_per_cell'][:], bins = 100, range = counts_range)
    ax[0,1].axvline(x=counts_threshold, color='r', linestyle='dashed', linewidth=2)
    ax[0,1].grid(False)
    ax[0,1].set_title('Histogram of Counts per Cell')
    ax[0,1].set_xlabel('Log Counts per Cell')
    ax[0,1].set_ylabel('Frequency (# of Cells)')
    ax[0,1].set_xscale('log')
    fig.suptitle(title)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    sc.pl.scatter(adata, x='n_counts', y='n_genes', size=20,ax=ax[1,0], title='Counts vs. Genes',show=False)
    ax[1,0].set_xlabel('Counts')
    ax[1,0].set_ylabel('Genes')
    ax[1,0].grid(False)
    ax[1,1].axis('off')
    plt.show()
    return fig

def normalize_data(adata,count,batch_key = None,**kwargs):
    '''This function normalizes the data to a total number of counts, does a log(x+1) transformation, and sets a 
    raw attribute of the anndata object. It also sets the highly-variable genes attribute of the anndata 
    observation parameters.
    
    Inputs:
    adata: AnnData object
    count: total number of counts to normalize with.
    **kwargs: any other arguments to normalize the total with (applied to sc.pp.normalize_total fxn)
    
    Outputs:
    adata: AnnData object with normalized and log(x+1) transformed counts, with a .raw attribute that holds the 
    normalized counts, and a highly_variable_genes attribute that store the top highly variable genes
    '''
    sc.pp.normalize_total(adata, target_sum=count, **kwargs)
    sc.pp.log1p(adata)
    adata.raw=adata
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key=batch_key)
    return adata

def merge_genes(adata, all_genes):
    '''This function "densifies" the anndata object. It preserves highly variable genes and all genes of
    interest.
    
    Input: 
    
    adata: the AnnData object
    all_genes: list of lists of genes of interest
    
    Output:
    
    adata: AnnData object with only highly-variable genes and all genes of interest.
    '''
    joint_genes=adata.var.highly_variable
    for i in all_genes:
        for j in i:
            joint_genes = joint_genes|adata.var_names.str.startswith(j)
    adata=adata[:,joint_genes]
    return adata

def scale_data(adata):
    '''This function regresses out the AnnData object againist total counts per cell, and scales the 
    gene expression matrix so that each gene has zero mean and unit variance.
    
    Input:
    adata: AnnData object with ['n_total_counts_per_cell'] parameter in observations
    
    Output:
    AnnData object
    '''
    sc.pp.regress_out(adata, ['n_total_counts_per_cell'])
    sc.pp.scale(adata)
    return adata

# Delete this function here
def marker_gene_expression(anndata, marker_dict, gene_symbol_key=None, partition_key='leiden'):
    """
    A function to get mean z-score expressions of marker genes
     
    Inputs:
    
        anndata         - An AnnData object containing the data set and a partition
        marker_dict     - A dictionary with cell-type markers. The markers should be stores as anndata.var_names
        or an anndata.var field with the key given by the gene_symbol_key input
        gene_symbol_key - The key for the anndata.var field with gene IDs or names that correspond to the marker 
                          genes
        partition_key   - The key for the anndata.obs field where the cluster IDs are stored. The default is
                          'leiden' 
    """

    #Test inputs
    if partition_key not in anndata.obs.columns.values:
        print('KeyError: The partition key was not found in the passed AnnData object.')
        print('   Have you done the clustering? If so, please tell pass the cluster IDs with the AnnData object!')
        raise

    if (gene_symbol_key != None) and (gene_symbol_key not in anndata.var.columns.values):
        print('KeyError: The provided gene symbol key was not found in the passed AnnData object.')
        print('   Check that your cell type markers are given in a format that your anndata object knows!')
        raise
        
    if gene_symbol_key:
        gene_ids = anndata.var[gene_symbol_key]
    else:
        gene_ids = anndata.var_names

    clusters = anndata.obs[partition_key].cat.categories
    n_clust = len(clusters)
    marker_exp = pd.DataFrame(columns=clusters)
    marker_exp['cell_type'] = pd.Series({}, dtype='str')
    marker_names = []
    
    z_scores = sc.pp.scale(anndata, copy=True)

    i = 0
    for group in marker_dict:
        # Find the corresponding columns and get their mean expression in the cluster
        for gene in marker_dict[group]:
            ens_idx = np.in1d(gene_ids, gene) #Note there may be multiple mappings
            if np.sum(ens_idx) == 0:
                continue
            else:
                z_scores.obs[ens_idx[0]] = z_scores.X[:,ens_idx].mean(1) #works for both single and multiple mapping
                ens_idx = ens_idx[0]

            clust_marker_exp = z_scores.obs.groupby(partition_key)[ens_idx].apply(np.mean).tolist()
            clust_marker_exp.append(group)
            marker_exp.loc[i] = clust_marker_exp
            marker_names.append(gene)
            i+=1

    #Replace the rownames with informative gene symbols
    marker_exp.index = marker_names

    return(marker_exp)

def gene_expression(anndata, marker_list, gene_symbol_key=None, partition_key='leiden'):
    """A function to get mean z-score expressions of marker genes
     
    Inputs:
    
    anndata         - An AnnData object containing the data set and a partition
    marker_dict     - A dictionary with cell-type markers. The markers should be stores as anndata.var_names
    or an anndata.var field with the key given by the gene_symbol_key input
    gene_symbol_key - The key for the anndata.var field with gene IDs or names that correspond to the marker 
                      genes
    partition_key   - The key for the anndata.obs field where the cluster IDs are stored. The default is
                     'leiden' 
     
    Output:
    
    marker_exp: a dataframe with average z-score gene expression across the different partitions
    """

    #Test inputs
    if partition_key not in anndata.obs.columns.values:
        print('KeyError: The partition key was not found in the passed AnnData object.')
        print('   Have you done the clustering? If so, please tell pass the cluster IDs with the AnnData object!')
        raise

    if (gene_symbol_key != None) and (gene_symbol_key not in anndata.var.columns.values):
        print('KeyError: The provided gene symbol key was not found in the passed AnnData object.')
        print('   Check that your cell type markers are given in a format that your anndata object knows!')
        raise
        
    if gene_symbol_key:
        gene_ids = anndata.var[gene_symbol_key]
    else:
        gene_ids = anndata.var_names

    clusters = anndata.obs[partition_key].cat.categories
    n_clust = len(clusters)
    marker_exp = pd.DataFrame(columns=clusters)
    marker_names = []

    i = 0
    
    for gene in marker_list:
        ens_idx = np.in1d(gene_ids, gene) #Note there may be multiple mappings
        if np.sum(ens_idx) == 0:
            continue
        else:
            anndata.obs[ens_idx[0]] = anndata.X[:,ens_idx].mean(1) #works for both single and multiple mapping
            ens_idx = ens_idx[0]

        clust_marker_exp = anndata.obs.groupby(partition_key)[ens_idx].apply(np.mean).tolist()
        marker_exp.loc[i] = clust_marker_exp
        marker_names.append(gene)
        i+=1

    #Replace the rownames with informative gene symbols
    marker_exp.index = marker_names

    return(marker_exp)

def gene_expression_norm(anndata, marker_list, gene_symbol_key=None, partition_key='leiden'):
    """A function to get normalized expressions of marker genes
     
    Inputs:
    
    anndata         - An AnnData object containing the data set and a partition
    marker_dict     - A dictionary with cell-type markers. The markers should be stores as anndata.var_names or 
                      an anndata.var field with the key given by the gene_symbol_key input
    gene_symbol_key - The key for the anndata.var field with gene IDs or names that correspond to the marker 
                      genes
    partition_key   - The key for the anndata.obs field where the cluster IDs are stored. The default is
                      'leiden' 
    
    Output:
    
    marker_exp: a dataframe with average normalized gene expression across the different partitions
    """

    #Test inputs
    if partition_key not in anndata.obs.columns.values:
        print('KeyError: The partition key was not found in the passed AnnData object.')
        print('   Have you done the clustering? If so, please tell pass the cluster IDs with the AnnData object!')
        raise

    if (gene_symbol_key != None) and (gene_symbol_key not in anndata.var.columns.values):
        print('KeyError: The provided gene symbol key was not found in the passed AnnData object.')
        print('   Check that your cell type markers are given in a format that your anndata object knows!')
        raise
        
    if gene_symbol_key:
        gene_ids = anndata.var[gene_symbol_key]
    else:
        gene_ids = anndata.var_names

    clusters = anndata.obs[partition_key].cat.categories
    n_clust = len(clusters)
    marker_exp = pd.DataFrame(columns=clusters)
    marker_names = []
    
    i = 0
    
    ann_data = anndata
    
    for gene in marker_list:
        ens_idx = np.in1d(gene_ids, gene) #Note there may be multiple mappings
        if np.sum(ens_idx) == 0:
            continue
        else:
            ann_data.obs[ens_idx[0]] = ann_data.raw[:, ann_data.var.index].X[:,ens_idx].mean(1) #works for both single and multiple mapping
            ens_idx = ens_idx[0]

        clust_marker_exp = ann_data.obs.groupby(partition_key)[ens_idx].apply(np.mean).tolist()
        marker_exp.loc[i] = clust_marker_exp
        marker_names.append(gene)
        i+=1        
        
    #Replace the rownames with informative gene symbols
    marker_exp.index = marker_names

    return(marker_exp)

#Define cluster score for all markers
def evaluate_partition(anndata, marker_dict, gene_symbol_key=None, partition_key='leiden'):
    ''' This function gives a cell-type score for each partition key (i.e. Leiden clusters)
    Inputs:
    
    anndata         - An AnnData object containing the data set and a partition
    marker_dict     - A dictionary with cell-type markers. The markers should be stores as anndata.var_names or 
                    an anndata.var field with the key given by the gene_symbol_key input
    gene_symbol_key - The key for the anndata.var field with gene IDs or names that correspond to the marker 
                      genes
    partition_key   - The key for the anndata.obs field where the cluster IDs are stored. The default is
                      'leiden'
                      
    Outputs:
    
    A dataframe with a score for each cell type.
    '''

    if partition_key not in anndata.obs.columns.values:
        print('KeyError: The partition key was not found in the passed AnnData object.')
        print('   Have you done the clustering? If so, please tell pass the cluster IDs with the AnnData object!')
        raise

    if (gene_symbol_key != None) and (gene_symbol_key not in anndata.var.columns.values):
        print('KeyError: The provided gene symbol key was not found in the passed AnnData object.')
        print('   Check that your cell type markers are given in a format that your anndata object knows!')
        raise
        

    if gene_symbol_key:
        gene_ids = anndata.var[gene_symbol_key]
    else:
        gene_ids = anndata.var_names

    clusters = np.unique(anndata.obs[partition_key])
    n_clust = len(clusters)
    n_groups = len(marker_dict)
    
    marker_res = np.zeros((n_groups, n_clust))
    z_scores = sc.pp.scale(anndata, copy=True)

    i = 0
    for group in marker_dict:
        # Find the corresponding columns and get their mean expression in the cluster
        j = 0
        for clust in clusters:
            cluster_cells = np.in1d(z_scores.obs[partition_key], clust)
            marker_genes = np.in1d(gene_ids, marker_dict[group])
            marker_res[i,j] = z_scores.X[np.ix_(cluster_cells,marker_genes)].mean()
            j += 1
        i+=1

    variances = np.nanvar(marker_res, axis=0)
    if np.all(np.isnan(variances)):
        print("No variances could be computed, check if your cell markers are in the data set.")
        print("Maybe the cell marker IDs do not correspond to your gene_symbol_key input or the var_names")
        raise

    marker_res_df = pd.DataFrame(marker_res, columns=clusters, index=marker_dict.keys())

    #Return the median of the variances over the clusters
    return(marker_res_df)

def silhouette_analysis(range_n_clusters, X, metric='cosine', method="complete"):
    '''This function computes silhouette score for simple cosine clustering. 
    
    Inputs:
    
    range_n_clusters: list of a range of clusters, from 2 to (n-1) where n is the total number of samples.
    X: an input gene expression matrix
    metric: scipy cluster hierarchy defined distance metric
    method: scipy cluster hierarchy defined method
    
    Output:
    
    scores: a list of lists of silhouette scores, where each entry in the outside list holds a list of 
    silhouette scores for range_n_clusters.
    
    '''
    scores = []
    for n_clusters in range_n_clusters:
        i=0
        cluster_avg = []
        while i < 1: # can change this if you want more replicates
            clusterer = KMeans(n_clusters=n_clusters, random_state=10)
            cluster_labels = clusterer.fit_predict(X)

    # The silhouette_score gives the average value for all the samples.
    # This gives a perspective into the density and separation of the formed
    # clusters
            silhouette_avg = silhouette_score(X, cluster_labels, metric=metric)
            cluster_avg.append(silhouette_avg)
            i+=1
        scores.append((n_clusters, cluster_avg))
    return scores

def silhouette_plots(adata, pathway_names, pathway_genes, metric='cosine', norm = False, ax=None,key="leiden"):
    '''This function gives silhouette scores and scatterplots for different numbers of clusters for a specific
    pathway.
    
    Inputs: 
    adata : AnnData object
    pathway_names : string name of the pathways you're evaluating
    pathway_genes : a list of pathway genes.
    metric: scipy cluster hierarchy defined distance metric
    norm : whether or not to used z-score or normalized data. z-score is default.
    ax: a matplotlib axis object to plot onto
    key: a partition key in the .obs attribute of the adata object we want to cluster on.
    
    Outputs:
    
    fig: a scatterplot with the silhouette scores for 2-(n-1) clusters
    scores: a dataframe with silhouette scores for every value in range_n_clusters
    
    '''
    range_n_clusters = list(range(2,len(adata.obs[key].unique())))
    if ax == None:
        fig, ax = plt.subplots(figsize=(4,3.85))
    else:
        fig = Figure()
    scores = pd.DataFrame()
    if norm == False:
        df = gene_expression(adata, pathway_genes)
    else:
        df = gene_expression_norm(adata, pathway_genes)
    df = df.T
    df.reset_index(inplace=True)
    X = df[df.columns[1:]].values
    score = silhouette_analysis(range_n_clusters, X, metric=metric)
    scores[pathway_names] = [np.mean(s[1]) for s in score]
    ax.plot(range_n_clusters, scores, '-o', markersize=0.5);
    ax.set_title(pathway_names)
    ax.tick_params(axis="x", labelsize=10)
    ax.set_xlabel('Number of Clusters')
    ax.set_ylabel('Silhouette Score')
    new_dict = dict(zip(scores.index.values, range_n_clusters))
    scores.rename(new_dict)
    fig.tight_layout()
    return fig, scores

def heatmap(adata, pathway_genes, num_clust, name, metric='cosine', method="complete", norm = False, key="leiden",
                     figsize=(6,6)):
    '''We group the leiden clusters based on similarity of expression of 
    specific genes in a pathway. 
    
    Inputs:
    
    adata : the AnnData gene expression matrix
    pathway_genes : a list of the genes in the pathway
    num_clust : the optimal number of clusters based on silhouette score on
    a given metric (automatically cosine) distance
    name : the name with which we want to label the clusters of this pathway
    leg_axes : we can change the coordinates of the legend
    
    Outputs: 
    
    AnnData object labeled with the pathway clusters.
    g.fig Return at Index 0: Clustermap Figure you can later save
    df Return at Index 1: A dataframe of all the gene expression values
    '''
    if norm:
        df = gene_expression_norm(adata, pathway_genes,partition_key=key)
    else:
        df = gene_expression(adata, pathway_genes,partition_key=key)
        
    d = sch.distance.pdist(df.transpose(), metric=metric)
    L=sch.linkage(d)
    linkage = sch.fcluster(L, num_clust,'maxclust')
    str_linkage = []
    for i in linkage:
        str_linkage.append(str(i))
    new_dict = dict(zip(df.columns, str_linkage))

    adata.obs[name] = adata.obs[key].replace(new_dict).astype('category')
    if norm:
        df = gene_expression_norm(adata, pathway_genes,partition_key=key)
    else:
        df = gene_expression(adata, pathway_genes,partition_key=key)
    cols = {}
    for j in list(df.columns):
        cols[j] = str(adata[adata.obs[key] == j].obs[name][0])
    cols=pd.Series(data=cols, name='Clusters')
    labels = adata.obs[name].unique()
    labels = list(map(str, labels))
    cmap = plt.get_cmap('Paired')
    colors = cmap(np.linspace(0, 1, len(labels)))
    lut1 = dict(zip(labels, colors))
    cols_to_return = []
    keys_for_colors = list(lut1.keys())
    for k in keys_for_colors:
        cols_to_return.append(lut1[k])
    adata.uns[name+'_colors'] = cols_to_return
    row_colors1 = cols.map(lut1)
    g = sns.clustermap(df, metric=metric, row_cluster=False, cmap='viridis',
                      col_linkage=L, col_colors=row_colors1,figsize=figsize, linewidths=1, linecolor='black');
    ax = g.ax_heatmap
    g.fig.suptitle((name + ' with ' + str(num_clust) + ' clusters'), y=1.0,x=0.5,fontsize='large') 
    ax.set_xlabel('Clustering', x=0.5)
    return g.fig, df

def exp_across_clusters(df, ax = None):
    '''Plots a bar chart of expression summed across different Leiden clusters.
    This is useful to visualize which clusters we can remove from our heatmaps.
    
    Input: 
    
    df: A dataframe whose rows are different genes, and whose columns are the Leiden
    cluster labels.
    
    Output: a bar chart with the total expression for each cluster.'''
    
    cols = df.sum(axis=0)
    return cols.plot.bar(title='Expression Sum Across Leiden Clusters', ax = ax);

def exp_across_genes(df, ax = None):
    '''Plots a bar chart of expression summed across different genes.
    This is useful to visualize which genes we can remove from our heatmaps.
    
    Input: 
    
    df: A dataframe whose rows are different genes, and whose columns are the Leiden
    cluster labels.
    
    Output: a bar chart with the total expression for each cluster.'''
    
    rows = df.sum(axis=1)
    return rows.plot.bar(title='Expression Sum Across Genes', ax=ax);

def exp_above_threshold(df, axis, threshold, fig_ax = None):
    '''This function visualizes counts of gene expression above a certain threshold.
    
    Inputs: 
    
    df: A pandas dataframe whose rows are genes, and whose columns are Leiden 
    cluster labels.
    axis: 0=rows, 1=columns
    threshold: define a threshold value for gene expression.
    
    Output: a bar chart
    '''
    ax = df.iloc[:,:].ge(threshold).sum(axis).plot.bar(
        title=('Histogram of Counts Above Threshold Value'), ax=fig_ax )
    ax.set_xlabel('Sum along axis '+ str(axis))
    ax.set_ylabel('Counts Above Threshold Value')
    return ax

def diff_exp_counts_df(adata, dict_diff_exp_genes, key_exp_genes, groupby='leiden'):
    '''This function computs a dataframe with the number and identities of differentially
    expressed genes between each pair of group labels in an AnnData object. This should be
    used as input for a plot of differentially expressed genes.
    
    Inputs:
    
    adata: AnnData object
    dict_exp_genes: dictionary of differentially expressed genes for specific pathways. The
    key-value pairs are pathway names and a dictionary whose keys are tuples of observations 
    to group by and the values are the differentially expressed genes and their p-values 
    themselves. To get this dictionary run the diff_exp_between_clusters function BEFOREHAND.
    key_exp_genes: the key name in the dict_exp_genes dictionary, the specific pathway whose
    differentially expressed genes you want to view.
    
    Outputs:
    
    [1] A dataframe that serves as input for a Bokeh heatmap plot. This dataframe specifies
    the number, identity, and p-values of differentially expressed genes for a specific pair 
    of cluster groupings.
    [2] A list specifying order of the row labels for a Bokeh heatmap plot.
    [3] A list specifying order of the column labels for a Bokeh heatmap plot.
    '''
    d1 = {}
    num_clust = len(adata.obs[groupby].unique())
    for i in range(0,num_clust):
        d1[i] = []

    for i in sorted(list(dict_diff_exp_genes[key_exp_genes].keys())):
        d1[i[0]].append(dict_diff_exp_genes[key_exp_genes][i])
    counts_d1 = {}
    genes_d1 = {}
    for i in range(0, num_clust):
        l = [0]*(i+1) 
        l = l + list(len(j) for j in d1[i])
        counts_d1[str(i)] = l
        k = [[]]*(i+1)
        k = k + list(j for j in d1[i])
        genes_d1[str(i)] = k
    counts_df1 = pd.DataFrame.from_dict(counts_d1)
    counts_df1.index = [str(i) for i in range(0,19)]
    genes_d1 = pd.DataFrame.from_dict(genes_d1)
    genes_d1.index = [str(i) for i in range(0,19)]
    heatmap_sns = sns.clustermap(counts_df1, cmap="viridis")
    plt.close()
    colname_list = [counts_df1.columns[col_id] for col_id in 
    heatmap_sns.dendrogram_col.reordered_ind]
    rowname_list = [counts_df1.index[row_id] for row_id in 
    heatmap_sns.dendrogram_row.reordered_ind]
    counts = counts_df1.stack().reset_index()
    counts.columns = ['x_val', 'y_val', 'counts']
    genes_df = genes_d1.stack().reset_index()
    genes_df.columns = ['x_val', 'y_val', 'genes']
    counts['genes'] = genes_df['genes']
    genes = []
    for j in list(counts.index):
        new_g = []
        for i in counts.loc[j]['genes']:
            if i[1] < 0.009:
                new_str = ": {:.3e}".format(i[1])
                new_g.append(i[0]+ new_str)
            else:
                new_g.append(i[0]+": {:.3}".format(i[1]))
        genes.append(new_g)
    counts['genes_formatted_for_display'] = genes
    counts.genes_formatted_for_display = counts.genes_formatted_for_display.transform(lambda x: "<br>".join(x))
    return counts, rowname_list, colname_list

def diff_exp_between_clusters(adata, dict_pathway_genes,key='leiden', method='wilcoxon'):
    '''This function outputs the differentially expressed genes and their p-values for 
    specific pathways in a given AnnData object. It computes these differentially expressed
    genes for pairwise combinations of a specific AnnData object observation. 
    
    Input:
    
    adata: AnnData object
    dict_pathway_genes: Dictionary whose key-value pairs are pathway names and a list of genes
    corresponding to the pathway.
    key: what anndata.obs grouping to use. 
    
    Returns
    
    A dictionary of dictionaries whose keys are pathway names and the values are dictionaries 
    whose keys are tuples of AnnData observation key groupings, and values are lists of tuples 
    specifying (pathway gene, p-value). The function also prints out a Bonferroni correction factor.
    '''
    # Just don't want to print out ranking genes information every time we run it.
    sc.settings.verbosity = 1 
    # Get the different number of combinations of leiden pairs.
    comb = list(combinations(adata.obs[key].unique(), 2))
    
    diff_exp={}
    for i in dict_pathway_genes.keys():
        diff_exp[i] = {}
        
    # Iterating through each possible combination for Leiden cluster pairs and doing the 
    # gene ranking. 
    for c in comb:
        results=[]
        sc.tl.rank_genes_groups(adata, groupby=key, groups=[c[0], c[1]], 
                            method=method, n_genes=adata.shape[1]);
        for j in diff_exp.keys():
            gene_mask1 = [gene in dict_pathway_genes[j] for gene in adata.uns['rank_genes_groups']['names'][str(c[0])]]
            gene_mask2 = [gene in dict_pathway_genes[j] for gene in adata.uns['rank_genes_groups']['names'][str(c[1])]]
            names1 = adata.uns['rank_genes_groups']['names'][gene_mask1]
            pvals1 = adata.uns['rank_genes_groups']['pvals'][gene_mask1]
            names2 = adata.uns['rank_genes_groups']['names'][gene_mask2]
            pvals2 = adata.uns['rank_genes_groups']['pvals'][gene_mask2]
            for i in range(0, len(names1)):
                results.append((names1[i][0], pvals1[i][0]))
            for i in range(0, len(names2)):
                results.append((names2[i][1], pvals2[i][1]))
            diff_exp[j][c] = results
            results=[]
    sc.settings.verbosity = 3  # resetting the verbosity.
    print("Since we are conducting the same test multiple times, our signficance level \
    needs to be divided by (if we do a Bonferroni Correction): "+str(len(comb)))
    return diff_exp

def trial_diff_exp_between_clusters(c, adata, dict_pathway_genes,key='leiden', method='wilcoxon'):
    '''This function outputs the differentially expressed genes and their p-values for 
    specific pathways in a given AnnData object. It computes these differentially expressed
    genes for pairwise combinations of a specific AnnData object observation. 
    
    Input:
    
    adata: AnnData object
    dict_pathway_genes: Dictionary whose key-value pairs are pathway names and a list of genes
    corresponding to the pathway.
    key: what anndata.obs grouping to use. 
    
    Returns
    
    A dictionary of dictionaries whose keys are pathway names and the values are dictionaries 
    whose keys are tuples of AnnData observation key groupings, and values are lists of tuples 
    specifying (pathway gene, p-value). The function also prints out a Bonferroni correction factor.
    '''
    # Just don't want to print out ranking genes information every time we run it.
    sc.settings.verbosity = 1 
    # Get the different number of combinations of leiden pairs.
    #comb = list(combinations(range(0,len(adata.obs[key].unique())), 2))
    
    diff_exp={}
    for i in dict_pathway_genes.keys():
        diff_exp[i] = {}
        
    # Iterating through each possible combination for Leiden cluster pairs and doing the 
    # gene ranking. 
    results=[]
    sc.tl.rank_genes_groups(adata, groupby=key, groups=[c[0], c[1]], 
                        method=method, n_genes=adata.shape[1]);
    for j in diff_exp.keys():
        gene_mask1 = [gene in dict_pathway_genes[j] for gene in adata.uns['rank_genes_groups']['names'][str(c[0])]]
        gene_mask2 = [gene in dict_pathway_genes[j] for gene in adata.uns['rank_genes_groups']['names'][str(c[1])]]
        names1 = adata.uns['rank_genes_groups']['names'][gene_mask1]
        pvals1 = adata.uns['rank_genes_groups']['pvals'][gene_mask1]
        names2 = adata.uns['rank_genes_groups']['names'][gene_mask2]
        pvals2 = adata.uns['rank_genes_groups']['pvals'][gene_mask2]
        for i in range(0, len(names1)):
            results.append((names1[i][0], pvals1[i][0]))
        for i in range(0, len(names2)):
            results.append((names2[i][1], pvals2[i][1]))
        diff_exp[j] = results
        results=[]
    sc.settings.verbosity = 3  # resetting the verbosity.
    #print("Since we are conducting the same test multiple times, our signficance level \
    #needs to be divided by (if we do a Bonferroni Correction): "+str(len(comb)))
    return diff_exp

def weight_dist(ge_df, weight_df, metric = "cosine"):
    '''Returns the flattened distance matrix of distances between leiden clusters based on 
    gene expression.
    
    Inputs:
    
    ge_df: A dataframe with rows as genes and columns as cluster labels.
    weight_df: A dataframe with lists of weights specifying differentially expressed genes
    for all combinations of clusters specified in ge_df.
    
    Output:
    
    return_list: A 1-D condensed distance matrix for distances between all pairs of cluster labels.
    
    '''
    ge_df = ge_df.transpose()
    combs = list(combinations(ge_df.index, 2))
    return_list = np.zeros(len(combs))
    if metric == "cosine":
        for idx,i in enumerate(combs):
            w1 = np.dot(ge_df.loc[i[0]], weight_df.loc[i[0]][i[1]])
            w2 = np.dot(ge_df.loc[i[1]], weight_df.loc[i[0]][i[1]])
            if 1.0 not in weight_df.loc[i[0]][i[1]]:
                dist = 0
            elif (w1==0.0 or w2==0.0):
                dist = 0
            else:
                dist = cosine(ge_df.loc[i[0]], ge_df.loc[i[1]], weight_df.loc[i[0]][i[1]])
            return_list[idx] = dist
    if metric == "euclidean":
        for idx,i in enumerate(combs):
            w1 = np.dot(ge_df.loc[i[0]], weight_df.loc[i[0]][i[1]])
            w2 = np.dot(ge_df.loc[i[1]], weight_df.loc[i[0]][i[1]])
            if 1.0 not in weight_df.loc[i[0]][i[1]]:
                dist = 0
            elif (w1==0.0 or w2==0.0):
                dist = 0
            else:
                dist = euclidean(ge_df.loc[i[0]], ge_df.loc[i[1]], weight_df.loc[i[0]][i[1]])
            return_list[idx] = dist
    return return_list

def weighted_silhouette_plot(adata, pathway_names,
                             pathway_genes, weight_df, key = "leiden", metric = "cosine",
                             norm = True, ax=None):
    '''This function gives silhouette scores and scatterplots for clusters based on our self-defined distance
    metric.
    
    Inputs:
    adata : AnnData object
    pathway_names : string name of the pathways you're evaluating
    pathway_genes : a list of pathway genes.
    weight_df: a dataframe with a list of dataframes with weights for distance between two clusters based on p-
    values
    key: a partition key in the .obs attribute of the adata object we want to cluster on.
    norm : whether or not to used z-score or normalized data. z-score is default.
    ax: a matplotlib axis object to plot onto
    
    Outputs:
    
    fig: a scatterplot with the silhouette scores for 2-(n-1) clusters
    sil_score_list: a dataframe with silhouette scores for every value in range_n_clusters
    
    '''
    fig, ax = plt.subplots(figsize=(5,5))
    range_n_clusters = list(range(2,len(adata.obs[key].unique())))
    new_range_n_clusters = []
    sil_score_list = []
    if norm:
        df = gene_expression_norm(adata, pathway_genes, partition_key = key)
    else:
        df = gene_expression(adata, pathway_genes, partition_key = key)
    d = weight_dist(df, weight_df, metric = metric)
    L=sch.linkage(d, metric='euclidean', method='complete')
    X = squareform(d)
    for i in range_n_clusters:
        linkage = sch.fcluster(L, i,'maxclust')
        if not np.array_equal(linkage,np.array([1]*len(adata.obs[key].unique()))):
            new_range_n_clusters.append(i)
            sil_score = silhouette_score(X, linkage, metric="precomputed", method='complete')
            sil_score_list.append([sil_score])
            sil_score = []
    ax.plot(new_range_n_clusters, sil_score_list, '-o', markersize=0.75);
    ax.tick_params(axis="x", labelsize=8);
    return fig, sil_score_list

def weighted_heatmap(adata, pathway_genes, pathway_pvals, num_clust, name, key = "leiden",norm = True,
                     leg_axes = (1.3, 1.3), leg_cols = 1, figsize=(10,6),
                    metric = "euclidean"):
    '''We group the leiden clusters based on similarity of expression of specific genes in a pathway. We weight 
    this similarity based on differentially expressed genes between all leiden clusters.
    
    Inputs :
    
    adata : the AnnData gene expression matrix
    pathway_genes : a list of the genes in the pathway
    pathway_pvals: a lower traingular dataframe that gives weights for differentially expressed genes for all
    combinations of labels.
    num_clust : the optimal number of clusters based on silhouette score on
    a given metric (automatically cosine) distance
    name : the name with which we want to label the clusters of this pathway
    leg_axes : we can change the coordinates of the legend
    
    Returns: 
    
    AnnData object labeled with the pathway clusters.
    Return at Index 0: Clustermap Figure you can later save
    Return at Index 1: A dataframe of all the gene expression values
    '''
    if norm:
        df = gene_expression_norm(adata, pathway_genes, partition_key = key)
    else:
        df = gene_expression(adata, pathway_genes, partition_key = key)
    d = weight_dist(df, pathway_pvals, metric = metric)
    L=sch.linkage(d, metric="euclidean", method='complete')
    linkage = sch.fcluster(L, num_clust,'maxclust')
    str_linkage = []
    for i in linkage:
        str_linkage.append(str(i))
    new_dict = dict(zip(list(df.columns), str_linkage))
    adata.obs[name] = adata.obs[key].replace(new_dict)
    if norm:
        df = gene_expression_norm(adata, pathway_genes, partition_key = key)
    else:
        df = gene_expression(adata, pathway_genes, partition_key = key)
    cols = {}
    for j in list(df.columns):
        cols[j] = str(adata[adata.obs[key] == j].obs[name][0])
    cols=pd.Series(data=cols, name='Clusters')
    labels = adata.obs[name].unique()
    labels = list(map(str, labels))
    cmap = plt.get_cmap('Paired')
    colors = cmap(np.linspace(0, 1, len(labels)))
    lut1 = dict(zip(labels, colors))
    cols_to_return = []
    keys_for_colors = [int(i) for i in lut1.keys()]
    keys_for_colors.sort()
    for k in keys_for_colors:
        cols_to_return.append(lut1[str(k)]) 
    adata.uns[name+'_colors'] = cols_to_return
    row_colors1 = cols.map(lut1)
    g = sns.clustermap(df, row_cluster=False, cmap='viridis',
                      col_linkage=L, col_colors=row_colors1,figsize=figsize, linewidths=1, linecolor='black');
    ax = g.ax_heatmap
    g.fig.suptitle((name + ' with ' + str(num_clust) + ' clusters'), y=1.0,x=0.5,fontsize='large') 
    ax.set_xlabel('Labels', x=0.5)
    return g.fig, df