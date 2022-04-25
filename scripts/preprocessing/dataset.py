# Set directories for the data, results, and figures.
datadir = '../../data/raw_data/'
figdir = '../figures/'
resdir = '../../data/processed_data/kharchenko_trunk_neural_crest/'

# Import the atlas dataset of counts
counts = pd.read_csv(datadir + "fig_version_countsSeurat_Apr22.csv", index_col = [0]).T

# Import metadata csv
meta = pd.read_csv(datadir + "meta_data_Aug2021.csv", index_col = [0])

# Import reference pathway list csv
all_pathways = pd.read_csv(datadir + "allPathways_listGenes_dec2021.tsv", delimiter="\t")

# Normalize the atlas counts to have 1e4 per cluster
integrated_counts = np.log1p(counts.div(counts.sum(axis=1),axis=0)*1e4)

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
    new_clip_sample = np.clip(new_sample,  a_min=0, a_max = gene_quantiles)
    # min.max scaler -- same as for the integrated dataset 
    new_sample = new_clip_sample.values 
    new_sample = new_sample.reshape(1,-1)

    test_sample =scaler.transform(new_sample)
    return test_sample

# plots rows as independent profiles 
# use only with a reduced number of data.points 
def barplot_profiles(top_cluster = [], pathway_genes = [], all_id_vars = ['cell_id'], hue_group = 'cell_id'): 
    cluster_tidy  =pd.melt(top_cluster, id_vars = all_id_vars, 
            value_vars = pathway_genes, var_name = 'gene', value_name = 'expression')

    plt.figure(figsize = (8,4), dpi = 300)
    g = sns.barplot(data = cluster_tidy, x ='gene',y = 'expression',  hue = hue_group,
               palette = "Set2")
    g.legend(loc='center left', bbox_to_anchor = (1, 0.5))
    plt.xticks(rotation = 90)
    plt.tight_layout() 
    
    
# finds correlation and cosine similarity for a target profile against a data.frame of profiles 
def find_profile_match(df = [], pathway_genes = [], target_profile = []):
    corr_coefs = df[pathway_genes].apply( lambda x: np.corrcoef(x, target_profile)[0,1], axis =1)
    cosine_sim = df[pathway_genes].apply( lambda x: distance.cosine(x, target_profile), axis =1)

    df['corr'] = corr_coefs
    df['cosine'] = 1-cosine_sim

    return df

# filter expressing cell types for this pathway 
# from the original integrated matrix: 
# the pipeline should be able to run from scratch starting here. 
def filter_on_celltypes(integrated_counts = [], pathway_genes =[], min_expr = 0.2, min_genes_on = 2):
    x = integrated_counts[pathway_genes]
    gene_quantiles = np.quantile(x, q= 0.99, axis = 0)

    x_clipped = quantile_saturation(x, gene_quantiles)

    # min max after saturation
    scaler= MinMaxScaler() 
    scaler.fit(x_clipped)
    # max value for each gene = 1
    df_clip = pd.DataFrame(scaler.transform(x_clipped), columns = pathway_genes)
    df_clip['cell_id'] = meta['cell_id'].values

    df_clip['ON'] = np.sum(df_clip[pathway_genes].values>min_expr, axis =1)>min_genes_on

    pathway_mat = df_clip[df_clip['ON']] # final data.frame min.max scaled 
    return pathway_mat 

# Process the integrated atlas so counts are MinMaxed.
x = integrated_counts.copy()
gene_quantiles = np.quantile(x, q=0.99, axis=0) # compute quantiles
scaler = MinMaxScaler()
x_clipped = quantile_saturation(x, gene_quantiles)
scaler.fit(x_clipped)

integrated_minmax = pd.DataFrame(scaler.transform(x_clipped), columns = x.columns)
integrated_minmax['cell_id'] = meta['cell_id'].values

