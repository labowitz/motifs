# motifs
Repository for "Combinatorial expression motifs in signaling pathways" paper by Granados, Kanrar, Elowitz (2023).

## Data and code availability

All code is available on GitHub at https://github.com/nivkanrar/motifs.
All data files and code outputs (.pdf, .csv, .h5ad, .RDS, etc.) are available on Google Drive https://drive.google.com/drive/folders/1tiJ0c8OAuk-Nh4IXGUILUVPLRvG2lKxQ?usp=sharing. 

## Folder structure 

./module/: low-level functions 

./data/raw_data/: includes un-processed single-cell datasets and metadata files. 

./data/processed_data/: individual processed single-cell datasets along with the integrated object. Generic files for color palettes, genes lists, and pathways lists.

./scripts/analysis/: scripts to generate the main figures in the manuscript.

./scripts/figures/: figures generated from scripts in ./scripts/analysis

./scripts/preprocessing/: scripts to do preprocessing and integration of data in ./data/raw_data/ folder, with results stored in ./data/processed_data/
## Key data files 

Path: ./data/processed_data/master_seurat.RDS 

Description: R data file containing a Seurat object with transcriptomic profiles for all 1206 cell types

Path: ./data/pathbank/pathway_df.csv

Description: csv file containing the pathway names and their components from PathBank after filtering

Path: ./data/raw_data/sc_data/

Description: folder with the raw data files (counts and metadata files) for each single-cell dataset part of the integrated atlas

## Main scripts 

Path: ./scripts/analysis/Figure_2C_3_4BCDE_S2ABC_S3D.R

Description: main script to run the motif analysis on the TGF-beta receptor genes. Generated Figures 2C, 3, 4BCDE, S2ABC, and S3D, pertaining to all analyses on TGF-beta receptor genes.

Path: ./module/module_new.R

Description: low-level functions for motif analysis, which are imported into other scripts. 

## Downloading the data from Google Drive
For large files, Google Drive compresses everything into multiple zip files. You can extract all zip files combined in Linux/Mac using the following: 
```
# move all zip files into a new directory
mkdir data
unzip ´*.zip´ -d data # will extract all zip files in the current directory
```

## How to run the code 
1. Set your file paths for data objects in ./scripts/analysis/imports.R 
- Data objects 
- Set the output folders 
- List of genes and pathways 
- Color palettes 

2. Run scripts/analysis/Figure_2C_4BCDE_S2ABC_S3D.R
- This script will generate all the plots for TGF-β as they appear in the manuscript
- This script shows a full example of the motif analysis for one pathway
- We adapted the same pipeline to run multiple pathways in parallel

## Installations
We use code written in both R and Python for our project.
### RStudio/R Code
It’s easiest to run our R code using RStudio. To isolate the package versions used in this project, create an RStudio project in the main directory for this project, and always run the code in this RStudio project.
R does not have a clear-cut way to share environments in a replicable way. We share the package versions we have installed in `./installed_packages.txt`. But due to package installation conflicts in R, it’s probably easiest to try a one-shot installation of the packages, which you can do by simply opening the `./module/module.R` script in RStudio and clicking on the prompt to install the unavailable packages.
### Python Code
We used the installations available on the Anaconda package manager for Python.
Install Python 3.7 Anaconda Shell Installer
Open terminal and follow the below steps.
```
conda create --name motifs -c conda-forge rpy2 (3.3.2)
conda install -c conda-forge pandas numpy=1.18.5 scipy scikit-learn jupyter gsl tzlocal simplegeneric natsort h5py tqdm patsy llvmlite numba networkx joblib numexpr pytables seaborn statsmodels
conda install -c conda-forge python-igraph leidenalg
pip install anndata anndata2ri fa2 gprofiler-official scanpy
R
install.packages(c('devtools', 'gam', 'RColorBrewer', 'BiocManager', 'plotly'))
update.packages(ask=F)
q()
conda install -c r r-xml
R
install.packages("data.table", type = "source", repos = "https://Rdatatable.gitlab.io/data.table")
BiocManager::install(c("scran", "MAST", "slingshot", "ComplexHeatmap", "Seurat", "tradeSeq", "DEsingle"), version="3.10")
```
