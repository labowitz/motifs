# motifs
Repository for "Combinatorial expression motifs in signaling pathways" paper by Granados, et al. (2022).

## Data and code availability

All code is available on GitHub at https://github.com/nivkanrar/motifs.
All data files and code outputs (.pdf, .csv, .h5ad, .RDS, etc.) are available on Google Drive https://drive.google.com/drive/folders/1tiJ0c8OAuk-Nh4IXGUILUVPLRvG2lKxQ?usp=sharing. 

## Folder structure 
./data/processed_data: individual processed single-cell datasets along with the integrated object. Generic files for color palettes, genes lists, and pathways lists.
 
./data/raw_data/: includes un-processed single-cell datasets and metadata files. 

./module/: low-level functions 

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
