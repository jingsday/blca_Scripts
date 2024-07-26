import scanpy as sc
import anndata as an
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanorama
from scipy.sparse import csr_matrix

from pathlib import Path
import os

sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

# Define the directory path
directory = '/home/jyang/Phd_project/project_UCD_blca/blca_DATA/blca_DATA_mouse_GSE174182_RAW'

# Get the list of files in the directory (non-recursive)
dirs = os.listdir(directory)

# Create an empty list to store the names
names_list = []

# Extract the unique names from the first 20 characters of the filenames
for x in dirs:
    name = x[:20]
    names_list.append(name)

# Remove duplicates by converting to a set and then back to a list
names_list = list(set(names_list))

# Print the unique names
print(names_list)

os.chdir(directory)

names_list=['GSM5288669_Sample-4_', 'GSM5288670_Sample-5_', 'GSM5288671_Sample-6_', 
            'GSM5288668_Sample-3_', 'GSM5288672_Sample-7_', 'GSM5288673_Sample-8_', 'GSM5288674_Sample-11_']

adata_list = []

# Loop over each sample and read in the AnnData object
for name in names_list:

    mtx =f"{name}filtered_matrix.mtx.gz"
    adata = sc.read_mtx(mtx)
    cells=pd.read_csv(f'{name}filtered_barcodes.tsv.gz',header=None)
    features=pd.read_csv(f'{name}filtered_features.tsv.gz',header=None,sep='\t')
    adata= adata.T

    #check the columns first to make sure they are the ones you need 
    adata.obs['CellID']= cells[0].tolist()
    adata.var['Gene']= features[1].tolist()
    adata.var.index= adata.var['Gene']
    adata.var_names_make_unique() 

    sc.pp.filter_cells(adata, min_genes=300)
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .97)
    adata = adata[adata.obs.n_genes_by_counts < upper_lim]

    adata.obs['source'] = name[:10]
      
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000, inplace=True)
    
    adata_list.append(adata)

    # Assuming adatas is a list of AnnData objects
adatas_csr = []
for ad in adata_list:
    # Convert the X attribute to CSR format if it's not already in CSR format
    if not isinstance(ad.X, csr_matrix):
        ad.X = ad.X.tocsr()
    adatas_csr.append(ad)

# Now run scanorama.correct_scanpy with the CSR formatted AnnData objects
adatas_cor = scanorama.correct_scanpy(adatas_csr, return_dimred=True)

batch_names = [adata.obs['source'].iloc[0] for adata in adatas_cor]

adatas_cor_full = adatas_cor[0].concatenate(adatas_cor[1:], batch_key='source', batch_categories=batch_names)                                  

adatas_cor_full

os.listdir('/home/jyang/Phd_project/project_UCD_blca/blca_OUTPUT/blca_OUTPUT_scanorama')

adatas_cor_full.write_h5ad('/home/jyang/Phd_project/project_UCD_blca/blca_OUTPUT/blca_OUTPUT_scanorama_cor_full.h5ad')

sc.pp.neighbors(adatas_cor_full, use_rep="X_scanorama")
sc.tl.umap(adatas_cor_full)

#Leiden resolution 0.1
sc.tl.leiden(
    adatas_cor_full, key_added="leiden01", n_iterations=2, flavor="igraph", directed=False,resolution=0.1
)

sc.pl.umap(
    adatas_cor_full, color=["leiden01", "source"], palette=sc.pl.palettes.default_20,save='blca_scanorama_umap_leiden01_source.svg'
)
sc.tl.rank_genes_groups(adatas_cor_full, "leiden01", method="t-test")
sc.pl.rank_genes_groups(adatas_cor_full, n_genes=25, sharey=False,save='blca_scanorama_markers_leiden01.svg')


sc.pl.umap(adatas_cor_full,color=['leiden01','Upk2','Upk1a','Upk1b','Cdh1','Upk3a','Upk3b','Ivl'],save='blca_scanorama_umap_leiden01_markers.svg')
sc.pl.violin(adatas_cor_full,['Upk2','Upk1a','Upk1b','Cdh1','Upk3a','Upk3b','Ivl'],groupby='leiden01',save='blca_scanorama_violin_leiden01_markers.svg')

#Leiden resolution 0.05
sc.tl.leiden(
    adatas_cor_full, key_added="leiden005", n_iterations=2, flavor="igraph", directed=False,resolution=0.05
)

sc.pl.umap(
    adatas_cor_full, color=["leiden005", "source"], palette=sc.pl.palettes.default_20,save='blca_scanorama_umap_leiden005_source.svg'
)

sc.tl.rank_genes_groups(adatas_cor_full, "leiden005", method="t-test")
sc.pl.rank_genes_groups(adatas_cor_full, n_genes=25, sharey=False,save='blca_scanorama_markers_leiden005.svg')

sc.pl.umap(adatas_cor_full,color=['leiden005','Upk2','Upk1a','Upk1b','Cdh1','Upk3a','Upk3b','Ivl'],save='blca_scanorama_umap_leiden005_markers.svg')
sc.pl.violin(adatas_cor_full,['Upk2','Upk1a','Upk1b','Cdh1','Upk3a','Upk3b','Ivl'],groupby='leiden005',save='blca_scanorama_violin_leiden005_markers.svg')

#Leiden resolution 0.01
sc.tl.leiden(
    adatas_cor_full, key_added="leiden001", n_iterations=2, flavor="igraph", directed=False,resolution=0.01
)

sc.pl.umap(
    adatas_cor_full, color=["leiden001", "source"], palette=sc.pl.palettes.default_20,save='blca_scanorama_umap_leiden001.svg'
)

sc.tl.rank_genes_groups(adatas_cor_full, "leiden001", method="t-test")
sc.pl.rank_genes_groups(adatas_cor_full, n_genes=25, sharey=False,save='blca_scanorama_markers_leiden001.svg')

sc.pl.umap(adatas_cor_full,color=['leiden001','Upk2','Upk1a','Upk1b','Cdh1','Upk3a','Upk3b','Ivl'],save='blca_scanorama_umap_leiden001_markers.svg')
sc.pl.violin(adatas_cor_full,['Upk2','Upk1a','Upk1b','Cdh1','Upk3a','Upk3b','Ivl'],groupby='leiden001',save='blca_scanorama_violin_leiden001_markers.svg')

#cluster counts 

print(adatas_cor_full.obs['leiden01'].value_counts(ascending=False))

print(adatas_cor_full.obs['leiden005'].value_counts(ascending=False))

print(adatas_cor_full.obs['leiden001'].value_counts(ascending=False))