#!/gpfs/data/yarmarkovichlab/Frank/BayesTS/logit_gate_env/bin/python3.7

import os,sys
import pandas as pd
import numpy as np
import anndata as ad
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import pearsonr,ttest_ind,ttest_rel,spearmanr
import pickle
from scipy.sparse import csr_matrix
import scanpy as sc
from sctriangulate import *
from sctriangulate.preprocessing import *
from sctriangulate.colors import bg_greyed_cmap

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

# read in h5ad
adata_rna = sc.read('RCC_upload_final_raw_counts.h5ad')
adata_rna.var['mt'] = adata_rna.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata_rna, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
# sc.pp.filter_cells(adata_rna, min_genes=300)
# sc.pp.filter_cells(adata_rna, min_counts=500)
# adata_rna = adata_rna[adata_rna.obs.pct_counts_mt < 20, :]
sc.pp.log1p(adata_rna)

adata_rna = adata_rna[adata_rna.obs['broad_type']=='RCC',:]


sc.pl.umap(adata_rna,color=['CA9','HAVCR1'],cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
plt.savefig('markers_ccRCC.pdf',bbox_inches='tight')
plt.close()

df = adata_rna[:,['CA9','HAVCR1']].to_df()
df.to_csv('markers_ccRCC.txt',sep='\t')

fig,ax = plt.subplots()
ax.scatter(df['CA9'].values,df['HAVCR1'].values)
print(spearmanr(df['CA9'].values,df['HAVCR1'].values))
plt.savefig('markers_ccRCC_scatter.pdf',bbox_inches='tight')
plt.close()

