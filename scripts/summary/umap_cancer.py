#!/gpfs/data/yarmarkovichlab/Frank/BayesTS/logit_gate_env/bin/python3.7

import pandas as pd
import numpy as np
import sys,os
import subprocess
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from ast import literal_eval
from tqdm import tqdm
import math
import anndata as ad
import scanpy as sc
from sctriangulate import *
from sctriangulate.preprocessing import *

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

cancers = [
    'BRCA',
    'KIRC',
    'COAD',
    'STAD',
    'MESO',
    'LIHC',
    'ESCA',
    'CESC',
    'BLCA',
    'RT',
    'AML',
    'DLBC',
    'GBM',
    'NBL',
    'PAAD',
    'HNSC',
    'OV',
    'LUSC',
    'LUAD',
    'CHOL',
    'SKCM',
]  

n_samples = [
    1118,
    542,
    483,
    412,
    87,
    374,
    185,
    306,
    412,
    63,
    151,
    48,
    170,
    157,
    179,
    522,
    429,
    502,
    541,
    35,
    472
]


# use peptide abudance to cluster immunopeptidome samples
rootdir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas'
for c in cancers:
    final_all = pd.read_csv(os.path.join(rootdir,c,'antigen','fdr','final_enhanced_full.txt'),sep='\t')
    cond = [False if ('[]' in item) and ('(\'HLA-' not in item) else True for item in final_all['presented_by_each_sample_hla']]
    final_all = final_all.loc[cond,:]
    final_all = final_all.loc[final_all['typ']=='self_gene',:]
    final_all = final_all.loc[final_all['unique'],:]
    all_peptides = set(final_all['pep'].values)
    msms = pd.read_csv(os.path.join(rootdir,c,'antigen','fdr','msmsScans_all_add_tesorai.txt'),sep='\t')
    msms = msms.loc[(msms['Identified']=='+') & (msms['Sequence'].isin(all_peptides)),:]
    msms['uid'] = [i1 + ';' + i2 for i1,i2 in zip(msms['level_0'],msms['Raw file'])]





# use transcriptome to cluster TCGA
rootdir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas'
genes = set()
for c in cancers:
    final_all = pd.read_csv(os.path.join(rootdir,c,'antigen','fdr','final_enhanced_full.txt'),sep='\t')
    cond = [False if ('[]' in item) and ('(\'HLA-' not in item) else True for item in final_all['presented_by_each_sample_hla']]
    final_all = final_all.loc[cond,:]
    final_all = final_all.loc[final_all['typ']=='self_gene',:]
    final_all = final_all.loc[final_all['unique'],:]
    genes = genes.union(set(final_all['ensgs'].values))

gene_df_list = []
for c in cancers:
    gene_tpm = pd.read_csv(os.path.join(rootdir,c,'gene_tpm.txt'),sep='\t',index_col=0)
    common = list(genes.intersection(set(gene_tpm.index)))
    gene_tpm = gene_tpm.loc[common,:]
    gene_tpm = gene_tpm.loc[np.logical_not(gene_tpm.index.duplicated()),:]
    gene_df_list.append(gene_tpm)
df = pd.concat(gene_df_list,axis=1,keys=cancers)
df = df.dropna(axis=0,how='any')

mi = df.columns
mi_df = mi.to_frame(index=False)
df.columns = mi_df[1].values

df = df.T
adata = ad.AnnData(X=df.values,obs=pd.DataFrame(index=df.index,data={'cancer':mi_df[0].values}),var=pd.DataFrame(index=df.columns))  # 7058 × 13033
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata,flavor='seurat',n_top_genes=5000)
adata.raw = adata
adata = adata[:,adata.var['highly_variable']]
sc.pp.scale(adata,max_value=10)
sc.tl.pca(adata,n_comps=50)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
umap_dual_view_save(adata,cols=['cancer'])





