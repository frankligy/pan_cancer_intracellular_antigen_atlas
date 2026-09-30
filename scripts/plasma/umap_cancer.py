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
from sctriangulate.colors import *

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


# first test yield
gbm = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/plasma/plasma_gbm_PXD008127/antigen/other_alg/all_neoantigen.txt',sep='\t')
lung_public = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/plasma/lung_cancer_plasma_public/antigen/other_alg/all_neoantigen.txt',sep='\t')
lung_nyu = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/plasma/lung_cancer_nyu/antigen/other_alg/all_neoantigen.txt',sep='\t')

yields = []
samples = []
for s,sub_df in gbm.groupby(by='Raw file'):
    yields.append(len(set(sub_df['Sequence'])))
    samples.append(s)
gbm_yields, gbm_samples = zip(*sorted(zip(yields,samples),reverse=True))

lung_public_yields = [len(set(lung_public['Sequence']))]
lung_public_samples = ['lung_public']

yields = []
samples = []
for s,sub_df in lung_nyu.groupby(by='Raw file'):
    yields.append(len(set(sub_df['Sequence'])))
    samples.append(s)
lung_nyu_yields, lung_nyu_samples = zip(*sorted(zip(yields,samples),reverse=True))

total_yields = list(gbm_yields) + lung_public_yields + list(lung_nyu_yields)
total_samples = list(gbm_samples) + lung_public_samples + list(lung_nyu_samples)

mean_yields = np.quantile(total_yields,[0.25,0.5,0.75])

fig,ax = plt.subplots(figsize=(15,4.8))
ax.bar(x=np.arange(len(total_samples)),height=total_yields,color=['#6BAED6']*len(gbm_samples)+['#74A892']*len(lung_public_samples)+['#E38C78']*len(lung_nyu_samples))
ax.axhline(y=mean_yields[0],linestyle='--',color='#5B6573')
ax.axhline(y=mean_yields[1],linestyle='--',color='#5B6573')
ax.axhline(y=mean_yields[2],linestyle='--',color='#5B6573')
ax.set_ylabel('Number of unique peptides identified per sample')
ax.set_title(mean_yields)
ax.set_xticks(np.arange(len(total_samples)))
ax.set_xticklabels(total_samples,fontsize=2,rotation=90)
plt.savefig('yield.pdf',bbox_inches='tight')
plt.close()

# generate reactome
result = pd.read_csv('Reactome_Pathways_2024_table.txt',sep='\t').iloc[:11,:]
fig,ax = plt.subplots()
ax.barh(y=np.arange(result.shape[0]),width=np.flip(np.negative(np.log10(result['Adjusted P-value'].values))))
ax.set_xlabel('-log10(adjusted p-value)')
ax.set_yticks(np.arange(result.shape[0]))
ax.set_yticklabels(np.flip(result['Term'].values),fontsize=4)
plt.savefig('plasma_enrichr.pdf',bbox_inches='tight')
plt.close()
sys.exit('stop')

# use peptide abudance to cluster immunopeptidome samples
rootdir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas'
c_df_list = []
for c in cancers:
    print(c)
    final_all = pd.read_csv(os.path.join(rootdir,c,'antigen','0.01','final_enhanced_all.txt'),sep='\t')
    cond = [False if ('[]' in item) and ('(\'HLA-' not in item) else True for item in final_all['presented_by_each_sample_hla']]
    final_all = final_all.loc[cond,:]
    final_all = final_all.loc[final_all['typ']=='self_gene',:]
    final_all = final_all.loc[final_all['unique'],:]
    gene2data = {}
    for ensg,sub_df in tqdm(final_all.groupby(by='ensgs')):
        sample2data = {}
        for item in sub_df['detailed_intensity']:
            item = literal_eval(item)
            for tup in item:
                sample2data.setdefault(tup[0],[]).append(tup[1])
        sample2value = {k:np.median(v) for k,v in sample2data.items()}   # sample1:0.56
        gene2data[ensg] = sample2value
    c_df = pd.DataFrame.from_dict(gene2data,orient='columns')
    c_df_list.append(c_df)

# add serum 1
final_add = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/plasma/plasma_gbm_PXD008127/antigen/other_alg/final_enhanced.txt',sep='\t')
cond = [False if ('[]' in item) and ('(\'HLA-' not in item) else True for item in final_add['presented_by_each_sample_hla']]
final_add = final_add.loc[cond,:]
final_add = final_add.loc[final_add['typ']=='self_gene',:]
final_add = final_add.loc[final_add['unique'],:]
gene2data = {}
for ensg,sub_df in tqdm(final_add.groupby(by='ensgs')):
    sample2data = {}
    for item in sub_df['detailed_intensity']:
        item = literal_eval(item)
        for tup in item:
            sample2data.setdefault(tup[0],[]).append(tup[1])
    sample2value = {k:np.median(v) for k,v in sample2data.items()}   # sample1:0.56
    gene2data[ensg] = sample2value
c_df = pd.DataFrame.from_dict(gene2data,orient='columns')  # sample * gene
c_df.to_csv()
c_df_list.append(c_df)

# add serum 2
final_add = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/plasma/lung_cancer_nyu/antigen/other_alg/final_enhanced.txt',sep='\t')
cond = [False if ('[]' in item) and ('(\'HLA-' not in item) else True for item in final_add['presented_by_each_sample_hla']]
final_add = final_add.loc[cond,:]
final_add = final_add.loc[final_add['typ']=='self_gene',:]
final_add = final_add.loc[final_add['unique'],:]
gene2data = {}
for ensg,sub_df in tqdm(final_add.groupby(by='ensgs')):
    sample2data = {}
    for item in sub_df['detailed_intensity']:
        item = literal_eval(item)
        for tup in item:
            sample2data.setdefault(tup[0],[]).append(tup[1])
    sample2value = {k:np.median(v) for k,v in sample2data.items()}   # sample1:0.56
    gene2data[ensg] = sample2value
c_df = pd.DataFrame.from_dict(gene2data,orient='columns')  # sample * gene
c_df_list.append(c_df)

# add serum 3
final_add = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/plasma/lung_cancer_plasma_public/antigen/other_alg/final_enhanced.txt',sep='\t')
cond = [False if ('[]' in item) and ('(\'HLA-' not in item) else True for item in final_add['presented_by_each_sample_hla']]
final_add = final_add.loc[cond,:]
final_add = final_add.loc[final_add['typ']=='self_gene',:]
final_add = final_add.loc[final_add['unique'],:]
gene2data = {}
for ensg,sub_df in tqdm(final_add.groupby(by='ensgs')):
    sample2data = {}
    for item in sub_df['detailed_intensity']:
        item = literal_eval(item)
        for tup in item:
            sample2data.setdefault(tup[0],[]).append(tup[1])
    sample2value = {k:np.median(v) for k,v in sample2data.items()}   # sample1:0.56
    gene2data[ensg] = sample2value
c_df = pd.DataFrame.from_dict(gene2data,orient='columns')  # sample * gene
c_df_list.append(c_df)


df = pd.concat(c_df_list,axis=0,join='outer',keys=cancers+['gbm_public','luad_nyu','luad_public']).fillna(value=0).T # gene * sample

mi = df.columns
mi_df = mi.to_frame(index=False)
df.columns = mi_df[1].values


gene_lfc = pd.read_csv(os.path.join(rootdir,'GBM','gene_lfc.txt'),sep='\t',index_col=0)
ensg2symbol = gene_lfc['gene_symbol'].to_dict()
df = df.T # sample * gene
df.to_csv('umap_cancer_immuno_data.txt',sep='\t')
df = pd.read_csv('umap_cancer_immuno_data.txt',sep='\t',index_col=0)
adata = ad.AnnData(X=df.values,obs=pd.DataFrame(index=df.index,data={'cancer':mi_df[0].values}),var=pd.DataFrame(index=df.columns))  
sc.pp.highly_variable_genes(adata,flavor='seurat',n_top_genes=5000)
adata.raw = adata
adata = adata[:,adata.var['highly_variable']]
sc.pp.scale(adata,max_value=10)
sc.tl.pca(adata,n_comps=50)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
umap_dual_view_save(adata,cols=['cancer'])
adata = adata.raw.to_adata()
adata_ori = adata.copy()
adata = adata_ori[adata_ori.obs['cancer'].isin(['GBM','LUAD','gbm_public','luad_nyu','luad_public'])]

gene_to_ensg = {
    "ETV5":   "ENSG00000244405",
    "SOX4":   "ENSG00000124766",
    "MDM2":   "ENSG00000135679",
    "DDR1":   "ENSG00000204580",
    "CSF1":   "ENSG00000184371",
    "UBE2A":  "ENSG00000077721",
    "BST2":   "ENSG00000130303",
    "NPM1":   "ENSG00000181163",
    "BCAP31": "ENSG00000185825",
    "CDKN2A": "ENSG00000147889",
    "SART3":  "ENSG00000075856",
    "PA2G4":  "ENSG00000170515",
    "COTL1":  "ENSG00000103187",
    "ATIC":   "ENSG00000138363",
    "GPNMB":  "ENSG00000136235",
    "CDK4":   "ENSG00000135446",
    "SOX11":  "ENSG00000176887",
}
ensgs = list(set(gene_to_ensg.values()).intersection(set(adata.var_names)))
sc.pl.heatmap(adata,ensgs,groupby='cancer',swap_axes=True,dendrogram=True)
plt.savefig('umap_heatmap_taa.pdf',bbox_inches='tight')
plt.close()
sc.pl.violin(adata,keys=ensgs,groupby='cancer')
plt.savefig('umap_violin_taa.pdf',bbox_inches='tight')
plt.close()

gene_to_ensg = {
    # GBM-enriched
    "OLIG2":  "ENSG00000205927",
    "GFAP":   "ENSG00000131095",
    "SOX11":  "ENSG00000176887",
    "AQP4":   "ENSG00000171885",

    # LUAD-enriched
    "NKX2-1": "ENSG00000136352",
    "NAPSA":  "ENSG00000131400",
    "SFTPB":  "ENSG00000168878",
    "SFTPC":  "ENSG00000168484",
}
ensgs = list(set(gene_to_ensg.values()).intersection(set(adata.var_names)))
sc.pl.heatmap(adata,ensgs,groupby='cancer',swap_axes=True,dendrogram=True)
plt.savefig('umap_heatmap_enrich.pdf',bbox_inches='tight')
plt.close()
sc.pl.violin(adata,keys=ensgs,groupby='cancer')
plt.savefig('umap_violin_enrich.pdf',bbox_inches='tight')
plt.close()

adata = adata_ori[adata_ori.obs['cancer']!='luad_public',:]
sc.tl.rank_genes_groups(adata,'cancer')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
plt.savefig('umap_genes.pdf',bbox_inches='tight')
sys.exit('stop')



# use transcriptome to cluster TCGA
rootdir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas'
genes = set()
for c in cancers:
    final_all = pd.read_csv(os.path.join(rootdir,c,'antigen','0.01','final_enhanced_all.txt'),sep='\t')
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





