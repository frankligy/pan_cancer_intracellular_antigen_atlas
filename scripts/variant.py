#!/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/bin/python3.7

import pandas as pd
import numpy as np
import sys,os
from Bio.SeqIO.FastaIO import SimpleFastaParser
import subprocess
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from ast import literal_eval
from tqdm import tqdm
import pickle
import anndata as ad
from scipy.sparse import csr_matrix

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
    'SKCM'
]  

n_samples = [
    1118,
    542,
    483,
    412,
    87,
    374,
    412,
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
    541,
    502,
    35,
    472
]

root_atlas_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas'
driver_genes = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/variants/driver_genes.txt'


driver = set(pd.read_csv(driver_genes,sep='\t')['gene'])
data = []
for c in cancers:
    final_path = os.path.join(root_atlas_dir,c,'antigen','fdr','final_enhanced.txt')
    final = pd.read_csv(final_path,sep='\t')
    cond = [False if ('[]' in item) and ('(\'HLA-' not in item) else True for item in final['presented_by_each_sample_hla']]
    final = final.loc[cond,:]
    final = final.loc[final['typ']=='variant',:]
    final = final.loc[final['unique'],:]
    data.append(final)
final = pd.concat(data,axis=0,keys=cancers).reset_index(level=-2).rename(columns={'level_0':'cancer'})
final['gene'] = [item.split('|')[0] for item in final['source']]
final['is_driver'] = final['gene'].isin(driver).values
final['mutation'] = ['|'.join([item.split('|')[0],item.split('|')[1]]) for item in final['source']]
final.to_csv('all_variants.txt',sep='\t',index=None);sys.exit('stop')
# final['pep'].value_counts().to_csv('all_variant_peptide.txt',sep='\t')
# final['gene'].value_counts().to_csv('all_variant_gene.txt',sep='\t')
# final['mutation'].value_counts().to_csv('all_variant_mutation.txt',sep='\t')

# plot cancer
col1 = []
col2 = []
for c,sub_df in final.groupby(by='cancer'):
    col1.append(c)
    col2.append(sub_df.shape[0])
other_c = list(set(cancers).difference(set(col1)))
for c in other_c:
    col1.append(c)
    col2.append(0)
df = pd.DataFrame(data={'cancer':col1,'value':col2}).sort_values(by='value',ascending=True)
print(df);sys.exit('stop')
df.plot.barh(x='cancer', y='value', rot=0)
plt.savefig('variant_by_cancer.pdf',bbox_inches='tight')
plt.close()

# plot gene
col1 = []
col2 = []
for g,sub_df in final.groupby(by='gene'):
    col1.append(g)
    col2.append(sub_df.shape[0])
df = pd.DataFrame(data={'gene':col1,'value':col2}).sort_values(by='value',ascending=True)
df.plot.barh(x='gene', y='value', rot=0,fontsize=1)
plt.savefig('variant_by_gene.pdf',bbox_inches='tight')
plt.close()

# plot gene driver
col1 = []
col2 = []
final_sub = final.loc[final['is_driver'],:]
for g,sub_df in final_sub.groupby(by='gene'):
    col1.append(g)
    col2.append(sub_df.shape[0])
df = pd.DataFrame(data={'gene':col1,'value':col2}).sort_values(by='value',ascending=True)
df.plot.barh(x='gene', y='value', rot=0,fontsize=5)
plt.savefig('variant_by_gene_driver.pdf',bbox_inches='tight')
plt.close()





