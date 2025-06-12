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
from matplotlib_venn import venn2,venn3

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

root_atlas_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas'
driver_genes = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/variants/driver_genes.txt'

def categorize_sets(set1, set2):
    return [
        set1 - set2,  # Unique to set1
        set2 - set1,  # Unique to set2
        set1 & set2   # Overlap between set1 and set2 
    ]

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
col = []
for item in final['source']:
    if ';' not in item:
        col.append(item)
    else:
        add = None
        for s in item.split(';'):
            if ('variant' in s) or ('inframe' in s):
                add = s
                break
        col.append(add)
final['real_source'] = col
final['gene'] = [item.split('|')[0] for item in final['real_source']]
final['is_driver'] = final['gene'].isin(driver).values
final['mutation'] = ['|'.join([item.split('|')[0],item.split('|')[1]]) for item in final['real_source']]
final['mutation_type'] = [item.split('|')[-1] for item in final['real_source']]
final['recurrency'] = [item.split('|')[2] for item in final['real_source']]
final.to_csv('all_variants.txt',sep='\t',index=None)


df = pd.read_csv('final_all_ts_antigens.txt',sep='\t')
final = df.loc[df['typ']=='variant',:]


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

# plot mutation type
col1 = []
col2 = []
for mt,sub_df in final.groupby(by='mutation_type'):
    col1.append(mt)
    col2.append(sub_df.shape[0])
df = pd.DataFrame(data={'mutation_type':col1,'value':col2}).sort_values(by='value',ascending=True)
df.plot.barh(x='mutation_type',y='value',rot=0,fontsize=5)
plt.savefig('variant_by_mutation_type.pdf',bbox_inches='tight')
plt.close()


# compare with common database
pep_iedb = pd.read_csv('all_epitope_no_b_human_linear_mhc_i.tsv',sep='\t')['Epitope - Name'].values.tolist()
pep_tsnadb_snv = pd.read_csv('SNV-derived.txt',sep='\t')
pep_tsnadb_indel = pd.read_csv('INDEL-derived.txt',sep='\t')
pep_tsnadb_shared = pd.read_csv('tsnadb_shared_20.txt',sep='\t')['Peptide'].values.tolist()

# pep_tsnadb_snv = pep_tsnadb_snv.loc[~pep_tsnadb_snv['Tissue'].isin(['Prostate','Thyroid','Uterus']),:]
# pep_tsnadb_indel = pep_tsnadb_indel.loc[~pep_tsnadb_indel['Tissue'].isin(['Prostate','Thyroid','Uterus']),:]
pep_tsnadb = pep_tsnadb_snv['Peptide'].values.tolist() + pep_tsnadb_indel['Peptide'].values.tolist()

common1 = set(final['pep'].values).intersection(set(pep_iedb))
common2 = set(final['pep'].values).intersection(set(pep_tsnadb))
common3 = set(final['pep'].values).intersection(set(pep_tsnadb_shared))

final['in_iedb'] = [True if item in common1 else False for item in final['pep']]
final['in_tsnadb'] = [True if item in common2 else False for item in final['pep']]

tesorai_dict = {}
for c in cancers:
    tesorai_dict[c] = set(pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/NYU_Tesorai_all_searches/tesorai_peptide_fdr_{}.tsv'.format(c),sep='\t')['clean_sequence'].values.tolist())
cond = []
for c,p in zip(final['cancer'],final['pep']):
    if p in tesorai_dict[c]:
        cond.append(True)
    else:
        cond.append(False)
final['in_tesorai'] = cond
final.to_csv('annotated_variant_antigen.txt',sep='\t',index=None)

item3 = categorize_sets(set(pep_tsnadb),set(final['pep'].values))
length3 = [len(item) for item in item3]
venn2(subsets=length3,set_labels=('TSNAdb','ImmunoVerse'))
plt.savefig('immunoverse_tsnadb.pdf',bbox_inches='tight')
plt.close()

