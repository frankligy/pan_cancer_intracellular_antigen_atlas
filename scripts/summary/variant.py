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

def chop_normal_pep_db(fasta_path,output_path,mers,allow_duplicates):
    '''
    chop any normal human proteome to certain mers

    :param fasta_path: string, the path to the human protein fasta file
    :param output_path: string, the path to the output fasta file
    :param mers: list, like [9,10] will generate 9mer and 10mer
    :param allow_duplicates: boolean. whether allow duplicate or not

    Example::

        chop_normal_pep_db(fasta_path='human_uniprot_proteome.fasta',output_path='./human_uniprot_proteome_mer9_10.fasta',mers=[9,10],allow_duplicates=False)
    '''
    # for human genome in uniprot, 9-10mer, remove duplicates will decrease from 44,741,578 to 41,638,172
    if allow_duplicates:
        with open(fasta_path,'r') as in_handle, open(output_path,'w') as out_handle:
            for title,seq in tqdm(SimpleFastaParser(in_handle)):
                count = 0
                length = len(seq)
                for i in range(length):
                    for mer in mers:
                        if i+mer <= length:
                            out_handle.write('>{}_{}_{}\n'.format(title,mer,count))
                            out_handle.write('{}\n'.format(seq[i:i+mer]))
                            count += 1
    else:
        with open(fasta_path,'r') as in_handle, open(output_path,'w') as out_handle:
            existing = set()
            for title,seq in tqdm(SimpleFastaParser(in_handle)):
                count = 0
                length = len(seq)
                for i in range(length):
                    for mer in mers:
                        if i+mer <= length:
                            subseq = seq[i:i+mer]
                            if subseq not in existing:                                
                                out_handle.write('>{}_{}_{}\n'.format(title,mer,count))
                                out_handle.write('{}\n'.format(subseq))
                                existing.add(subseq)
                                count += 1  

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

final = pd.read_csv('./stats/final_all_ts_antigens.txt',sep='\t')
final = final.loc[final['typ']=='variant',:]
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

final['recurrency'] = final['recurrency'].astype('int')
tmp = final.loc[final['recurrency']>1,:]
print(len(tmp['pep'].unique()))
tmp = final.loc[final['recurrency']==1,:]
print(len(tmp['pep'].unique()))

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
df.to_csv('variant_by_cancer.txt',sep='\t',index=None)
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
df.to_csv('variant_by_gene.txt',sep='\t',index=None)
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
df.to_csv('variant_by_gene_driver.txt',sep='\t',index=None)
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
df.set_index(keys='mutation_type',inplace=True)
df.plot.pie(y='value',figsize=(5,5))
plt.savefig('variant_by_mutation_type_pie.pdf',bbox_inches='tight')
plt.close()
df['prop'] = df['value']/df['value'].sum()
df.to_csv('variant_by_mutation_type.txt',sep='\t')


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

tesorai_folder = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/NYU_Tesorai_all_searches'
tesorai_dict = {}
for c in cancers:
    # get all tesorai
    old_dir = os.getcwd()
    os.chdir(tesorai_folder)
    cmd = 'for f in tesorai_peptide_fdr_*.tsv; do echo $f; done | grep "{}"'.format(c)
    needed_files = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
    os.chdir(old_dir)

    # concat all seq
    t = []
    for tesorai_file in needed_files:
        tesorai = pd.read_csv(os.path.join(tesorai_folder,tesorai_file),sep='\t')
        a = tesorai['clean_sequence'].values.tolist()
        t.extend(a)
    tesorai_dict[c] = set(t)

cond = []
for c,p in zip(final['cancer'],final['pep']):
    if p in tesorai_dict[c]:
        cond.append(True)
    else:
        cond.append(False)
final['in_tesorai'] = cond
final.to_csv('annotated_variant_antigen.txt',sep='\t',index=None)
print(len(final['pep'].unique()))

tmp = final.loc[final['in_iedb'],:]
print(len(tmp['pep'].unique()))


# need to make sure pep_tsnadb are in the search space
titles = []
seqs = []
for c in cancers:
    p = os.path.join(root_atlas_dir,c,'db_fasta_tesorai','mutation.fasta')
    with open(p) as in_handle:
        for title,seq in SimpleFastaParser(in_handle):
            titles.append(title)
            seqs.append(seq)
with open('all_mutations_search_space.fasta','w') as f:
    for t,s in zip(titles,seqs):
        f.write('>{}\n{}\n'.format(t,s))

# for mer in [8,9,10,11,12,13,14,15]:
#     print(mer)
#     chop_normal_pep_db(fasta_path='all_mutations_search_space.fasta',
#                        output_path='all_mutations_search_space_{}.fasta'.format(mer),mers=[mer],allow_duplicates=False)


fasta_dic = {}
for mer in [8,9,10,11,12,13,14,15]:
    print(mer)
    lis = []
    with open('all_mutations_search_space_{}.fasta'.format(mer),'r') as in_handle:
        for title,seq in SimpleFastaParser(in_handle):
            lis.append(seq)
    lis = set(lis)
    fasta_dic[mer] = lis

valid_pep_tsnadb = []
for pep in pep_tsnadb:
    l = len(pep)
    lis = fasta_dic[l]
    if pep in lis:
        valid_pep_tsnadb.append(pep)


item3 = categorize_sets(set(valid_pep_tsnadb),set(final['pep'].values))
length3 = [len(item) for item in item3]
venn2(subsets=length3,set_labels=('TSNAdb','ImmunoVerse'))
plt.savefig('immunoverse_tsnadb.pdf',bbox_inches='tight')
plt.close()

