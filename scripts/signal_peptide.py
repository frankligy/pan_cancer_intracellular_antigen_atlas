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
human_canonical_proteome = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/gene/ensembl_protein.fasta'
MEMBRANE_GENE = '/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/NeoVerse/data/human_membrane_proteins_acc2ens.txt'


df = pd.read_csv(MEMBRANE_GENE,sep='\t')
membrane_ensg = set(df['Ens'].tolist())

with open(human_canonical_proteome,'r') as in_handle, open('membrane_protein.fasta','w') as f:
    for title,seq in SimpleFastaParser(in_handle):
        ensg,enst,gs = title.split('|')
        if ensg in membrane_ensg:
            f.write('>{}\n{}\n'.format(title,seq))


dic = {}
with open('membrane_protein.fasta','r') as in_handle:
    for title,seq in SimpleFastaParser(in_handle):
        dic[title] = seq


final = {}
ensg2symbol = {}
final2prop = {}
result = pd.read_csv('run_summary_sec_spi.signalp5',sep='\t',index_col=0)
for uid,item in zip(result.index,result['CS Position']):
    tmp = item.split('CS pos: ')[1].split('. ')[0]
    p1,p2 = tmp.split('-')
    p1,p2 = int(p1),int(p2)
    k = uid.split('|')[0]
    gs = uid.split('|')[2]
    signal_seq = dic[uid][:p1]
    prop = (p1,len(dic[uid])-p1)
    ensg2symbol[k] = gs
    final[k] = signal_seq
    final2prop[k] = prop


# in signal
data = []
normal = pd.read_csv('/gpfs/data/yarmarkovichlab/Aman/hla_ligand_atlas.txt',sep='\t')
for item1,item2 in zip(normal['peptide'],normal['ensg']):
    signal_seq = final.get(item2,None)
    if signal_seq is not None:
        if item1 in signal_seq:
            data.append([item1,item2,ensg2symbol[item2],'normal'])


cancer = pd.read_csv('all_unique_self_gene.txt',sep='\t')
for item1,item2,item3 in zip(cancer['pep'],cancer['ensgs'],cancer['cancer']):
    signal_seq = final.get(item2,None)
    if signal_seq is not None:
        if item1 in signal_seq:
            data.append([item1,item2,ensg2symbol[item2],item3])

df = pd.DataFrame(data=data,index=np.arange(len(data)),columns=['pep','ensg','gs','source'])
df = df.drop_duplicates()
df.to_csv('all_signal_peptides.txt',sep='\t',index=None)
df_in_signal = df

# not in signal
data = []
normal = pd.read_csv('/gpfs/data/yarmarkovichlab/Aman/hla_ligand_atlas.txt',sep='\t')
for item1,item2 in zip(normal['peptide'],normal['ensg']):
    signal_seq = final.get(item2,None)
    if signal_seq is not None:
        if item1 not in signal_seq:
            data.append([item1,item2,ensg2symbol[item2],'normal'])


cancer = pd.read_csv('all_unique_self_gene.txt',sep='\t')
for item1,item2,item3 in zip(cancer['pep'],cancer['ensgs'],cancer['cancer']):
    signal_seq = final.get(item2,None)
    if signal_seq is not None:
        if item1 not in signal_seq:
            data.append([item1,item2,ensg2symbol[item2],item3])

df = pd.DataFrame(data=data,index=np.arange(len(data)),columns=['pep','ensg','gs','source'])
df = df.drop_duplicates()
df.to_csv('all_signal_peptides_not.txt',sep='\t',index=None)
df_not_in_signal = df

# calculate potential
metrics = []
for k,v in final2prop.items():
    n_in_signal = len(df_in_signal.loc[df_in_signal['ensg']==k,:]['pep'].unique())
    n_not_in_signal = len(df_not_in_signal.loc[df_not_in_signal['ensg']==k,:]['pep'].unique())
    pot_in_singal = n_in_signal / v[0]
    pot_not_in_signal = n_not_in_signal / v[1]
    metrics.append([n_in_signal,n_not_in_signal,pot_in_singal,pot_not_in_signal])
metrics_df = pd.DataFrame(data=metrics,index=[ensg2symbol[item] for item in final2prop.keys()],columns=['n_in_signal','n_not_in_signal','potential_in_signal','potential_not_in_signal'])
metrics_df.to_csv('signal_potential.txt',sep='\t')
metrics_df = metrics_df.loc[(metrics_df['n_in_signal']>0) | (metrics_df['n_not_in_signal']>0),:]
metrics_df.to_csv('signal_potential_useful.txt',sep='\t')


metrics_df.sort_values(by='potential_in_signal',ascending=False,inplace=True)
y1 = []
y2 = []
for i1,i2 in zip(metrics_df['potential_in_signal'],metrics_df['potential_not_in_signal']):
    y1.append(i1)
    y2.append(i2)
x = np.arange(len(y1))
fig,ax = plt.subplots()
ax.plot(x,y1,lw=1,marker='o',markersize=1,label='signal_region')
ax.plot(x,y2,lw=1,marker='o',markersize=1,label='non_signal_region')
ax.set_ylabel('presentation_potential')
ax.set_xlabel('membrane_proteins')
ax.set_xticks([])
ax.legend()
plt.savefig('signal_peptide_potential.pdf',bbox_inches='tight')
plt.close()







    



