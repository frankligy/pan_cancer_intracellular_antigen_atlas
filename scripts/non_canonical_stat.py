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
bayests_xy_path = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/gene/full_results_XY_essential_tissues.txt'
gtex_median_path = '/gpfs/data/yarmarkovichlab/chordoma/NeoVerse_analysis/bulk-gex_v8_rna-seq_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct'
membrane_path = '/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/NeoVerse/data/human_membrane_proteins_acc2ens.txt'



splicing_df = pd.read_csv('all_splicing.txt',sep='\t')
splicing_df = splicing_df.loc[splicing_df['highest_abundance'].notna(),:]
nuorf_df = pd.read_csv('all_nuorf.txt',sep='\t')
nuorf_df = nuorf_df.loc[nuorf_df['highest_abundance'].notna(),:]
nuorf_df = nuorf_df.loc[nuorf_df['highest_score']>40,:]
variant_df = pd.read_csv('all_variants.txt',sep='\t')
variant_df = variant_df.loc[variant_df['highest_abundance'].notna(),:]
fusion_df = pd.read_csv('all_fusion.txt',sep='\t')
fusion_df = fusion_df.loc[fusion_df['highest_abundance'].notna(),:]
fusion_df = fusion_df.loc[~fusion_df['source'].str.contains('nc'),:]
ir_df = pd.read_csv('all_ir.txt',sep='\t')
ir_df = ir_df.loc[ir_df['highest_abundance'].notna(),:]


data = np.empty((21,5),dtype=np.float64)
for i,df in enumerate([fusion_df,variant_df,ir_df,splicing_df,nuorf_df]):
    tmp = []
    all_c = df['cancer'].unique().tolist()
    for c in cancers:
        if c not in all_c:
            tmp.append(0)
        else:
            tmp.append(df.loc[df['cancer']==c,:].shape[0])
    data[:,i] = np.array(tmp)

data = data / data.sum(axis=1).reshape(-1,1)

df = pd.DataFrame(data=data,index=cancers,columns=['fusion','variant','intron','splicing','nuORF'])
df.to_csv('non_canonical_stat_df.txt',sep='\t')
sys.exit('stop')


fig,ax = plt.subplots()
ax.bar(cancers,data[:,0],label='fusion')
ax.bar(cancers,data[:,1],bottom=data[:,0],label='variant')
ax.bar(cancers,data[:,2],bottom=data[:,:2].sum(axis=1),label='intron')
ax.bar(cancers,data[:,3],bottom=data[:,:3].sum(axis=1),label='splicing')
ax.bar(cancers,data[:,4],bottom=data[:,:4].sum(axis=1),label='nuORF')
ax.set_ylabel('Proportion')
ax.set_ylim([0,1.05])
ax.set_xticklabels(cancers,rotation=60)
ax.legend()
plt.savefig('non_canonical_prop.pdf',bbox_inches='tight')
plt.close()






    