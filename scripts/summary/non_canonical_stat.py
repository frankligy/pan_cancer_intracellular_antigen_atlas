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
bayests_xy_path = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/gene/full_results_XY_essential_tissues.txt'
gtex_median_path = '/gpfs/data/yarmarkovichlab/chordoma/NeoVerse_analysis/bulk-gex_v8_rna-seq_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct'
membrane_path = '/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/NeoVerse/data/human_membrane_proteins_acc2ens.txt'

safety_screen_df = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/code/post_safety_screen.txt',sep='\t')
safety_screen_df = safety_screen_df.loc[~safety_screen_df['cond_stringent'],:]
safety_screen_bl = list(set(safety_screen_df['pep'].values.tolist()))

splicing_df = pd.read_csv('all_splicing.txt',sep='\t')
splicing_df = splicing_df.loc[~splicing_df['pep'].isin(safety_screen_bl),:]
nuorf_df = pd.read_csv('all_nuorf.txt',sep='\t')
nuorf_df = nuorf_df.loc[~nuorf_df['pep'].isin(safety_screen_bl),:]
variant_df = pd.read_csv('all_variants_ts.txt',sep='\t')
fusion_df = pd.read_csv('all_fusion.txt',sep='\t')
fusion_df = fusion_df.loc[~fusion_df['pep'].isin(safety_screen_bl),:]
fusion_df = fusion_df.loc[~fusion_df['source'].str.contains('nc'),:]
ir_df = pd.read_csv('all_ir.txt',sep='\t')
ir_df = ir_df.loc[~ir_df['pep'].isin(safety_screen_bl),:]

self_translate_te_df = pd.read_csv('ts_te_antigen.txt',sep='\t')  # remember, after safety screen, do autonomy check and update 
real_autonomy_check = pd.read_csv('splicing_ir_dic/final.txt',sep='\t')
real_autonomy = set(real_autonomy_check.loc[real_autonomy_check['not_has_ss'] & real_autonomy_check['not_in_ir'],:]['pep'].values)
self_translate_te_df = self_translate_te_df.loc[self_translate_te_df['pep'].isin(real_autonomy),:]
te_all_df = pd.read_csv('te_all_antigens.txt',sep='\t')
orf2_taa = te_all_df.loc[te_all_df['source'].str.contains('L1_ORF2'),:]
self_translate_te_df = pd.concat([self_translate_te_df,orf2_taa],axis=0)
self_translate_te_df['typ'] = np.full(shape=self_translate_te_df.shape[0],fill_value='self_translate_te')
self_translate_te_df = self_translate_te_df.loc[~self_translate_te_df['pep'].isin(safety_screen_bl),:]

te_chimeric_df = pd.read_csv('te_all_antigens.txt',sep='\t')
te_chimeric_df = te_chimeric_df.loc[te_chimeric_df['typ']=='TE_chimeric_transcript',:]
original_all_self_translate =  pd.read_csv('ts_te_antigen.txt',sep='\t')
reclassified_te_chimeric = original_all_self_translate.loc[~original_all_self_translate['pep'].isin(real_autonomy),:]
te_chimeric_df = pd.concat([te_chimeric_df,reclassified_te_chimeric],axis=0)
te_chimeric_df['typ'] = np.full(shape=te_chimeric_df.shape[0],fill_value='TE_chimeric_transcript')
te_chimeric_df = te_chimeric_df.loc[~te_chimeric_df['pep'].isin(safety_screen_bl),:]


data = np.empty((21,7),dtype=np.float64)
for i,df in enumerate([fusion_df,variant_df,ir_df,splicing_df,nuorf_df,self_translate_te_df,te_chimeric_df]):
    tmp = []
    all_c = df['cancer'].unique().tolist()
    for c in cancers:
        if c not in all_c:
            tmp.append(0)
        else:
            tmp.append(df.loc[df['cancer']==c,:].shape[0])
    data[:,i] = np.array(tmp)

data = data / data.sum(axis=1).reshape(-1,1)

df = pd.DataFrame(data=data,index=cancers,columns=['fusion','variant','intron','splicing','nuORF','self_translating_TE','te_chimeric'])
df.to_csv('non_canonical_stat_df.txt',sep='\t')


fig,ax = plt.subplots()
ax.bar(cancers,data[:,0],label='fusion')
ax.bar(cancers,data[:,1],bottom=data[:,0],label='variant')
ax.bar(cancers,data[:,2],bottom=data[:,:2].sum(axis=1),label='intron')
ax.bar(cancers,data[:,3],bottom=data[:,:3].sum(axis=1),label='splicing')
ax.bar(cancers,data[:,4],bottom=data[:,:4].sum(axis=1),label='nuORF')
ax.bar(cancers,data[:,5],bottom=data[:,:5].sum(axis=1),label='self_translating')
ax.bar(cancers,data[:,6],bottom=data[:,:6].sum(axis=1),label='te_chimeric')
ax.set_ylabel('Proportion')
ax.set_ylim([0,1.05])
ax.set_xticklabels(cancers,rotation=60)
ax.legend()
plt.savefig('non_canonical_prop.pdf',bbox_inches='tight')
plt.close()






    