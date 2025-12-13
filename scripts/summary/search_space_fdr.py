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

cancers2immuno = {
    'BRCA':'breast_cancer',
    'KIRC':'kidney_clear_cell',
    'COAD':'colon_cancer',
    'STAD':'stomach_cancer',
    'MESO':'mesothelioma',
    'LIHC':'liver_cancer',
    'ESCA':'esophageal_cancer',
    'CESC':'cervical_cancer',
    'BLCA':'bladder_cancer',
    'RT':'rhabdoid_tumor',
    'AML':'AML',
    'DLBC':'DLBC',
    'GBM':'GBM',
    'NBL':'neuroblastoma',
    'PAAD':'pancreatic_cancer',
    'HNSC':'head_and_neck',
    'OV':'ovarian_cancer',
    'LUSC':'lung_LUSC',
    'LUAD':'lung_LUAD',
    'CHOL':'bile_duct_cancer',
    'SKCM':'melanoma'
}

def classify_source(item):
    if item.startswith('REV__'):
        identity = 'reverse'
    elif item.startswith('CON__'):
        identity = 'contaminant'
    elif 'missense_variant' in item or 'inframe' in item or 'frameshift_variant' in item:
        identity = 'variant'
    elif item.startswith('tr') or item.startswith('sp'):
        identity = 'pathogen'
    elif 'fusion' in item:
        identity = 'fusion' 
    elif 'intron' in item:
        identity = 'intron_retention'
    elif 'nuORF' in item:
        identity = 'nuORF'
    elif item.startswith('chr'):
        if 'TE_info' in item:
            identity = 'TE_chimeric_transcript'
        else:
            identity = 'splicing'
    elif '_dup' in item or item.endswith('sense'):
        identity = 'ERV'
    elif item.startswith('ENSG'):
        identity = 'self_gene'
    elif item.startswith('nc'):
        identity = 'nc_isoform'
    else:
        identity = 'unknown'
    return identity

root_atlas_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas'
immuno_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/immunopeptidome'

iedb = pd.read_csv('all_epitope_no_b_human_linear_mhc_i.tsv',sep='\t')
iedb_peptides = set(iedb['Epitope - Name'].unique())

cancers = ['BRCA','KIRC']

# search_space
cats = ['contaminant','variant','pathogen','fusion','intron_retention','nuORF','TE_chimeric_transcript','splicing','ERV','unknown','nc_isoform','self_gene']
data = np.empty((len(cancers),len(cats)),dtype=np.float64)
for i,c in enumerate(cancers):
    space_dict = {}
    combined_fasta = os.path.join(immuno_dir,cancers2immuno[c],'combined_{}_pan_cancer.fasta'.format(c))
    with open(combined_fasta,'r') as in_handle:
        for title,seq in SimpleFastaParser(in_handle):
            identity = classify_source(title)
            if identity in space_dict.keys():
                space_dict[identity] += len(seq)
            else:
                space_dict[identity] = 0
                space_dict[identity] += len(seq)

    for cat in cats:
        if cat not in space_dict.keys():
            space_dict[cat] = 0
    data[i,:] = pd.Series(space_dict).loc[cats].values

df = pd.DataFrame(data=data,index=cancers,columns=cats)
df['tmp'] = df['nc_isoform'] + df['self_gene']
df['self_gene'] = df['tmp']
df.drop(columns=['unknown','nc_isoform','tmp'],inplace=True)
increase_by = df.iloc[:,:9].sum(axis=1) / df.iloc[:,9]

data = df.values
fig,ax = plt.subplots()
ax.bar(cancers,data[:,0],label='contaminant')
ax.bar(cancers,data[:,1],bottom=data[:,0],label='variant')
ax.bar(cancers,data[:,2],bottom=data[:,:2].sum(axis=1),label='pathogen')
ax.bar(cancers,data[:,3],bottom=data[:,:3].sum(axis=1),label='fusion')
ax.bar(cancers,data[:,4],bottom=data[:,:4].sum(axis=1),label='intron_retention')
ax.bar(cancers,data[:,5],bottom=data[:,:5].sum(axis=1),label='nuORF')
ax.bar(cancers,data[:,6],bottom=data[:,:6].sum(axis=1),label='TE_chimeric_transcript')
ax.bar(cancers,data[:,7],bottom=data[:,:7].sum(axis=1),label='splicing')
ax.bar(cancers,data[:,8],bottom=data[:,:8].sum(axis=1),label='self_translating_TE')
ax.bar(cancers,data[:,9],bottom=data[:,:9].sum(axis=1),label='self_gene')
for i,c in enumerate(cancers):
    ax.text(x=i,y=data[:,:10].sum(axis=1)[i],s=round(increase_by.iloc[i],2))
ax.set_ylabel('value')
# ax.set_ylim([0,1.05])
ax.set_xticklabels(cancers,rotation=60)
ax.legend()
plt.savefig('search_space_total.pdf',bbox_inches='tight')
plt.close()

# 70 and iedb
# data_list = []
# for c in cancers:
#     msms1 = pd.read_csv(os.path.join(root_atlas_dir,c,'antigen','0.01','msmsScans_all_add_tesorai.txt'),sep='\t')
#     msms5 = pd.read_csv(os.path.join(root_atlas_dir,c,'antigen','0.05','msmsScans_all_add_tesorai.txt'),sep='\t')
#     # maxquant fdr iedb and 70
#     msms1_ = msms1.loc[msms1['Identified_vanilla']==True,:]
#     msms5_ = msms5.loc[msms5['Identified_vanilla']==True,:]
#     only_in_5 = set(msms5_['Sequence']).difference(set(msms1_['Sequence']))
#     iedb_common = iedb_peptides.intersection(only_in_5)
#     seventy1 = msms1_.loc[msms1_['Score']>70,:]
#     seventy5 = msms5_.loc[msms5_['Score']>70,:]
#     prop_seventy1 = seventy1.shape[0] / msms1_.shape[0]
#     prop_seventy5 = seventy5.shape[0] / msms5_.shape[0]
#     data_m = [len(iedb_common),prop_seventy1,prop_seventy5]

#     # rescore fdr iedb and 70
#     msms1_ = msms1.loc[msms1['Identified_rescore']==True,:]
#     msms5_ = msms5.loc[msms5['Identified_rescore']==True,:]
#     only_in_5 = set(msms5_['Sequence']).difference(set(msms1_['Sequence']))
#     iedb_common = iedb_peptides.intersection(only_in_5)
#     seventy1 = msms1_.loc[msms1_['Score']>70,:]
#     seventy5 = msms5_.loc[msms5_['Score']>70,:]
#     prop_seventy1 = seventy1.shape[0] / msms1_.shape[0]
#     prop_seventy5 = seventy5.shape[0] / msms5_.shape[0]
#     data_r = [len(iedb_common),prop_seventy1,prop_seventy5]

#     data_list.append(data_m+data_r)

# df = pd.DataFrame.from_records(data=data_list,index=cancers,columns=['maxquant_iedb_only_5','maxquant_hc_1','maxquant_hc_5','rescore_iedb_only_5','rescore_hc_1','rescore_hc_5'])
# df1 = df.loc[:,['maxquant_hc_1','maxquant_hc_5']]
# df1.plot.bar(rot=0)
# plt.savefig('maxquant_hc.pdf',bbox_inches='tight')
# plt.close()

# df1 = df.loc[:,['rescore_hc_1','rescore_hc_5']]
# df1.plot.bar(rot=0)
# plt.savefig('rescore_hc.pdf',bbox_inches='tight')
# plt.close()

# df1 = df.loc[:,['maxquant_iedb_only_5','rescore_iedb_only_5']]
# df1.plot.bar(rot=0)
# plt.savefig('iedb_miss.pdf',bbox_inches='tight')
# plt.close()

# anchor
data = []
anchors = {
    'KIRC':['ATFLGSLTGK', 'KLIAGLIFLK', 'DLSRRDVSL']
}

for k,vs in anchors.items():
    msms5 = pd.read_csv(os.path.join(root_atlas_dir,k,'antigen','0.05','msmsScans_all_add_tesorai.txt'),sep='\t')
    for v in vs:
        sub = msms5.loc[(msms5['Sequence']==v) & (msms5['Identified']=='+'),:]
        # maxquant range
        sub2 = sub.loc[sub['Identified_vanilla']==True,:]
        min_m = sub2['qval_vanilla'].min()
        max_m = sub2['qval_vanilla'].max()
        # rescore range
        sub2 = sub.loc[sub['Identified_rescore']==True,:]
        min_r = sub2['qval_rescore'].min()
        max_r = sub2['qval_rescore'].max()
        data.append((k,v,min_m,max_m,min_r,max_r))
df = pd.DataFrame(data=data,columns=['cancer','antigen','maxquant_qval_min','maxquant_qval_max','rescore_qval_min','rescore_qval_max'])
df.to_csv('anchor_qval_range.txt',sep='\t',index=None)



