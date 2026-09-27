#!/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/bin/python3.7

import os,sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import re
import subprocess
from tqdm import tqdm
from io import StringIO
import multiprocessing as mp
import math
import bisect
from matplotlib_venn import venn2,venn3

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

cancers = [
    'BRCA','KIRC','COAD','STAD','MESO','LIHC','ESCA','CESC','BLCA','RT','AML','DLBC','GBM','NBL','PAAD','HNSC','OV','LUSC','LUAD','CHOL','SKCM'
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

def categorize_sets(set1, set2, set3):
    return [
        set1 - set2 - set3,  # Unique to set1
        set2 - set1 - set3,  # Unique to set2
        (set1 & set2) - set3,  # Overlap between set1 and set2 but not in set3
        set3 - set1 - set2,  # Unique to set3
        (set3 & set1) - set2,  # Overlap between set3 and set1 but not in set2
        (set3 & set2) - set1,  # Overlap between set3 and set2 but not in set1
        set1 & set2 & set3  # Overlap among all three
    ]

def compare_result(tesorai,df):
    col = []
    for item1,item2 in zip(tesorai['filename'],tesorai['scan_id']):
        col.append(item1.split('.')[0] + ',' + str(item2))
    tesorai['uid'] = col
    tesorai.set_index(keys='uid',inplace=True)

    uid2sequence = tesorai['clean_sequence'].to_dict()
    uid2rt = tesorai['retention_time'].to_dict()
    uid2pre_mz = tesorai['precursor_mz'].to_dict()
    uid2pre_intensity = tesorai['intensity'].to_dict()
    uid2pre_scan = tesorai['precursor_scan_id'].to_dict()
    uid2proteins = tesorai['possible_protein_ids'].to_dict()
    uid2score = tesorai['score'].to_dict()

    df['tesorai_sequence'] = df['uid'].map(uid2sequence).values
    df['tesorai_rt'] = df['uid'].map(uid2rt).values
    df['tesorai_pre_mz'] = df['uid'].map(uid2pre_mz).values
    df['tesorai_pre_intensity'] = df['uid'].map(uid2pre_intensity).values
    df['tesorai_pre_scan'] = df['uid'].map(uid2pre_scan).values
    df['tesorai_proteins'] = df['uid'].map(uid2proteins).values
    df['tesorai_score'] = df['uid'].map(uid2score).values

    df['same_seq'] = [i1==i2 for i1,i2 in zip(df['Sequence'],df['tesorai_sequence'])]
    df['tesorai_length'] = [len(item) if isinstance(item,str) else None for item in df['tesorai_sequence']]

    selected_columns = [
        'Raw file',
        'Scan number',
        'Retention time',
        'Identified',
        'Sequence',
        'm/z',
        'Precursor full scan number',
        'Precursor intensity',
        'Proteins',
        'Score',
        'PEP',
        'uid',
        'Identified_vanilla',
        'Identified_rescore',
        'tesorai_sequence',
        'tesorai_rt',
        'tesorai_pre_mz',
        'tesorai_pre_intensity',
        'tesorai_pre_scan',
        'tesorai_proteins',
        'tesorai_score',
        'same_seq',
        'tesorai_length'
    ]

    df = df.loc[:,selected_columns]
    return df



# data_all = []
# data_hq = []
# for c in cancers:
#     print(c)
#     # get all txt
#     immuno_dir = os.path.join('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/immunopeptidome',cancers2immuno[c])
#     cmd1 = 'find . -type f -name "msmsScans_new.txt"'
#     cmd2 = 'find . -type f -name "accumulatedMsmsScans_new.txt"'
#     old_dir = os.getcwd()
#     os.chdir(immuno_dir)
#     if c != 'MESO':
#         all_txt = subprocess.run(cmd1,shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#     else:
#         all_txt = subprocess.run(cmd2,shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#     all_study = [item.split('/')[1] for item in all_txt]
#     all_txt = [os.path.join(immuno_dir,item[2:]) for item in all_txt]
#     os.chdir(old_dir)

#     # combine
#     df_list = []
#     for txt,study in zip(all_txt,all_study):
#         msms = pd.read_csv(txt,sep='\t')
#         df_list.append(msms)
#     df = pd.concat(df_list,axis=0,keys=all_study).reset_index(level=-2)
#     # assert df.shape[0] == len(set(df['uid']))

#     # remove rev, some are maxquant-labelled rev, some has REV in proteins
#     col1 = []
#     col2 = []
#     col3 = []
#     for i1,i2,i3,i4 in zip(df['Identified'],df['Identified_vanilla'],df['Identified_rescore'],df['Proteins']):
#         if 'REV' not in i4:  
#             col1.append(i1)
#             col2.append(i2)
#             col3.append(i3)
#         else:
#             col1.append(None)
#             col2.append(False)
#             col3.append(False)
#     df['Identified'] = col1
#     df['Identified_vanilla'] = col2
#     df['Identified_rescore'] = col3

#     # compare
#     tesorai = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/NYU_Tesorai_all_searches/tesorai_peptide_fdr_{}.tsv'.format(c),sep='\t')
#     if c != 'MESO':
#         now_df = compare_result(tesorai,df)

    
#     # plot peptide 8-11 and additional 12-15 and ambiguity
#     fig = plt.figure(figsize=(15,6))
#     gs = mpl.gridspec.GridSpec(nrows=1,ncols=2,width_ratios=(0.5,0.5),wspace=0.2)
#     ax1 = fig.add_subplot(gs[0])
#     ax2 = fig.add_subplot(gs[1])

#     # particular for MESO
#     if c == 'MESO':
#         maxquant_peptide = set(df.loc[(df['Identified_vanilla']),:]['Sequence'])
#         rescore_peptide = set(df.loc[(df['Identified_rescore']==True),:]['Sequence'])
#         tesorai['tesorai_length'] = [len(item) for item in tesorai['clean_sequence']]
#         tesorai_peptide_8_11 = set(tesorai.loc[tesorai['tesorai_length'].isin([8,9,10,11]),:]['clean_sequence'])
#         tesorai_peptide_12_15 = set(tesorai.loc[tesorai['tesorai_length'].isin([12,13,14,15]),:]['clean_sequence'])
#         item7 = categorize_sets(maxquant_peptide,rescore_peptide,tesorai_peptide_8_11)
#         length7 = [len(item) for item in item7]
#         data_all.append([c]+length7)
#         venn3(subsets = tuple(length7),set_labels=('maxquant','rescore','tesorai_8_11'),ax=ax1)
#         ax1.set_title('{}\n#12-15:{}\n#ambiguity_psm:{},{}'.format(c,len(tesorai_peptide_12_15),0,0))

#         maxquant_peptide = set(df.loc[(df['Identified_vanilla']) & (df['Score']>70),:]['Sequence'])
#         rescore_peptide = set(df.loc[(df['Identified_rescore']==True) & (df['Score']>70),:]['Sequence'])
#         tesorai['tesorai_length'] = [len(item) for item in tesorai['clean_sequence']]
#         tesorai_peptide_8_11 = set(tesorai.loc[(tesorai['tesorai_length'].isin([8,9,10,11])) & (tesorai['score']>5),:]['clean_sequence'])
#         tesorai_peptide_12_15 = set(tesorai.loc[(tesorai['tesorai_length'].isin([12,13,14,15])) & (tesorai['score']>5),:]['clean_sequence'])
#         item7 = categorize_sets(maxquant_peptide,rescore_peptide,tesorai_peptide_8_11)
#         length7 = [len(item) for item in item7]
#         data_hq.append([c]+length7)
#         venn3(subsets = tuple(length7),set_labels=('maxquant','rescore','tesorai_8_11'),ax=ax2)
#         ax2.set_title('{}\n#12-15:{}\n#ambiguity_psm:{},{}'.format(c,len(tesorai_peptide_12_15),0,0))

#     else:
#         maxquant_peptide = set(now_df.loc[(now_df['Identified_vanilla']),:]['Sequence'])
#         rescore_peptide = set(now_df.loc[(now_df['Identified_rescore']==True),:]['Sequence'])
#         tesorai_peptide_8_11 = set(now_df.loc[now_df['tesorai_length'].isin([8,9,10,11]),:]['tesorai_sequence'])
#         tesorai_peptide_12_15 = set(now_df.loc[now_df['tesorai_length'].isin([12,13,14,15]),:]['tesorai_sequence'])
#         ambiguity = now_df.loc[(now_df['Identified']=='+') & (now_df['tesorai_sequence'].notna()) & (~now_df['same_seq']),:]
#         ambiguity1 = ambiguity.loc[(ambiguity['Score']>70) & (ambiguity['tesorai_score']>5),:]
#         ambiguity2 = ambiguity.loc[(ambiguity['Score']<70) & (ambiguity['tesorai_score']<5),:]
#         item7 = categorize_sets(maxquant_peptide,rescore_peptide,tesorai_peptide_8_11)
#         length7 = [len(item) for item in item7]
#         data_all.append([c]+length7)
#         venn3(subsets = tuple(length7),set_labels=('maxquant','rescore','tesorai_8_11'),ax=ax1)
#         ax1.set_title('{}\n#12-15:{}\n#ambiguity_psm:{},{}'.format(c,len(tesorai_peptide_12_15),ambiguity1.shape[0],ambiguity2.shape[0]))

#         maxquant_peptide = set(now_df.loc[(now_df['Identified_vanilla']) & (now_df['Score']>70),:]['Sequence'])
#         rescore_peptide = set(now_df.loc[(now_df['Identified_rescore']==True) & (now_df['Score']>70),:]['Sequence'])
#         tesorai_peptide_8_11 = set(now_df.loc[(now_df['tesorai_length'].isin([8,9,10,11])) & (now_df['tesorai_score']>5),:]['tesorai_sequence'])
#         tesorai_peptide_12_15 = set(now_df.loc[(now_df['tesorai_length'].isin([12,13,14,15])) & (now_df['tesorai_score']>5),:]['tesorai_sequence'])
#         item7 = categorize_sets(maxquant_peptide,rescore_peptide,tesorai_peptide_8_11)
#         length7 = [len(item) for item in item7]
#         data_hq.append([c]+length7)
#         venn3(subsets = tuple(length7),set_labels=('maxquant','rescore','tesorai_8_11'),ax=ax2)
#         ax2.set_title('{}\n#12-15:{}\n#ambiguity_psm:{},{}'.format(c,len(tesorai_peptide_12_15),ambiguity1.shape[0],ambiguity2.shape[0]))

#     plt.savefig('venn3_plot_{}.pdf'.format(c),bbox_inches='tight')
#     plt.close()

# final_data_all = pd.DataFrame.from_records(data_all,columns=['cancer','unique_to_maxquant','unique_to_rescore','m_r_not_t','unique_to_tesorai','t_m_not_r','t_r_not_m','all'])
# final_data_all.to_csv('final_data_all.txt',sep='\t',index=None)
# final_data_hq = pd.DataFrame.from_records(data_hq,columns=['cancer','unique_to_maxquant','unique_to_rescore','m_r_not_t','unique_to_tesorai','t_m_not_r','t_r_not_m','all'])
# final_data_hq.to_csv('final_data_hq.txt',sep='\t',index=None)


# plot all
final_data_all = pd.read_csv('final_data_all.txt',sep='\t',index_col=0)
a = final_data_all.iloc[:,[0,2,4,6]].sum(axis=1) 
b = final_data_all.iloc[:,[0,2,4,6,1,5]].sum(axis=1) 
c = final_data_all.iloc[:,[0,2,4,6,3,5]].sum(axis=1) 
to_plot = pd.concat([a,b,c],axis=1)
to_plot.columns = ['maxquant','rescore','tesorai']

# no break
to_plot.plot.bar(rot=90)
plt.savefig('final_data_all.pdf',bbox_inches='tight')
plt.close()

# break
fig,axes = plt.subplots(nrows=2,ncols=1,sharex=True,gridspec_kw={'height_ratios':(0.7,0.3),'hspace':0.1})
to_plot.plot.bar(rot=90,ax=axes[0])
to_plot.plot.bar(rot=90,ax=axes[1])
axes[0].set_ylim((10000,350000))
axes[1].set_ylim((0,10000))
axes[0].spines['bottom'].set_visible(False)
axes[1].spines['top'].set_visible(False)
axes[0].tick_params(bottom=False)
d = 0.015
axes[0].plot((-d,d),(-d,d),transform=axes[0].transAxes,clip_on=False,color='k')
axes[0].plot((1-d,1+d),(-d,d),transform=axes[0].transAxes,clip_on=False,color='k')
axes[1].plot((-d,d),(1-d,1+d),transform=axes[1].transAxes,clip_on=False,color='k')
axes[1].plot((1-d,1+d),(1-d,1+d),transform=axes[1].transAxes,clip_on=False,color='k')
plt.savefig('final_data_all_break.pdf',bbox_inches='tight')
plt.close()

# stats
to_plot['rescore_increase_by'] = (to_plot['rescore'] - to_plot['maxquant']) / to_plot['maxquant']
to_plot['tesorai_increase_by'] = (to_plot['tesorai'] - to_plot['maxquant']) / to_plot['maxquant']
to_plot.to_csv('increase_by_data_all.txt',sep='\t')


# plot hq
final_data_hq = pd.read_csv('final_data_hq.txt',sep='\t',index_col=0)
a = final_data_hq.iloc[:,[0,2,4,6]].sum(axis=1) 
b = final_data_hq.iloc[:,[0,2,4,6,1,5]].sum(axis=1) 
c = final_data_hq.iloc[:,[0,2,4,6,3,5]].sum(axis=1) 
to_plot = pd.concat([a,b,c],axis=1)
to_plot.columns = ['maxquant','rescore','tesorai']
to_plot.plot.bar(rot=90)
plt.savefig('final_data_hq.pdf',bbox_inches='tight')
plt.close()

fig,axes = plt.subplots(nrows=2,ncols=1,sharex=True,gridspec_kw={'height_ratios':(0.7,0.3),'hspace':0.1})
to_plot.plot.bar(rot=90,ax=axes[0])
to_plot.plot.bar(rot=90,ax=axes[1])
axes[0].set_ylim((2000,250000))
axes[1].set_ylim((0,2000))
axes[0].spines['bottom'].set_visible(False)
axes[1].spines['top'].set_visible(False)
axes[0].tick_params(bottom=False)
d = 0.015
axes[0].plot((-d,d),(-d,d),transform=axes[0].transAxes,clip_on=False,color='k')
axes[0].plot((1-d,1+d),(-d,d),transform=axes[0].transAxes,clip_on=False,color='k')
axes[1].plot((-d,d),(1-d,1+d),transform=axes[1].transAxes,clip_on=False,color='k')
axes[1].plot((1-d,1+d),(1-d,1+d),transform=axes[1].transAxes,clip_on=False,color='k')
plt.savefig('final_data_hq_break.pdf',bbox_inches='tight')
plt.close()
to_plot['rescore_increase_by'] = (to_plot['rescore'] - to_plot['maxquant']) / to_plot['maxquant']
to_plot['tesorai_increase_by'] = (to_plot['tesorai'] - to_plot['maxquant']) / to_plot['maxquant']
to_plot.to_csv('increase_by_data_hq.txt',sep='\t')



    


