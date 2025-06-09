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
from scipy.stats import ttest_ind

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

rootdir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/immunopeptidome/melanoma'
peptides = {
    'PMEL_splicing':'KTWDQVPFSV',
    'ImmTAC':'YLEPGPVTA',
    'canonical_KTW':'KTWGQYWQV',
    'canonical_ITD':'ITDQVPFSV'
}
# peptides = {
#     'non_canonical':'SEQDPENRAW',
#     'canonical':'SEQDPEKAWGA'
# }
# peptides = {
#     'non_canonical':'QAQLEVSVQY',
#     'canonical':'DADFDEVQY'
# }
col = 'percentile'

# all_txts = subprocess.run('find {} -type f -name "msmsScans_new.txt"'.format(rootdir),shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
# tumor_dfs = []
# tumor_raws = []
# for txt in tqdm(all_txts):
#     msms = pd.read_csv(txt,sep='\t')
#     msms = msms.loc[msms['Identified']=='+',:]
#     for raw,sub_df in msms.groupby(by='Raw file'):
#         each_raw_data = []
#         for p,sub_df2 in sub_df.groupby(by='Sequence'):
#             intensity = sub_df2['Precursor intensity'].values.max()
#             each_raw_data.append((p,intensity))
#         each_raw_df = pd.DataFrame.from_records(data=each_raw_data,columns=['peptide','intensity'])
#         each_raw_df = each_raw_df.loc[each_raw_df['intensity'].notna(),:]
#         upper = np.quantile(each_raw_df['intensity'].values,0.75)
#         each_raw_df['norm'] = np.log2(each_raw_df['intensity'].values/upper)
#         # also add percentile
#         each_raw_df = each_raw_df.sort_values(by='intensity')
#         each_raw_df['percentile'] = [(i+1)/each_raw_df.shape[0] for i in range(each_raw_df.shape[0])]
#         # add to list
#         tumor_dfs.append(each_raw_df)
#         tumor_raws.append(raw)
# final = pd.concat(tumor_dfs,axis=0,keys=tumor_raws).reset_index(level=-2).rename(columns={'level_0':'raw_file'})
# final['log_intensity'] = np.log2(final['intensity'].values)
# final.to_csv('final_melanoma.txt',sep='\t',index=None)


final = pd.read_csv('final_melanoma.txt',sep='\t')
data = {}
df_list = []
for name,seq in peptides.items():
    df = final.loc[final['peptide']==seq,:]
    data[name] = df[col].values
    df_list.append(df)
selected_df = pd.concat(df_list,axis=0)
selected_df.to_csv('check.txt',sep='\t')

# compare_data = {}
# for raw,sub_df in selected_df.groupby(by='raw_file'):
#     if len(sub_df['peptide'].unique()) == 4:
#         for name,seq in peptides.items():
#             a = sub_df.loc[sub_df['peptide']==seq,col].values[0]
#             compare_data.setdefault(name,[]).append(a)

# ratio1 = np.median(compare_data['PMEL_splicing']) / np.median(compare_data['ImmTAC'])
# ratio2 = np.median(compare_data['PMEL_splicing']) / np.median(compare_data['canonical_KTW'])
# ratio3 = np.median(compare_data['PMEL_splicing']) / np.median(compare_data['canonical_ITD'])
# print(ratio1,ratio2,ratio3,np.median(compare_data['PMEL_splicing']),np.median(compare_data['ImmTAC']),np.median(compare_data['canonical_KTW']),np.median(compare_data['canonical_ITD']))

# compare_data = {}
# for raw,sub_df in selected_df.groupby(by='raw_file'):
#     for name,seq in peptides.items():
#         try:
#             a = sub_df.loc[sub_df['peptide']==seq,col].values[0]
#         except:
#             a = 0
#         compare_data.setdefault(name,[]).append(a)
# compare_data_new = {}
# for k,v in compare_data.items():
#     compare_data_new[k] = [item for item in v if item > 0]
# compare_data = compare_data_new

# # print(ttest_ind(compare_data['non_canonical'],compare_data['canonical']))
# fig,ax = plt.subplots()
# sns.boxplot(list(compare_data.values()),ax=ax)
# ax.set_xticks((0,1,2,3))
# ax.set_xticklabels(list(compare_data.keys()))
# ax.set_ylabel(col)
# # ax.set_ylim((-2.5,4))
# plt.savefig('abundance_{}_4_full.pdf'.format(col),bbox_inches='tight')
# plt.close()

# plot by rank
highest_list = []
for raw,sub_df in selected_df.groupby(by='raw_file'):
    highest = sub_df.sort_values(by=col,ascending=False)['peptide'].iloc[0]
    highest_list.append(highest)
a,b = np.unique(highest_list,return_counts=True)
dic = {}
for i1,i2 in zip(a,b):
    dic[i1] = i2
dic['ITDQVPFSV'] = 0
dic = dict(sorted(dic.items(), key=lambda x: x[1],reverse=True))
fig,ax = plt.subplots()
ax.bar(x=np.arange(len(dic)),height=dic.values())
ax.set_ylim([0,30])
plt.savefig('abudance_champion.pdf',bbox_inches='tight')
plt.close()
