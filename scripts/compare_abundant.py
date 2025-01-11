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
col = 'norm'

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
#         tumor_dfs.append(each_raw_df)
#         tumor_raws.append(raw)
# final = pd.concat(tumor_dfs,axis=0,keys=tumor_raws).reset_index(level=-2).rename(columns={'level_0':'raw_file'})
# final['log_intensity'] = np.log2(final['intensity'].values)
# final.to_csv('final_breast_cancer.txt',sep='\t',index=None)


final = pd.read_csv('final_melanoma.txt',sep='\t')
data = {}
df_list = []
for name,seq in peptides.items():
    df = final.loc[final['peptide']==seq,:]
    data[name] = df[col].values
    df_list.append(df)
selected_df = pd.concat(df_list,axis=0)

compare_data = {}
for raw,sub_df in selected_df.groupby(by='raw_file'):
    if len(sub_df['peptide'].unique()) == 3:
        for name,seq in peptides.items():
            a = sub_df.loc[sub_df['peptide']==seq,col].values[0]
            compare_data.setdefault(name,[]).append(a)



# print(ttest_ind(compare_data['non_canonical'],compare_data['canonical']))
fig,ax = plt.subplots()
sns.boxplot(list(compare_data.values()),ax=ax)
ax.set_xticks((0,1,2))
ax.set_xticklabels(list(compare_data.keys()))
ax.set_ylabel(col)
ax.set_ylim((-2.5,4))
plt.savefig('abundance_{}.pdf'.format(col),bbox_inches='tight')
plt.close()
