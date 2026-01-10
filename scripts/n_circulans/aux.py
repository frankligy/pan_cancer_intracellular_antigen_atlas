#!/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/bin/python3.7

import pandas as pd
import numpy as np
import sys,os
import subprocess
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import math

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

df = pd.read_csv('SraRunTable.csv')

# # AIR
# sub_df = df.loc[df['Specimen_Name'].str.upper().str.contains('AIR'),:]
# for srr in sub_df['Run']:
#     print(srr)
#     ori_r1 = '{}_1.fastq.gz'.format(srr)
#     ori_r2 = '{}_2.fastq.gz'.format(srr)
#     des_r1 = '{}_1.fastq'.format(srr)
#     des_r2 = '{}_2.fastq'.format(srr)
#     cmd = 'gunzip -c ./raw/{} > ./AIR/{}'.format(ori_r1,des_r1)
#     subprocess.run(cmd,shell=True)
#     cmd = 'gunzip -c ./raw/{} > ./AIR/{}'.format(ori_r2,des_r2)
#     subprocess.run(cmd,shell=True)

# # PORT
# sub_df = df.loc[df['Specimen_Name'].str.upper().str.contains('PORT'),:]
# for srr in sub_df['Run']:
#     print(srr)
#     ori_r1 = '{}_1.fastq.gz'.format(srr)
#     ori_r2 = '{}_2.fastq.gz'.format(srr)
#     des_r1 = '{}_1.fastq'.format(srr)
#     des_r2 = '{}_2.fastq'.format(srr)
#     cmd = 'gunzip -c ./raw/{} > ./PORT/{}'.format(ori_r1,des_r1)
#     subprocess.run(cmd,shell=True)
#     cmd = 'gunzip -c ./raw/{} > ./PORT/{}'.format(ori_r2,des_r2)
#     subprocess.run(cmd,shell=True)

# # FTO
# cond = []
# for item in df['Specimen_Name']:
#     item = item.upper()
#     if 'RTO' in item or 'LTO' in item or 'ROV' in item or 'LOV' in item or 'RFT' in item or 'LFT' in item:
#         cond.append(True)
#     else:
#         cond.append(False)
# sub_df = df.loc[cond,:]
# for srr in sub_df['Run']:
#     print(srr)
#     ori_r1 = '{}_1.fastq.gz'.format(srr)
#     ori_r2 = '{}_2.fastq.gz'.format(srr)
#     des_r1 = '{}_1.fastq'.format(srr)
#     des_r2 = '{}_2.fastq'.format(srr)
#     cmd = 'gunzip -c ./raw/{} > ./FTO/{}'.format(ori_r1,des_r1)
#     subprocess.run(cmd,shell=True)
#     cmd = 'gunzip -c ./raw/{} > ./FTO/{}'.format(ori_r2,des_r2)
#     subprocess.run(cmd,shell=True)

# # PG
# cond = []
# for item in df['Specimen_Name']:
#     item = item.upper()
#     if '-PG' in item or '-RPG' in item or '-LPG' in item:
#         cond.append(True)
#     else:
#         cond.append(False)
# sub_df = df.loc[cond,:]
# for srr in sub_df['Run']:
#     print(srr)
#     ori_r1 = '{}_1.fastq.gz'.format(srr)
#     ori_r2 = '{}_2.fastq.gz'.format(srr)
#     des_r1 = '{}_1.fastq'.format(srr)
#     des_r2 = '{}_2.fastq'.format(srr)
#     cmd = 'gunzip -c ./raw/{} > ./PG/{}'.format(ori_r1,des_r1)
#     subprocess.run(cmd,shell=True)
#     cmd = 'gunzip -c ./raw/{} > ./PG/{}'.format(ori_r2,des_r2)
#     subprocess.run(cmd,shell=True)

# # PG
# cond = []
# for item in df['Specimen_Name']:
#     item = item.upper()
#     if '-CX' in item:
#         cond.append(True)
#     else:
#         cond.append(False)
# sub_df = df.loc[cond,:]
# for srr in sub_df['Run']:
#     print(srr)
#     ori_r1 = '{}_1.fastq.gz'.format(srr)
#     ori_r2 = '{}_2.fastq.gz'.format(srr)
#     des_r1 = '{}_1.fastq'.format(srr)
#     des_r2 = '{}_2.fastq'.format(srr)
#     cmd = 'gunzip -c ./raw/{} > ./CX/{}'.format(ori_r1,des_r1)
#     subprocess.run(cmd,shell=True)
#     cmd = 'gunzip -c ./raw/{} > ./CX/{}'.format(ori_r2,des_r2)
#     subprocess.run(cmd,shell=True)

# # NTC
# cond = []
# for item in df['Specimen_Name']:
#     item = item.upper()
#     if 'NTC' in item:
#         cond.append(True)
#     else:
#         cond.append(False)
# sub_df = df.loc[cond,:]
# for srr in sub_df['Run']:
#     print(srr)
#     ori_r1 = '{}_1.fastq.gz'.format(srr)
#     ori_r2 = '{}_2.fastq.gz'.format(srr)
#     des_r1 = '{}_1.fastq'.format(srr)
#     des_r2 = '{}_2.fastq'.format(srr)
#     cmd = 'gunzip -c ./raw/{} > ./NTC/{}'.format(ori_r1,des_r1)
#     subprocess.run(cmd,shell=True)
#     cmd = 'gunzip -c ./raw/{} > ./NTC/{}'.format(ori_r2,des_r2)
#     subprocess.run(cmd,shell=True)

# # DB
# cond = []
# for item in df['Specimen_Name']:
#     item = item.upper()
#     if 'DB' in item:
#         cond.append(True)
#     else:
#         cond.append(False)
# sub_df = df.loc[cond,:]
# for srr in sub_df['Run']:
#     print(srr)
#     ori_r1 = '{}_1.fastq.gz'.format(srr)
#     ori_r2 = '{}_2.fastq.gz'.format(srr)
#     des_r1 = '{}_1.fastq'.format(srr)
#     des_r2 = '{}_2.fastq'.format(srr)
#     cmd = 'gunzip -c ./raw/{} > ./DB/{}'.format(ori_r1,des_r1)
#     subprocess.run(cmd,shell=True)
#     cmd = 'gunzip -c ./raw/{} > ./DB/{}'.format(ori_r2,des_r2)
#     subprocess.run(cmd,shell=True)



'''analysis'''
data = []
stream_of_data = []
conditions = ['FTO','CX','PG','PORT','AIR','NTC','DB']
ns = [369,152,122,81,130,111,36]
taxa_level = 'Genus'
name = 'Clostridium'
for c in conditions:
    taxa = pd.read_csv('{}/df_taxa.txt'.format(c),sep='\t',index_col=0)
    taxa = taxa.loc[taxa[taxa_level]==name,:]
    all_asv = taxa.index.tolist()
    count = pd.read_csv('{}/df.txt'.format(c),sep='\t',index_col=0)
    count = count.loc[:,all_asv]
    s = count.sum(axis=1)
    stream_of_data.append(s.values.tolist())
    for item in s:
        data.append((c,item))
df = pd.DataFrame.from_records(data=data,columns=['condition','count'])
fig,ax = plt.subplots()
sns.stripplot(data=df,x='condition',y='count',ax=ax)
ax.boxplot(x=stream_of_data,positions=np.arange(len(conditions)),patch_artist=False,showfliers=False)
l = ['{}\nn={}'.format(i1,i2) for i1,i2 in zip(conditions,ns)]
ax.set_xticklabels(l)
plt.savefig('{}_{}.pdf'.format(taxa_level,name),bbox_inches='tight')
plt.close()
    

