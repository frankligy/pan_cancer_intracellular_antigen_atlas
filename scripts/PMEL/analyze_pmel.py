#!/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/bin/python3.7

import pandas as pd
import numpy as np
import sys,os
import subprocess
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import anndata as ad

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

# calculate the lfc between tumor and skin
dic = {
    'spliced':'chr12:55957671-55958473',
    'can1':'chr12:55957671-55957923',
    'can2':'chr12:55958084-55958473'
}

mapping = {
    'spliced':'ENSG00000185664:E9.5-E11.1',
    'can1':'ENSG00000185664:E9.5-E10.1',
    'can2':'ENSG00000185664:E10.2-E11.1'
}



# skcm
fc_dict1 = {}
tumor = pd.read_csv('../../atlas/SKCM/splicing_all.txt',sep='\t')
data = []
for k,v in dic.items():
    t = tumor.loc[tumor['Unnamed: 0']==v,'count'].values
    for item in t:
        data.append((k,'melanoma',item))
    df = pd.read_csv('{}_normal_data.txt'.format(mapping[k]),sep='\t',index_col=0)
    n1 = df.loc[df['tissue']=='Skin - Not Sun Exposed (Suprapubic)','read_count'].values
    n2 = df.loc[df['tissue']=='Skin - Sun Exposed (Lower leg)','read_count'].values
    for item in n1:
        data.append((k,'skin_not_sun_exposed',item))
    for item in n2:
        data.append((k,'skin_sun_exposed',item))
    fc1 = t.mean()/n1.mean()
    fc2 = t.mean()/n2.mean()
    fc_dict1[k] = fc2
print(fc_dict1)
final = pd.DataFrame.from_records(data,columns=['junction','region','value'])
fig,ax = plt.subplots()
sns.boxplot(final,x='junction',y='value',hue='region',ax=ax)
ax.set_yscale('log')
plt.savefig('pmel_log_skcm.pdf',bbox_inches='tight')
plt.close()

# uvm
fc_dict2 = {}
tumor = pd.read_csv('../../atlas/UVM/splicing_all.txt',sep='\t')
data = []
for k,v in dic.items():
    t = tumor.loc[tumor['Unnamed: 0']==v,'count'].values
    for item in t:
        data.append((k,'melanoma',item))
    df = pd.read_csv('{}_normal_data.txt'.format(mapping[k]),sep='\t',index_col=0)
    n1 = df.loc[df['tissue']=='Skin - Not Sun Exposed (Suprapubic)','read_count'].values
    n2 = df.loc[df['tissue']=='Skin - Sun Exposed (Lower leg)','read_count'].values
    for item in n1:
        data.append((k,'skin_not_sun_exposed',item))
    for item in n2:
        data.append((k,'skin_sun_exposed',item))
    fc1 = t.mean()/n1.mean()
    fc2 = t.mean()/n2.mean()
    fc_dict2[k] = fc2
print(fc_dict2)
final = pd.DataFrame.from_records(data,columns=['junction','region','value'])
fig,ax = plt.subplots()
sns.boxplot(final,x='junction',y='value',hue='region',ax=ax)
ax.set_yscale('log')
plt.savefig('pmel_log_uvm.pdf',bbox_inches='tight')
plt.close()

# uvm vs skcm
df_data = []
for k,v in fc_dict1.items():
    df_data.append((k,v,'SKCM'))
for k,v in fc_dict2.items():
    df_data.append((k,v,'UVM'))
df = pd.DataFrame.from_records(data=df_data,columns=['junction','value','cancer'])
fig,ax = plt.subplots()
sns.pointplot(df,x='cancer',y='value',hue='junction',ax=ax)
plt.savefig('pmel_compare.pdf',bbox_inches='tight')
plt.close()






