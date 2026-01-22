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
from scipy.stats import spearmanr, pearsonr


mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

# nature
df_patient = pd.read_excel('41586_2023_6706_MOESM3_ESM.xlsx',sheet_name='Patient')
df_pdx = pd.read_excel('41586_2023_6706_MOESM3_ESM.xlsx',sheet_name='PDX')

# NBL
df = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas/NBL/antigen/fdr/msmsScans_all_add_tesorai.txt',sep='\t')
df = df.loc[df['final_identity']=='self_gene',:]
df = df.loc[df['Identified']=='+',:]

# sample
tested_pair = {
    'SK-N-AS':'SK-N-AS',
    'NB-1691':'NB-1691',
    'NB-1771':'NB-1771',
    # 'NB-Ebc1':'NB-Ebc1',
    'NB-SD':'NB_PDX_NB-SD',
    'COG-N-415x':'NB_PDX_COG-N-415x',
    'COG-N-440x':'NB_PDX_COG-N-440x',
    'COG-N-471x':'NB_PDX_COG-N-471x'
}

# tested_pair = {
#     'PALVKK':'NB1_PALVKK',
#     'PANXJL':'NB2_PANXJL',
#     'PAPBJE':'NB3_PAPBJE',
#     'PAPCTS':'NB4_PAPCTS',
#     'PAPVRN':'NB5_PAPVRN',
#     'PARKGJ':'NB6_PARKGJ',
#     'PARZCJ':'NB7_PARZCJ',
#     'PASEGA':'NB8_PASEGA'
# }

# test
for sample1,sample2 in tested_pair.items():
    tmp = df_pdx.loc[df_pdx['Line']==sample1,:]
    dic1 = pd.Series(index=tmp['Peptide'].values,data=tmp['Average area'].values).to_dict()
    tmp = df.loc[df['Raw file'].str.startswith(sample2),:]
    tmp['Precursor intensity'] = tmp['Precursor intensity'].fillna(value=0).values
    dic2 = {}
    for pep,sub_df in tmp.groupby(by='Sequence'):
        a = sub_df['Precursor intensity'].mean()
        dic2[pep] = a

    common = set(dic1.keys()).intersection(set(dic2.keys()))
    x = []
    y = []
    for pep in common:
        added_x = dic1[pep]
        added_y = dic2[pep]
        if added_x > 1 and added_y > 1:
            x.append(added_x)
            y.append(added_y)
    fig,ax = plt.subplots()
    ax.scatter(x,y,s=3)
    ax.set_xscale('log')
    ax.set_yscale('log')
    s1,p1 = spearmanr(x,y)
    s2,p2 = pearsonr(x,y)
    ax.set_title('{}_{}_{}_{}'.format(sample1,s1,s2,len(x)))
    ax.set_xlabel('log10(Thermo_AUC)')
    ax.set_ylabel('log10(MaxQuant_max_precursor_intensity)')
    plt.savefig('{}.pdf'.format(sample1),bbox_inches='tight')
    plt.close()




    


