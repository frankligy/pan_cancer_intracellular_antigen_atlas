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

root_atlas_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas'

df = pd.read_csv('./stats/final_all_ts_antigens.txt',sep='\t')
p2t = pd.Series(index=df['pep'].values,data=df['typ'].values).to_dict()
p2s = pd.Series(index=df['pep'].values,data=df['source'].values).to_dict()
vc = df['pep'].value_counts()
data = []
for p in vc.index:
    tmp = df.loc[df['pep']==p,:]
    cs = tmp['cancer'].values.tolist()
    data.append((p,','.join(cs),len(cs),p2t[p],p2s[p]))
result = pd.DataFrame.from_records(data,columns=['peptide','cancers','recurrency','category','source'])
result.to_csv('self_gene_common/all_common_by_cancer.txt',sep='\t',index=None)