#!/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/bin/python3.7

import pandas as pd
import numpy as np
import sys,os
import subprocess
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import re

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

# 30 tissues
# website miss stomach
# compared to the peptide.py which is the folder name aman defined, here needs to change a few to match the readme

tissues = [
    'AdrenalGland',
    'Aorta',
    'Bladder',
    'BoneMarrow',
    'Brain',
    'Cerebellum',
    'Colon',
    'Esophagus',
    'Gallbladder',
    'Heart',
    'Kidney',
    'Liver',
    'Lung',
    'LymphNode',
    'Mamma',
    'Muscle',
    'Myelon',
    'Ovary',
    'Pancreas',
    'Prostate',
    'Skin',
    'SmallIntestine',
    'Spleen',
    'Stomach',
    'Testis',
    'Thymus',
    'Thyroid',
    'Tongue',
    'Trachea',
    'Uterus',
]

df = pd.read_csv('all_raw_files.txt',sep='\t',header=None)
df.columns = ['id','file','ftp','type','other']

anno = {
    'DN02':'A*11:01; A*68:01; B*15:01; B*35:03; C*03:03; C*04:01',
    'DN03':'A*01:01; A*11:01; B*15:01; B*35:01; C*03:03; C*04:01',
    'DN04':'A*02:01; A*23:01; B*27:05; B*50:01; C*02:02; C*06:02',
    'DN05':'A*01:01; A*11:01; B*07:02; B*49:01; C*07:01; C*07:02',
    'DN06':'A*03:01; A*68:02; B*07:02; B*14:02; C*07:02; C*08:02',
    'DN08':'A*32:01; A*68:01; B*15:01; B*44:02; C*03:03; C*07:04',
    'DN09':'A*24:02; A*30:01; B*13:02; B*35:08; C*04:01; C*06:02',
    'DN11':'A*01:01; A*69:01; B*37:01; B*49:01; C*06:02; C*07:01',
    'DN12':'A*02:01; A*11:01; B*15:01; B*35:01; C*03:04; C*04:01',
    'DN13':'A*02:05; A*11:01; B*40:02; B*58:01; C*02:02; C*07:01',
    'DN14':'A*02:01; A*68:02; B*14:02; B*27:05; C*02:02; C*08:02',
    'DN15':'A*01:01; A*02:01; B*08:01; B*44:02; C*07:01; C*07:04',
    'DN16':'A*01:01; A*24:02; B*08:01; B*41:01; C*07:01; C*17:01',
    'DN17':'A*03:01; A*24:02; B*35:03; B*45:01; C*04:01; C*16:01',
    'DN278':'A*02:01; B*44:02; C*05:01',
    'DN281':'A*11:01; A*26:01; B*08:01; B*35:01; C*07:02; C*04:01',
    'DN1':'A*03:01; A*29:02; B*07:02; B*44:03; C*07:02; C*16:01',
    'DN3':'A*24:02; A*25:01; B*18:01; B*41:01; C*12:03; C*17:01',
    'DN4':'A*02:01; A*26:08; B*15:01; B*44:02; C*03:04; C*05:01',
    'DN5':'A*01:01; A*03:01; B*07:06; B*07:02; C*07:02; C*15:05',
    'DN6':'A*01:01; A*25:01; B*13:02; B*39:01; C*06:02; C*12:03'
}

pat = re.compile(r'-(DN\d+)_')

data = []
for t in tissues:
    sub = df.loc[df['file'].str.contains(t),:]
    for f in sub['file']:
        uid = re.search(pat,f).group(1)
        hla = anno[uid]
        data.append((t,f,uid,hla))
final = pd.DataFrame.from_records(data,columns=['tissue','file','uid','hla'])
final.to_csv('final.txt',sep='\t',index=None)



