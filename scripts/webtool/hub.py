#!/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/bin/python3.7

import pandas as pd
import numpy as np
import sys,os
from tqdm import tqdm
import mygene
import subprocess
import re
import argparse
import pickle
import pysam
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import bisect
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from collections import Counter
import anndata as ad
from scipy.sparse import csr_matrix
from collections import Counter
from ast import literal_eval
import math
import argparse
import json
from copy import deepcopy

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


mapping = {
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
    'SKCM':'melanoma',
    'BALL':'B_ALL',
    'TALL':'T_ALL',
    'CML':'CML',
    'CLL':'CLL',
    'MCL':'MCL',
    'FL':'follicular_lymphoma',
    'meningioma':'meningioma',
    'UCEC':'endometrioid_cancer',
    'ependymoma':'ependymoma',
    'EWS':'erwin_sarcoma',
    'MM':'multiple_myeloma',
    'PRAD':'prostate_cancer',
    'OS':'osteosarcoma',
    'UCS':'uterine_serous_carcinoma',
    'schwannoma':'schwannoma'
}

root_dir = '/gpfs/data/yarmarkovichlab/public/ImmunoVerse/ImmunoVerse_Hub'

df_list = []
for c in mapping.keys():
    print(c)
    meta = pd.read_csv(os.path.join(root_dir,'{}_metadata.txt'.format(c)),sep='\t')
    sbatch = os.path.join(root_dir,'{}_download.sbatch'.format(c))
    col = []
    for row in meta.itertuples():
        f = row.sample
        s = row.study
        if f.endswith('RAW'):
            print(f)
        cmd = 'grep "{}" "{}"'.format(f,sbatch)
        r = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]   
        if len(r) >= 1:
            r = [item.lstrip('#').lstrip(' ').rstrip(' ') for item in r if '://' in item]
            l = ';'.join(r)
            col.append(l)
        else:
            # raw RAW
            cmd = 'grep "{}" "{}"'.format(f.replace(r'.raw',r'.RAW'),sbatch)
            r = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1] 
            if len(r) >= 1:
                r = [item.lstrip('#').lstrip(' ').rstrip(' ') for item in r if '://' in item]
                l = ';'.join(r)
                col.append(l)
            else:
                # % #
                cmd = 'grep "{}" "{}"'.format(f.replace(r'#',r'').replace(r'%',r'_'),sbatch)
                r = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1] 
                if len(r) >= 1:
                    r = [item.lstrip('#').lstrip(' ').rstrip(' ') for item in r if '://' in item]
                    l = ';'.join(r)
                    col.append(l)
                else:
                    col.append(None)
            
    meta['download_link'] = col
    df_list.append(meta)
all_meta = pd.concat(df_list,axis=0,keys=list(mapping.keys())).reset_index(level=-2).rename(columns={'level_0':'cancer'})
all_meta.to_csv(os.path.join(root_dir,'sdrf.txt'),sep='\t',index=None)

with open(os.path.join(root_dir,'sdrf_notes.txt'),'w') as f:
    f.write('1. For files without download link, usually it is due to whole study was uploaded as a compressed zip file, or single file, when compressed, has a different name, please resort to each download.sbatch file to get the link and fix these discrepencies\n')
    f.write('2. PRIDE seems to work better using https instead of ftp, in case that happens, just replace ftp with https or vice versa\n')
    f.write('3. Massive is constantly changing the root server URL, when error occurs, search MSV ID and gran the latest server URL\n')
    f.write('4. Certain files were stored as RAW on server, I uniformally rename them for lower case raw file\n')
    f.write('5. URL can not have special character such as "#", "%", etc, which can explain some discrepencies between link and file name\n')

    






