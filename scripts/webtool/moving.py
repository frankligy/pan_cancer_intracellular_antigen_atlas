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

# you need to move figures, final and meta

root_atlas_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas'
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


# move final
df = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/codes/summary/stats/final_all_ts_antigens.txt',sep='\t')
wl = {}
for c,sub_df in df.groupby(by='cancer'):
    wl[c] = set(sub_df['pep'])
des_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/codes/webtool/app/static'
if not os.path.exists(des_dir):
    os.mkdir(des_dir)
for c in cancers:
    final_path = os.path.join(root_atlas_dir,c,'antigen','fdr','final_enhanced.txt')
    final = pd.read_csv(final_path,sep='\t')
    valid_peps = wl[c]
    final = final.loc[final['pep'].isin(valid_peps),:]
    after = os.path.join(des_dir,'{}_final_enhanced.txt'.format(c))
    final.to_csv(after,sep='\t',index=None)


# move meta
des_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/codes/webtool/app/static'
if not os.path.exists(des_dir):
    os.mkdir(des_dir)
root_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/immunopeptidome'
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
    'SKCM':'melanoma'
}
for c in cancers:
    final_path = os.path.join(root_dir,mapping[c],'metadata.txt')
    previous = final_path
    after = os.path.join(des_dir,'{}_metadata.txt'.format(c))
    subprocess.run('cp {} {}'.format(previous,after),shell=True)


# special
freq = '/gpfs/data/yarmarkovichlab/medulloblastoma/neoverse_folder/NeoVerse_final_output_new/antigens/US_HLA_frequency.csv'
subprocess.run('cp {} {}'.format(freq,des_dir),shell=True)

deepimmuno_result = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/codes/webtool/all_deepimmuno_immunogenicity.txt'
subprocess.run('cp {} {}'.format(deepimmuno_result,des_dir),shell=True)


# move figures
des_dir = '/gpfs/data/yarmarkovichlab/public/ImmunoVerse/assets'
if not os.path.exists(des_dir):
    os.mkdir(des_dir)
for c in cancers:
    print(c)
    dir_path = os.path.join(root_atlas_dir,c,'assets')
    old_dir = os.getcwd()
    os.chdir(dir_path)
    all_pngs = subprocess.run('for f in *.png; do echo $f; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
    os.chdir(old_dir)
    for png in tqdm(all_pngs):
        previous = os.path.join(dir_path,png)
        after = os.path.join(des_dir,'{}_{}'.format(c,png))
        subprocess.run('cp {} {}'.format(previous,after),shell=True)
