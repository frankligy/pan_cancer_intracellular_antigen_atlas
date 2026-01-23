#!/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/bin/python3.7

import pandas as pd
import numpy as np
import sys,os
import subprocess
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import anndata as ad
from Bio.SeqIO.FastaIO import SimpleFastaParser
from tqdm import tqdm
import itertools

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

df = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas/nbl_neg/antigen/0.01/final_enhanced.txt',sep='\t')
cond = [False if ('[]' in item) and ('(\'HLA-' not in item) else True for item in df['presented_by_each_sample_hla']]
df = df.loc[cond,:]
df = df.loc[df['typ']=='pathogen',:]
df = df.loc[df['source'].str.contains('NIACI'),:]

fasta_dic = {}
for mer in [8,9,10,11,12,13,14,15]:
    lis = []
    with open('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/codes/summary/ensembl_protein_{}.fasta'.format(mer),'r') as in_handle:
        for title,seq in SimpleFastaParser(in_handle):
            lis.append(seq)
    lis = set(lis)
    fasta_dic[mer] = lis

col = []
for pep,typ in tqdm(zip(df['pep'],df['typ']),total=df.shape[0]):
    if typ == 'self_gene':
        col.append(False)
    else:
        if 'I' in pep or 'L' in pep:
            options = []
            for c in pep:
                if c == 'I' or c == 'L':
                    options.append(['I','L'])
                else:
                    options.append([c])
            expanded = [''.join(p) for p in itertools.product(*options)]
            l = len(pep)
            lis = fasta_dic[l]
            flag = False
            for item in expanded:
                if item in lis:
                    flag = True
                    break
            col.append(flag)
                
        else:
            col.append(False)

df['is_ambiguous_IL'] = col
df = df.loc[~df['is_ambiguous_IL'],:]
df.to_csv('nbl_neg_niaci.txt',sep='\t',index=None)


