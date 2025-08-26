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
import anndata as ad
from scipy.sparse import csr_matrix
import re

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

gtf_path = '/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/NeoVerse/data/Homo_sapiens.GRCh38.110.gtf'
gtf = pd.read_csv(gtf_path,sep='\t',skiprows=5,header=None)
gtf.columns = ['chrom','source','feature','start','end','score','strand','phase','attrs']
enst2ensg = {}
enst2coord = {}
enst2strand = {}
enst2chrom = {}
enst2dic = {}  
enst2cds_seq = {}

pat = re.compile(r'transcript_id "(ENST\d+)";')
col = []
for item in gtf['attrs']:
    match = re.search(pat,item)
    if match:
        tid = match.group(1)
    else:
        tid = None
    col.append(tid)
gtf['tid'] = col    

pat = re.compile(r'gene_id "(ENSG\d+)";')
col = []
for item in gtf['attrs']:
    match = re.search(pat,item)
    if match:
        gid = match.group(1)
    else:
        gid = None
    col.append(gid)
gtf['gid'] = col  

# filter by chrom
valid_chrom = set(['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12',
                   'chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM'])

gtf['chrom'] = ['chr' + str(item) for item in gtf['chrom']]
gtf = gtf.loc[gtf['chrom'].isin(valid_chrom),:]


# fill enst2cds_seq
with open('Homo_sapiens.GRCh38.cds.all.fa','r') as in_handle:
    for title,seq in SimpleFastaParser(in_handle):
        enst = title.split('.')[0]
        enst2cds_seq[enst] = seq


# main
for tid,sub_df in tqdm(gtf.groupby(by='tid')):
    if (not tid is None) and ('CDS' in sub_df['feature'].values):  
        sub_df = sub_df.loc[sub_df['feature']=='CDS',:]
        strand = sub_df['strand'].iloc[0]
        chrom = sub_df['chrom'].iloc[0]
        ensg = sub_df['gid'].iloc[0]
        enst2ensg[tid] = ensg
        enst2strand[tid] = strand
        enst2chrom[tid] = chrom
        coord = []  # [(a,b),(c,d)]
        dic = {} # {1:(pos_of_first_base_in_that_cds,pos_of_first_base_in_total_cds,number_of_already_accumulated_residues,last_stop_phase)}
        
        accumulated_cds = 1  # the first nt of this cds, should be this number of pos in total cds 
        accumulated_residue = 0   # after iterating each cds, how many residues have already been completed
        stop_phase = 0   # ensembl web format, 1 means 1 additional nt will go for the junction codon

        for i,row in enumerate(sub_df.itertuples()):   
            coord.append((row.start,row.end))
            if strand == '+':
                dic[i+1] = (int(row.start),accumulated_cds,accumulated_residue,stop_phase)
            elif strand == '-':   # becasue always record the start pos (in the sense of translation) of each cds, also this gtf neg strand not sorted by coord, manual check to make sense
                dic[i+1] = (int(row.end),accumulated_cds,accumulated_residue,stop_phase)
            span = int(row.end) - int(row.start) + 1  # how many nt of this cds 
            accumulated_cds += span
            span = span + stop_phase  # add the stop_phase to factor in given nt
            n_residue = span // 3
            phase = span % 3
            accumulated_residue += n_residue
            stop_phase = phase

        try:
            assert stop_phase == 0
        except:   # seems to be CDS incomplete
            coord = None
            dic = None
        enst2coord[tid] = coord
        enst2dic[tid] = dic



# assemble
df = pd.Series(enst2ensg,name='ensg').to_frame()
df['strand'] = df.index.map(enst2strand).values
df['chrom'] = df.index.map(enst2chrom).values
df['coord'] = df.index.map(enst2coord).values
df['dic'] = df.index.map(enst2dic).values
df['cds_seq'] = df.index.map(enst2cds_seq).values

# annotate canonical
can = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/gene/canonical_enhanced.txt',sep='\t',index_col=0)
all_can = set(can.index.tolist())
df['is_canonical'] = [True if item in all_can else False for item in df.index]

df.to_csv('hg38_cds_final.txt',sep='\t')


