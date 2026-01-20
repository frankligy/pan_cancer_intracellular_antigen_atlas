#!/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/bin/python3.7

import pandas as pd
import numpy as np
import sys,os
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import SeqIO
import subprocess
import matplotlib.pyplot as plt
import matplotlib as mpl
from tqdm import tqdm
import pycircos
from Bio.Seq import Seq

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


df = pd.read_excel('scitranslmed.ads6335_data_files_s1_to_s15.v2/ads6335_data_file_s12.xlsx').set_index(keys='HopkinsID')
hid2tumor = df['project_id'].to_dict()
hid2st = df['sample_type'].to_dict()

count = pd.read_csv('scitranslmed.ads6335_data_files_s1_to_s15.v2/ads6335_data_file_s4.csv',index_col=0)
with open('scitranslmed.ads6335_data_files_s1_to_s15.v2/data_s3_columns.txt','w') as f:
    for s in count.columns:
        f.write('{}\n'.format(s))
count = count.loc[:,['Niallia circulans','Campylobacter ureolyticus','Clostridium intestinale','Helicobacter pylori','Fusobacterium nucleatum','Alphapapillomavirus 9','Alphapapillomavirus 7']]
count['tumor'] = [hid2tumor[item] for item in count.index]
count['sample_type'] = [hid2st[item] for item in count.index]
count.to_csv('check.txt',sep='\t')

