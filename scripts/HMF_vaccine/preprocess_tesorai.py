#!/gpfs/data/yarmarkovichlab/Frank/pan_cancer/antigen_portal/spectrum_env/bin/python3.8

import sys
import os
import numpy as np
import pandas as pd
from pyteomics import mzml
from pyteomics import mgf
import matplotlib.pyplot as plt
import matplotlib as mpl
import multiprocessing as mp
import subprocess
from tqdm import tqdm
import json
import argparse 
from pyteomics.mass import calculate_mass
import re
from ast import literal_eval
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from io import StringIO
import bisect


mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

# preprocess

# df = pd.read_csv('../pan_cancer_melanoma_vaccine_outputs_tesorai_melanoma_new_peptide_fdr_005.tsv',sep='\t')
# df['filename'] = [item.split('.raw')[0] for item in df['filename']]
# df.to_csv('../immuno/pan_cancer/combined/txt/other_alg.txt',sep='\t',index=None)

# generate hla_dic, it should be "pan_cancer@A375_OE3":["A*02:01","B*18:01","B*15:01","C*07:01","C*03:04"]
meta = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/immunopeptidome/melanoma/metadata.txt',sep='\t')
for item1,item2 in zip(meta['sample'],meta['HLA']):
    part1 = '"' + "pan_cancer" + "@" + item1.split('.raw')[0] + '"' + ':' 
    try:
        part2 = '[' + ','.join(['"' + item + '"' for item in item2.split(',')]) + ']'
    except:
        part2 = '[]'
    assemble = part1 + part2 + ','
    print(assemble)






