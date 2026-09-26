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


mapping = {
    'NYU2507_rerun':'8-21-2024_20240821_TT_Evo_IC_2507_rerun_2923',
    'NYU2479_rerun':'8-21-2024_20240821_TT_Evo_IC_2479_rerun_2924'
}

df = pd.read_csv('../lung_cancer_quantified_psm_fdr.tsv',sep='\t')
os.makedirs('../immuno/dummy/combined/txt')
df['filename'] = [mapping[item.split('.')[0]] for item in df['filename']]
df = df.loc[~df['is_decoy'],:]
df.to_csv('../immuno/dummy/combined/txt/other_alg.txt',sep='\t',index=None)
