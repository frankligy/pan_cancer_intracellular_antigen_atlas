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

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

anno = pd.read_csv('offset_table.txt',sep='\t',index_col=0)
anno = anno.loc[anno['offset'].notna(),:]

for row in anno.itertuples():
    srr = row.Index
    with open('./result/{}/outputDir/offset.correction.parameters.txt'.format(srr),'w') as f:
        for item in row.offset.split(';'):
            l,o = item.split(',')
            f.write('{}\t{}\n'.format(l,o))

        










