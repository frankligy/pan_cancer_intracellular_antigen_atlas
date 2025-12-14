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

anno = pd.read_csv('ely_et_al_anno.txt',sep='\t',index_col=0)
delete = ['SRR15513220','SRR15513224','SRR15513225','SRR15513226']
second = ['SRR1528686','SRR1528687','SRR1528688','SRR1528689']

for srr in anno.index:
    if srr in delete:
        continue
    elif srr in second:
        with open('./result/{}/outputDir/offset.correction.parameters.txt'.format(srr),'w') as f:
            f.write('26\t12\n27\t12\n28\t12\n')
    else:
        with open('./result/{}/outputDir/offset.correction.parameters.txt'.format(srr),'w') as f:
            f.write('25\t12\n26\t12\n27\t12\n')









