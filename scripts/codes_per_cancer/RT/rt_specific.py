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
import argparse
import json

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'



mutation = pd.read_csv('../../variants/RT/neoantigen_rt_curated.txt',sep='\t')
with open('../../atlas/RT/db_fasta/mutation.fasta','w') as f:
    for row in mutation.itertuples():
        f.write('>{}|{}|69|0.5|ENSG0000|chr{}:1-999|A/G|missense_variant\n{}\n'.format(row.gene,row.position,row.chromosome,row.best_peptide))