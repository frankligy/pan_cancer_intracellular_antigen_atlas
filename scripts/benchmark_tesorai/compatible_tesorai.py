#!/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/bin/python3.7

import os,sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import re
import subprocess
from tqdm import tqdm
from io import StringIO
import multiprocessing as mp
import math
import bisect
from matplotlib_venn import venn2,venn3

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

# tesorai_peptide_fdr_CANCER_suffix.tsv
# tech not mixed


mapping = {
    'frank_brca_rescue_1_pep_fdr.tsv':'tesorai_peptide_fdr_BRCA_rescue_1.tsv',
    'frank_kirc_rescue_1_pep_fdr.tsv':'tesorai_peptide_fdr_KIRC_rescue_1.tsv',
    'frank_kirc_rescue_2_pep_fdr.tsv':'tesorai_peptide_fdr_KIRC_rescue_2.tsv'
}

for k,v in mapping.items():
    df = pd.read_csv('./rescue_raw/{}'.format(k),sep='\t',index_col=0)
    df = df.loc[~df['is_decoy'],:]
    df['filename'] = [item.split('.zip')[0] + '.d' if item.endswith('.zip') else item for item in df['filename']]
    df.drop(columns=['job_id','is_decoy','retention_time_normalized','precursor_charge'],inplace=True)
    df['possible_protein_ids'] = [';;'.join(item.split(';')) for item in df['possible_protein_ids']]
    df.rename(columns={'precursor_intensity':'intensity'},inplace=True)
    df.to_csv('./{}'.format(v),sep='\t',index=None)