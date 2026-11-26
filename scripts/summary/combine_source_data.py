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
from scipy.stats import ttest_ind,rankdata,mannwhitneyu


mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


fs = [
    'source_data_fig1.xlsx',
    'source_data_fig2.xlsx',
    'source_data_fig3.xlsx',
    'source_data_fig4.xlsx',
    'source_data_fig5.xlsx',
    'source_data_fig6.xlsx',
    'source_data_fig7.xlsx',
    'source_data_figE1.xlsx',
    'source_data_figE2.xlsx',
    'source_data_figE3.xlsx',
    'source_data_figE4.xlsx',
    'source_data_figE5.xlsx',
    'source_data_figE6.xlsx',
    'source_data_figE7.xlsx',
    'source_data_figE8.xlsx',
    'source_data_figE9.xlsx',
    'source_data_figE10.xlsx'
]


with pd.ExcelWriter('./table_source_data/combined_source_data.xlsx') as writer:
    for f in fs:
        fig = f.split('.')[0].split('_')[-1]
        df_dict = pd.read_excel('./table_source_data/{}'.format(f),sheet_name=None)
        for k,v in df_dict.items():
            new_k = '{}_{}'.format(fig,k)
            v.to_excel(writer,sheet_name=new_k,index=None)
