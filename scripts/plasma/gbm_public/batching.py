#!/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/bin/python3.7

import pandas as pd
import numpy as np
import sys,os
import subprocess
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import anndata as ad

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

meta_path = '/gpfs/data/yarmarkovichlab/plasma_gbm_PXD008127/immuno/metadata.txt'
immuno_dir_path = '/gpfs/data/yarmarkovichlab/plasma_gbm_PXD008127/immuno'
study_to_batch = 'PXD008127'

meta = pd.read_csv(meta_path,sep='\t')
meta = meta.loc[meta['study']==study_to_batch,:]
total_b = meta['batch'].unique()

# create folder
for b in total_b:
    sub_folder_name = '{}_{}'.format(study_to_batch,b)
    os.mkdir(os.path.join(immuno_dir_path,sub_folder_name))

# move
for b,raw in zip(meta['batch'],meta['sample']):
    sub_folder_name = '{}_{}'.format(study_to_batch,b)
    old_name = '{}/{}/{}'.format(immuno_dir_path,study_to_batch,raw)
    new_name = '{}/{}_{}/{}'.format(immuno_dir_path,study_to_batch,b,raw)
    os.rename(old_name,new_name)


