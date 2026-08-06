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

df = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/proteome/pdc-client_v1.0.8_Ubuntu_x64/PDC_file_manifest_11302025_103803.tsv',sep='\t')
df.iloc[:3].to_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/proteome/pdc-client_v1.0.8_Ubuntu_x64/manifest_test1.tsv',sep='\t',index=None)