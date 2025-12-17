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

cancers2immuno = {
    'BRCA':'breast_cancer',
    'KIRC':'kidney_clear_cell',
    'COAD':'colon_cancer',
    'STAD':'stomach_cancer',
    'MESO':'mesothelioma',
    'LIHC':'liver_cancer',
    'ESCA':'esophageal_cancer',
    'CESC':'cervical_cancer',
    'BLCA':'bladder_cancer',
    'RT':'rhabdoid_tumor',
    'AML':'AML',
    'DLBC':'DLBC',
    'GBM':'GBM',
    'NBL':'neuroblastoma',
    'PAAD':'pancreatic_cancer',
    'HNSC':'head_and_neck',
    'OV':'ovarian_cancer',
    'LUSC':'lung_LUSC',
    'LUAD':'lung_LUAD',
    'CHOL':'bile_duct_cancer',
    'SKCM':'melanoma'
}

# # make sure no suffix added
# all_tesorai = subprocess.run("for file in tesorai_peptide_fdr_*.tsv; do echo $file; done",shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
# for tesorai in all_tesorai:
#     df = pd.read_csv(tesorai,sep='\t')
#     all_file = set([item.split('.')[0] for item in df['filename']])
#     c = tesorai.split('_')[3]
#     if c in cancers2immuno.keys():
#         meta = pd.read_csv(os.path.join('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/immunopeptidome',cancers2immuno[c],'metadata.txt'),sep='\t')
#         doced_file = set([item.split('.')[0] for item in meta['sample']])
#         try:
#             assert len(all_file.difference(doced_file)) == 0
#         except AssertionError:
#             print(tesorai,print(all_file.difference(doced_file)))

# # fix suffix
# df = pd.read_csv('tesorai_peptide_fdr_PAAD_rescue_1.tsv',sep='\t')
# col = []
# for item in df['filename']:
#     b = item.split('.raw')[0]
#     bb = '_'.join(b.split('_')[:-1])
#     assemble = '{}.raw'.format(bb)
#     col.append(assemble)
# df['filename'] = col
# df.to_csv('tesorai_peptide_fdr_PAAD_rescue_1.tsv',sep='\t',index=None)



# make sure no RAW no ZIP
all_tesorai = subprocess.run("for file in tesorai_peptide_fdr_*.tsv; do echo $file; done",shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
for tesorai in all_tesorai:
    df = pd.read_csv(tesorai,sep='\t')
    cond = np.any(df['filename'].str.endswith('RAW').values)
    print(tesorai,cond)


# tesorai_peptide_fdr_CANCER_suffix.tsv
# tech not mixed

mapping = {
    'frank_brca_rescue_1_pep_fdr.tsv':'tesorai_peptide_fdr_BRCA_rescue_1.tsv',
    'frank_kirc_rescue_1_pep_fdr.tsv':'tesorai_peptide_fdr_KIRC_rescue_1.tsv',
    'frank_kirc_rescue_2_pep_fdr.tsv':'tesorai_peptide_fdr_KIRC_rescue_2.tsv',
    'frank_coad_rescue_1_pep_fdr.tsv':'tesorai_peptide_fdr_COAD_rescue_1.tsv',
    'frank_paad_rescue_1_pep_fdr.tsv':'tesorai_peptide_fdr_PAAD_rescue_1.tsv',
    'frank_hnsc_rescue_1_pep_fdr.tsv':'tesorai_peptide_fdr_HNSC_rescue_1.tsv',
    'frank_hnsc_rescue_2_pep_fdr.tsv':'tesorai_peptide_fdr_HNSC_rescue_2.tsv'
}

for k,v in mapping.items():
    df = pd.read_csv('./rescue_raw/{}'.format(k),sep='\t',index_col=0)
    df = df.loc[~df['is_decoy'],:]
    df = df.loc[df['possible_protein_ids'].notna(),:]
    df['filename'] = [item.split('.zip')[0] + '.d' if item.endswith('.zip') else item for item in df['filename']]
    df['filename'] = [item.split('.RAW')[0] + '.raw' if item.endswith('.RAW') else item for item in df['filename']]
    df.drop(columns=['job_id','is_decoy','retention_time_normalized','precursor_charge'],inplace=True)
    df['possible_protein_ids'] = [';;'.join(item.split(';')) for item in df['possible_protein_ids']]
    df.rename(columns={'precursor_intensity':'intensity'},inplace=True)
    df.to_csv('./{}'.format(v),sep='\t',index=None)


