#!/gpfs/data/yarmarkovichlab/Frank/BayesTS/logit_gate_env/bin/python3.7

import pandas as pd
import numpy as np
import sys,os
import subprocess
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from ast import literal_eval
from tqdm import tqdm
import math
import anndata as ad
import scanpy as sc
from sctriangulate import *
from sctriangulate.preprocessing import *
from sctriangulate.colors import *

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


# # POLR2K_LYLETRSEF_A2402
# f = '05.22.2024 pMHC LYLETRSEF A2402 Biotin (638519801220274399).csv'
# dir_path = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/codes/summary/refold/pMHC_refolding_Rawdata'
# df = pd.read_csv(os.path.join(dir_path,f),sep='\t',encoding='utf-16',skiprows=2)
# fig,ax = plt.subplots()
# ax.plot(df['ml'].values,df['mAU'].values,c='k',lw=1,ls='-')
# ax.set_xlabel('Retention Volumn (mL)')
# ax.set_ylabel('A280 Intensity (mAU)')
# ax.set_title('POLR2K_LYLETRSEF_A2402')
# plt.savefig('POLR2K_LYLETRSEF_A2402.pdf',bbox_inches='tight')
# plt.close()

# # ZNF749_RYLPSSVFL_A2402
# f = '05.17.2024 pMHC RYLPSS-- HLA2402 First Run(638515361501072969).csv'
# dir_path = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/codes/summary/refold/pMHC_refolding_Rawdata'
# df = pd.read_csv(os.path.join(dir_path,f),sep='\t',encoding='utf-16',skiprows=2)
# fig,ax = plt.subplots()
# ax.plot(df['ml'].values,df['mAU'].values,c='k',lw=1,ls='-')
# ax.set_xlabel('Retention Volumn (mL)')
# ax.set_ylabel('A280 Intensity (mAU)')
# ax.set_title('ZNF749_RYLPSSVFL_A2402')
# plt.savefig('ZNF749_RYLPSSVFL_A2402.pdf',bbox_inches='tight')
# plt.close()

# # C19orf48_AYPASLQTL_A2402
# f = '05.10.2024 pMHC AYPASLQTL HLA 2402 wBiotin First Run(638509392801906071).csv'
# dir_path = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/codes/summary/refold/pMHC_refolding_Rawdata'
# df = pd.read_csv(os.path.join(dir_path,f),sep='\t',encoding='utf-16',skiprows=2)
# fig,ax = plt.subplots()
# ax.plot(df['ml'].values,df['mAU'].values,c='k',lw=1,ls='-')
# ax.set_xlabel('Retention Volumn (mL)')
# ax.set_ylabel('A280 Intensity (mAU)')
# ax.set_title('C19orf48_AYPASLQTL_A2402')
# plt.savefig('C19orf48_AYPASLQTL_A2402.pdf',bbox_inches='tight')
# plt.close()

# # TMEM203_STIRVLSGY_A2601
# f = '01.17.2025 pMHC A2601 TMEM203_3UTR(STIR..(638727155330116255).csv'
# dir_path = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/codes/summary/refold/pMHC_refolding_Rawdata'
# df = pd.read_csv(os.path.join(dir_path,f),sep='\t',encoding='utf-16',skiprows=2)
# fig,ax = plt.subplots()
# ax.plot(df['ml'].values,df['mAU'].values,c='k',lw=1,ls='-')
# ax.set_xlabel('Retention Volumn (mL)')
# ax.set_ylabel('A280 Intensity (mAU)')
# ax.set_title('TMEM203_STIRVLSGY_A2601')
# plt.savefig('TMEM203_STIRVLSGY_A2601.pdf',bbox_inches='tight')
# plt.close()

# PMEL_spliced_A0201
f = '02.25.2025 pMHC A0201 PMEL Splicing 4Rxn(638760831189883078).csv'
dir_path = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/codes/summary/refold/pMHC_refolding_Rawdata'
df = pd.read_csv(os.path.join(dir_path,f),sep='\t',encoding='utf-16',skiprows=2)
fig,ax = plt.subplots()
ax.plot(df['ml.2'].values,df['mAU.1'].values,c='k',lw=1,ls='-')
ax.set_xlabel('Retention Volumn (mL)')
ax.set_ylabel('A280 Intensity (mAU)')
ax.set_title('PMEL_spliced_A0201')
plt.savefig('PMEL_spliced_A0201.pdf',bbox_inches='tight')
plt.close()

# # SLC45A2_spliced_A1101
# f = '07.26.2024 A1101 FTD..(638576074283955483).csv'
# dir_path = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/codes/summary/refold/pMHC_refolding_Rawdata'
# df = pd.read_csv(os.path.join(dir_path,f),sep='\t',encoding='utf-16',skiprows=2)
# fig,ax = plt.subplots()
# ax.plot(df['ml.2'].values,df['mAU.1'].values,c='k',lw=1,ls='-')
# ax.set_xlabel('Retention Volumn (mL)')
# ax.set_ylabel('A280 Intensity (mAU)')
# ax.set_title('SLC45A2_spliced_A1101')
# plt.savefig('SLC45A2_spliced_A1101.pdf',bbox_inches='tight')
# plt.close()

# # CMV_LLD_A0201
# f = '08.19.2025 pMHC A0201 CMV_LLD 4 Rxn(638912016158616656).csv'
# dir_path = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/codes/summary/refold/pMHC_refolding_Rawdata'
# df = pd.read_csv(os.path.join(dir_path,f),sep='\t',encoding='utf-16',skiprows=2)
# fig,ax = plt.subplots()
# ax.plot(df['ml'].values,df['mAU'].values,c='k',lw=1,ls='-')
# ax.set_xlabel('Retention Volumn (mL)')
# ax.set_ylabel('A280 Intensity (mAU)')
# ax.set_title('CMV_LLD_A0201')
# plt.savefig('CMV_LLD_A0201.pdf',bbox_inches='tight')
# plt.close()
