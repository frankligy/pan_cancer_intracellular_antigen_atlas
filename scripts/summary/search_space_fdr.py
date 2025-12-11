#!/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/bin/python3.7

import pandas as pd
import numpy as np
import sys,os
from Bio.SeqIO.FastaIO import SimpleFastaParser
import subprocess
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from ast import literal_eval
from tqdm import tqdm
import pickle

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

cancers = [
    'BRCA',
    'KIRC',
    'COAD',
    'STAD',
    'MESO',
    'LIHC',
    'ESCA',
    'CESC',
    'BLCA',
    'RT',
    'AML',
    'DLBC',
    'GBM',
    'NBL',
    'PAAD',
    'HNSC',
    'OV',
    'LUSC',
    'LUAD',
    'CHOL',
    'SKCM'
]  

n_samples = [
    1118,
    542,
    483,
    412,
    87,
    374,
    185,
    306,
    412,
    63,
    151,
    48,
    170,
    157,
    179,
    522,
    429,
    502,
    541,
    35,
    472
]

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

def classify_source(item):
    if item.startswith('REV__'):
        identity = 'reverse'
    elif item.startswith('CON__'):
        identity = 'contaminant'
    elif 'missense_variant' in item or 'inframe' in item or 'frameshift_variant' in item:
        identity = 'variant'
    elif item.startswith('tr') or item.startswith('sp'):
        identity = 'pathogen'
    elif 'fusion' in item:
        identity = 'fusion' 
    elif 'intron' in item:
        identity = 'intron_retention'
    elif 'nuORF' in item:
        identity = 'nuORF'
    elif item.startswith('chr'):
        if 'TE_info' in item:
            identity = 'TE_chimeric_transcript'
        else:
            identity = 'splicing'
    elif '_dup' in item or item.endswith('sense'):
        identity = 'ERV'
    elif item.startswith('ENSG'):
        identity = 'self_gene'
    elif item.startswith('nc'):
        identity = 'nc_isoform'
    else:
        identity = 'unknown'
    return identity

root_atlas_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas'
immuno_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/immunopeptidome'

iedb = pd.read_csv('all_epitope_no_b_human_linear_mhc_i.tsv',sep='\t')
iedb_peptides = set(iedb['Epitope - Name'].unique())

c = 'BRCA'

# search_space
space_dict = {}
combined_fasta = os.path.join(immuno_dir,cancers2immuno[c],'combined_{}_pan_cancer.fasta'.format(c))
with open(combined_fasta,'r') as in_handle:
    for title,seq in SimpleFastaParser(in_handle):
        identity = classify_source(title)
        if identity in space_dict.keys():
            space_dict[identity] += len(seq)
        else:
            space_dict[identity] = 0

# maxquant fdr iedb and 70
msms1 = pd.read_csv(os.path.join(root_atlas_dir,c,'antigen','0.01','msmsScans_all_add_tesorai.txt'),sep='\t')
msms1 = msms1.loc[msms1['Identified_vanilla']==True,:]

msms5 = pd.read_csv(os.path.join(root_atlas_dir,c,'antigen','0.05','msmsScans_all_add_tesorai.txt'),sep='\t')
msms5 = msms5.loc[msms5['Identified_vanilla']==True,:]

only_in_5 = set(msms5['Sequence']).difference(set(msms1['Sequence']))

iedb_common = iedb_peptides.intersection(only_in_5)
print(len(iedb_common))
seventy1 = msms1.loc[msms1['Score']>70,:]
seventy5 = msms5.loc[msms5['Score']>70,:]
prop_seventy1 = seventy1.shape[0] / msms1.shape[0]
prop_seventy5 = seventy5.shape[0] / msms5.shape[0]
print(prop_seventy1,prop_seventy5)

# rescore fdr iedb and 70
msms1 = pd.read_csv(os.path.join(root_atlas_dir,c,'antigen','0.01','msmsScans_all_add_tesorai.txt'),sep='\t')
msms1 = msms1.loc[msms1['Identified_rescore']==True,:]

msms5 = pd.read_csv(os.path.join(root_atlas_dir,c,'antigen','0.05','msmsScans_all_add_tesorai.txt'),sep='\t')
msms5 = msms5.loc[msms5['Identified_rescore']==True,:]

only_in_5 = set(msms5['Sequence']).difference(set(msms1['Sequence']))

iedb_common = iedb_peptides.intersection(only_in_5)
print(len(iedb_common))
seventy1 = msms1.loc[msms1['Score']>70,:]
seventy5 = msms5.loc[msms5['Score']>70,:]
prop_seventy1 = seventy1.shape[0] / msms1.shape[0]
prop_seventy5 = seventy5.shape[0] / msms5.shape[0]
print(prop_seventy1,prop_seventy5)


