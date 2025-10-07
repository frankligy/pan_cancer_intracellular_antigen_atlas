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
import math
import argparse
import json
from copy import deepcopy

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

root_atlas_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas'
root_des_dir = '/gpfs/data/yarmarkovichlab/public/ImmunoVerse/molecular_catalogue'
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

added_cancer = [
    'ALL_BALL_P1',
    'KICH',
    'WT',
    'CCSK',
    'KIRP',
    'LGG',
    'READ',
    'UCS',
    'SARC',
    'TGCT',
    'THYM',
    'THCA',
    'UVM',
    'OS',
    'PRAD'
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


# raw ms
root_des_dir = '/gpfs/data/yarmarkovichlab/public/ImmunoVerse/raw_MS_result'
root_immuno_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/immunopeptidome'
## tesorai
tesorai_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/NYU_Tesorai_all_searches'
for c in cancers:
    tesorai_path = os.path.join(tesorai_dir,'tesorai_peptide_fdr_{}.tsv'.format(c))
    subprocess.run("cp {} {}".format(tesorai_path,os.path.join(root_des_dir,'tesorai')),shell=True)
## rescore
for c in cancers:
    if not os.path.exists(os.path.join(root_des_dir,'rescore',c)):
        os.mkdir(os.path.join(root_des_dir,'rescore',c))
    atlas_dir = os.path.join(root_atlas_dir,c)
    rescore_dir = os.path.join(atlas_dir,'rescore')
    old_dir = os.getcwd()
    os.chdir(rescore_dir)
    psms = subprocess.run('find . -mindepth 2 -maxdepth 2 -type f -name "ms2rescore.mokapot.psms.txt" -exec echo {} \;',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
    for psm in psms:
        d = psm.split('/')[1]
        if not os.path.exists(os.path.join(root_des_dir,'rescore',c,d)):
            os.mkdir(os.path.join(root_des_dir,'rescore',c,d))
        subprocess.run("cp {} {}".format(psm,os.path.join(root_des_dir,'rescore',c,d)),shell=True)
    os.chdir(old_dir)
## maxquant
for c in cancers:
    if not os.path.exists(os.path.join(root_des_dir,'maxquant',c)):
        os.mkdir(os.path.join(root_des_dir,'maxquant',c))
    immuno_dir = os.path.join(root_immuno_dir,cancers2immuno[c])
    old_dir = os.getcwd()
    os.chdir(immuno_dir)
    ### msms.txt
    msms_all = subprocess.run('find . -mindepth 4 -maxdepth 4 -type f -name "msms.txt" -exec echo {} \;',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
    for msms in msms_all:
        d = msms.split('/')[1]
        if not os.path.exists(os.path.join(root_des_dir,'maxquant',c,d)):
            os.mkdir(os.path.join(root_des_dir,'maxquant',c,d))
        subprocess.run("cp {} {}".format(msms,os.path.join(root_des_dir,'maxquant',c,d)),shell=True)
    ### msmsScans.txt
    msms_all = subprocess.run('find . -mindepth 4 -maxdepth 4 -type f -name "msmsScans.txt" -exec echo {} \;',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
    for msms in msms_all:
        d = msms.split('/')[1]
        if not os.path.exists(os.path.join(root_des_dir,'maxquant',c,d)):
            os.mkdir(os.path.join(root_des_dir,'maxquant',c,d))
        subprocess.run("cp {} {}".format(msms,os.path.join(root_des_dir,'maxquant',c,d)),shell=True)
    os.chdir(old_dir)
## tesorai normal
source = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/NYU_Tesorai_all_searches/tesorai_peptide_fdr_healthy.tsv'
destination = '/gpfs/data/yarmarkovichlab/public/ImmunoVerse/raw_MS_result/HLA_ligand_atlas/tesorai'
subprocess.run('cp {} {}'.format(source,destination),shell=True)
## maxquant normal
result_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/immuno'
old_dir = os.getcwd()
os.chdir(result_dir)
all_tissues = subprocess.run("for f in *; do echo $f; done",shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
os.chdir(old_dir)
root_des_dir_normal = '/gpfs/data/yarmarkovichlab/public/ImmunoVerse/raw_MS_result/HLA_ligand_atlas'
for t in all_tissues:
    if not os.path.exists(os.path.join(root_des_dir_normal,'maxquant',t)):
        os.mkdir(os.path.join(root_des_dir_normal,'maxquant',t))
    immuno_dir = os.path.join(result_dir,t)
    old_dir = os.getcwd()
    os.chdir(immuno_dir)
    ### msms.txt
    msms_all = subprocess.run('find . -mindepth 4 -maxdepth 4 -type f -name "msms.txt" -exec echo {} \;',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
    for msms in msms_all:
        d = msms.split('/')[1]
        if not os.path.exists(os.path.join(root_des_dir_normal,'maxquant',t,d)):
            os.mkdir(os.path.join(root_des_dir_normal,'maxquant',t,d))
        subprocess.run("cp {} {}".format(msms,os.path.join(root_des_dir_normal,'maxquant',t,d)),shell=True)
    ### msmsScans.txt
    msms_all = subprocess.run('find . -mindepth 4 -maxdepth 4 -type f -name "msmsScans.txt" -exec echo {} \;',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
    for msms in msms_all:
        d = msms.split('/')[1]
        if not os.path.exists(os.path.join(root_des_dir_normal,'maxquant',t,d)):
            os.mkdir(os.path.join(root_des_dir_normal,'maxquant',t,d))
        subprocess.run("cp {} {}".format(msms,os.path.join(root_des_dir_normal,'maxquant',t,d)),shell=True)
    os.chdir(old_dir)

# consolidated ms
root_des_dir = '/gpfs/data/yarmarkovichlab/public/ImmunoVerse/raw_MS_result'
root_atlas_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas'
for c in cancers:
    if not os.path.exists(os.path.join(root_des_dir,'consolidated',c)):
        os.mkdir(os.path.join(root_des_dir,'consolidated',c))
    df_path = os.path.join(root_atlas_dir,c,'antigen','fdr','msmsScans_all_add_tesorai.txt')
    cmd = 'cp {} {}'.format(df_path,os.path.join(root_des_dir,'consolidated',c))
    subprocess.run(cmd,shell=True)








# molecular_catalogue, I'll manually add mutation
for c in cancers + added_cancer:
    atlas_dir = os.path.join(root_atlas_dir,c)
    gene_tpm = os.path.join(atlas_dir,'gene_tpm.txt')
    gene_lfc = os.path.join(atlas_dir,'gene_lfc.txt')
    pathogen_all = os.path.join(atlas_dir,'pathogen_all.txt')
    pathogen_rec = os.path.join(atlas_dir,'pathogen_rec.txt')
    intron_all = os.path.join(atlas_dir,'intron_all.txt')
    intron_peptide_all = os.path.join(atlas_dir,'intron_peptide_all.txt')
    intron_rec = os.path.join(atlas_dir,'intron_rec.txt')
    splicing_all = os.path.join(atlas_dir,'splicing_all.txt')
    splicing_rec = os.path.join(atlas_dir,'splicing_rec.txt')
    if not c in ['NBL','OS','RT','CCSK','WT','ALL_BALL_P1']:
        fusion_all = os.path.join(atlas_dir,'fusion.txt')
    else:
        fusion_all = os.path.join(atlas_dir,'fusion_all.txt')
    if not c in ['NBL','OS','RT','CCSK','WT','ALL_BALL_P1']:
        fusion_recurrent = os.path.join(atlas_dir,'fusion_recurrent.txt')
    else:
        fusion_recurrent = os.path.join(atlas_dir,'fusion_rec.txt')
    if not c in ['NBL','OS','RT','CCSK','WT','ALL_BALL_P1']:
        hla_types = os.path.join(atlas_dir,'hla_types.txt')
    else:
        hla_types = None
    if not c in ['RT','CCSK','ALL_BALL_P1']:
        mutation_rec = os.path.join(atlas_dir,'mutation_rec.txt')
    else:
        mutation_rec = None
    te_all = os.path.join(atlas_dir,'tumor_erv.h5ad')
    te_rec = os.path.join(atlas_dir,'ERV.txt')
    te_good = os.path.join(atlas_dir,'good_erv.txt')

    need_to_cp = [gene_tpm,gene_lfc,pathogen_all,pathogen_rec,intron_all,intron_peptide_all,intron_rec,splicing_all,splicing_rec,
                  fusion_all,fusion_recurrent, hla_types, mutation_rec, te_all, te_rec, te_good]

    des_dir = os.path.join(root_des_dir,c)
    if not os.path.exists(des_dir):
        os.mkdir(des_dir)

    for f in need_to_cp:
        if f is not None:
            subprocess.run('cp {} {}'.format(f,des_dir),shell=True)

# search space
root_des_dir = '/gpfs/data/yarmarkovichlab/public/ImmunoVerse/search_space'
for c in cancers + added_cancer:
    db_fasta_dir = os.path.join(root_atlas_dir,c,'db_fasta')
    des_dir = os.path.join(root_des_dir,c)
    if not os.path.exists(des_dir):
        os.mkdir(des_dir)
    subprocess.run('cp -R {} {}'.format(db_fasta_dir,des_dir),shell=True)


