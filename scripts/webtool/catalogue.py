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