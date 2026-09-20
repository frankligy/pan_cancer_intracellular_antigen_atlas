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
    'PRAD',
    'UCEC'
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

# # raw ms
# root_des_dir = '/gpfs/data/yarmarkovichlab/public/ImmunoVerse/raw_MS_result'
# root_immuno_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/immunopeptidome'
# ## tesorai
# tesorai_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/NYU_Tesorai_all_searches'
# for c in cancers:
#     old_dir = os.getcwd()
#     os.chdir(tesorai_dir)
#     cmd = 'for f in tesorai_peptide_fdr_*.tsv; do echo $f; done | grep "{}"'.format(c)
#     needed_files = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#     os.chdir(old_dir)
#     for tesorai_file in needed_files:
#         tesorai_path = os.path.join(tesorai_dir,tesorai_file)
#         subprocess.run("cp {} {}".format(tesorai_path,os.path.join(root_des_dir,'tesorai')),shell=True)
# ## rescore
# for c in cancers:
#     if not os.path.exists(os.path.join(root_des_dir,'rescore',c)):
#         os.mkdir(os.path.join(root_des_dir,'rescore',c))
#     atlas_dir = os.path.join(root_atlas_dir,c)
#     rescore_dir = os.path.join(atlas_dir,'rescore')
#     old_dir = os.getcwd()
#     os.chdir(rescore_dir)
#     psms = subprocess.run('find . -mindepth 2 -maxdepth 2 -type f -name "ms2rescore.mokapot.psms.txt" -exec echo {} \;',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#     for psm in psms:
#         d = psm.split('/')[1]
#         if not os.path.exists(os.path.join(root_des_dir,'rescore',c,d)):
#             os.mkdir(os.path.join(root_des_dir,'rescore',c,d))
#         subprocess.run("cp {} {}".format(psm,os.path.join(root_des_dir,'rescore',c,d)),shell=True)
#     os.chdir(old_dir)
# ## maxquant
# for c in cancers:
#     if not os.path.exists(os.path.join(root_des_dir,'maxquant',c)):
#         os.mkdir(os.path.join(root_des_dir,'maxquant',c))
#     immuno_dir = os.path.join(root_immuno_dir,cancers2immuno[c])
#     old_dir = os.getcwd()
#     os.chdir(immuno_dir)
#     ### msms.txt
#     msms_all = subprocess.run('find . -mindepth 4 -maxdepth 4 -type f -name "msms.txt" -exec echo {} \;',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#     for msms in msms_all:
#         d = msms.split('/')[1]
#         if not os.path.exists(os.path.join(root_des_dir,'maxquant',c,d)):
#             os.mkdir(os.path.join(root_des_dir,'maxquant',c,d))
#         subprocess.run("cp {} {}".format(msms,os.path.join(root_des_dir,'maxquant',c,d)),shell=True)
#     ### msmsScans.txt
#     msms_all = subprocess.run('find . -mindepth 4 -maxdepth 4 -type f -name "msmsScans.txt" -exec echo {} \;',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#     for msms in msms_all:
#         d = msms.split('/')[1]
#         if not os.path.exists(os.path.join(root_des_dir,'maxquant',c,d)):
#             os.mkdir(os.path.join(root_des_dir,'maxquant',c,d))
#         subprocess.run("cp {} {}".format(msms,os.path.join(root_des_dir,'maxquant',c,d)),shell=True)
#     os.chdir(old_dir)
#     ### evidence.txt
#     msms_all = subprocess.run('find . -mindepth 4 -maxdepth 4 -type f -name "evidence.txt" -exec echo {} \;',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#     for msms in msms_all:
#         d = msms.split('/')[1]
#         if not os.path.exists(os.path.join(root_des_dir,'maxquant',c,d)):
#             os.mkdir(os.path.join(root_des_dir,'maxquant',c,d))
#         subprocess.run("cp {} {}".format(msms,os.path.join(root_des_dir,'maxquant',c,d)),shell=True)
#     os.chdir(old_dir)
# ## tesorai normal
# source = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/NYU_Tesorai_all_searches/tesorai_peptide_fdr_normal.tsv'
# destination = '/gpfs/data/yarmarkovichlab/public/ImmunoVerse/raw_MS_result/HLA_ligand_atlas/tesorai'
# subprocess.run('cp {} {}'.format(source,destination),shell=True)
# ## maxquant normal
# result_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/immuno'
# old_dir = os.getcwd()
# os.chdir(result_dir)
# all_tissues = subprocess.run("for f in *; do echo $f; done",shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
# os.chdir(old_dir)
# root_des_dir_normal = '/gpfs/data/yarmarkovichlab/public/ImmunoVerse/raw_MS_result/HLA_ligand_atlas'
# for t in all_tissues:
#     if not os.path.exists(os.path.join(root_des_dir_normal,'maxquant',t)):
#         os.mkdir(os.path.join(root_des_dir_normal,'maxquant',t))
#     immuno_dir = os.path.join(result_dir,t)
#     old_dir = os.getcwd()
#     os.chdir(immuno_dir)
#     ### msms.txt
#     msms_all = subprocess.run('find . -mindepth 4 -maxdepth 4 -type f -name "msms.txt" -exec echo {} \;',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#     for msms in msms_all:
#         d = msms.split('/')[1]
#         if not os.path.exists(os.path.join(root_des_dir_normal,'maxquant',t,d)):
#             os.mkdir(os.path.join(root_des_dir_normal,'maxquant',t,d))
#         subprocess.run("cp {} {}".format(msms,os.path.join(root_des_dir_normal,'maxquant',t,d)),shell=True)
#     ### msmsScans.txt
#     msms_all = subprocess.run('find . -mindepth 4 -maxdepth 4 -type f -name "msmsScans.txt" -exec echo {} \;',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#     for msms in msms_all:
#         d = msms.split('/')[1]
#         if not os.path.exists(os.path.join(root_des_dir_normal,'maxquant',t,d)):
#             os.mkdir(os.path.join(root_des_dir_normal,'maxquant',t,d))
#         subprocess.run("cp {} {}".format(msms,os.path.join(root_des_dir_normal,'maxquant',t,d)),shell=True)
#     os.chdir(old_dir)
#     ### evidence.txt
#     msms_all = subprocess.run('find . -mindepth 4 -maxdepth 4 -type f -name "evidence.txt" -exec echo {} \;',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#     for msms in msms_all:
#         d = msms.split('/')[1]
#         if not os.path.exists(os.path.join(root_des_dir_normal,'maxquant',t,d)):
#             os.mkdir(os.path.join(root_des_dir_normal,'maxquant',t,d))
#         subprocess.run("cp {} {}".format(msms,os.path.join(root_des_dir_normal,'maxquant',t,d)),shell=True)
#     os.chdir(old_dir)


# # consolidated ms
# root_des_dir = '/gpfs/data/yarmarkovichlab/public/ImmunoVerse/raw_MS_result'
# root_atlas_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas'
# for c in cancers:
#     if not os.path.exists(os.path.join(root_des_dir,'consolidated',c)):
#         os.mkdir(os.path.join(root_des_dir,'consolidated',c))
#     df_path = os.path.join(root_atlas_dir,c,'antigen','0.01','msmsScans_all_add_tesorai.txt')
#     cmd = 'cp {} {}'.format(df_path,os.path.join(root_des_dir,'consolidated',c))
#     subprocess.run(cmd,shell=True)


# molecular_catalogue, I'll manually add mutation
variant_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/variants'
gene_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/gene'
custom_dic = {
    'BRCA':['TCGA-BRCA.mutect2_snv.tsv','gProfiler_hsapiens_brca.csv','manifest_BRCA_1118.tsv','TCGA-BRCA.htseq_counts.tsv'],
    'KIRC':['TCGA-KIRC.mutect2_snv.tsv','gProfiler_hsapiens_kirc.csv','manifest_kirc_542.tsv','TCGA-KIRC.htseq_counts.tsv'],
    'COAD':['TCGA-COAD.mutect2_snv.tsv','gProfiler_hsapiens_coad.csv','manifest_coad_524.tsv','TCGA-COAD.htseq_counts.tsv'],
    'STAD':['TCGA-STAD.somaticmutation_wxs.tsv','gProfiler_hsapiens_stad.csv','manifest_stad_448.tsv','TCGA-STAD.star_counts.tsv'],
    'MESO':['TCGA-MESO.mutect2_snv.tsv','gProfiler_hsapiens_meso.csv','manifest_MESO_87.tsv','TCGA-MESO.htseq_counts.tsv'],
    'LIHC':['TCGA-LIHC.mutect2_snv.tsv','gProfiler_hsapiens_lihc.csv','manifest_LIHC_424.tsv','TCGA-LIHC.htseq_counts.tsv'],
    'ESCA':['TCGA-ESCA.somaticmutation_wxs.tsv','gProfiler_hsapiens_esca.csv','manifest_esca_198.tsv','TCGA-ESCA.star_counts.tsv'],
    'CESC':['TCGA-CESC.mutect2_snv.tsv','gProfiler_hsapiens_cesc.csv','manifest_CESC_309.tsv','TCGA-CESC.htseq_counts.tsv'],
    'BLCA':['TCGA-BLCA.somaticmutation_wxs.tsv','gProfiler_hsapiens_blca.csv','manifest_blca_431.tsv','TCGA-BLCA.star_counts.tsv'],
    'RT':[None,None,'manifest_rt_69.tsv','TARGET-RT.star_counts.tsv'],
    'AML':['TCGA-LAML.mutect2_snv.tsv','gProfiler_hsapiens_aml.csv','manifest_AML_151.tsv','TARGET-AML.star_counts.tsv'],
    'DLBC':['TCGA-DLBC.mutect2_snv.tsv','gProfiler_hsapiens_dlbc.csv','manifest_DLBC_48.tsv','TCGA-DLBC.htseq_counts.tsv'],
    'GBM':['TCGA-GBM.mutect2_snv.tsv','gProfiler_hsapiens_gbm.csv','manifest_GBM_175.tsv','TCGA-GBM.htseq_counts.tsv'],
    'NBL':['TARGET-NBL.mutect2_snv.tsv','gProfiler_hsapiens_nbl.csv','manifest_NBL_161.tsv','TARGET-NBL.htseq_counts.tsv'],
    'PAAD':['TCGA-PAAD.mutect2_snv.tsv','gProfiler_hsapiens_paad.csv','manifest_PAAD_183.tsv','TCGA-PAAD.htseq_counts.tsv'],
    'HNSC':['TCGA-HNSC.mutect2_snv.tsv','gProfiler_hsapiens_hnsc.csv','manifest_HNSC_566.tsv','TCGA-HNSC.htseq_counts.tsv'],
    'OV':['TCGA-OV.mutect2_snv.tsv','gProfiler_hsapiens_ov.csv','manifest_ov_429.tsv','TCGA-OV.htseq_counts.tsv'],
    'LUSC':['TCGA-LUSC.somaticmutation_wxs.tsv','gProfiler_hsapiens_lusc.csv','manifest_lusc_553.tsv','TCGA-LUSC.star_counts.tsv'],
    'LUAD':['TCGA-LUAD.mutect2_snv.tsv','gProfiler_hsapiens_luad.csv','manifest_LUAD_541.tsv','TCGA-LUAD.htseq_counts.tsv'],
    'CHOL':['TCGA-CHOL.somaticmutation_wxs.tsv','gProfiler_hsapiens_chol.csv','manifest_chol_44.tsv','TCGA-CHOL.star_counts.tsv'],
    'SKCM':['TCGA-SKCM.mutect2_snv.tsv','gProfiler_hsapiens_skcm.csv','manifest_SKCM_472.tsv','TCGA-SKCM.htseq_counts.tsv'],
    'ALL_BALL_P1':[None,None,'manifest_all_ball_p1.tsv','TARGET-ALL-P1.star_counts.tsv'],
    'KICH':['TCGA-KICH.somaticmutation_wxs.tsv','gProfiler_hsapiens_KICH.csv','manifest_KICH-91.tsv','TCGA-KICH.star_counts.tsv'],
    'WT':['TARGET-WT.somaticmutation_wxs.tsv','gProfiler_hsapiens_wt.csv','manifest_WT_136.tsv','TARGET-WT.star_counts.tsv'],
    'CCSK':[None,None,'manifest_ccsk_13.tsv','TARGET-CCSK.star_counts.tsv'],
    'KIRP':['TCGA-KIRP.somaticmutation_wxs.tsv','gProfiler_hsapiens_KIRP.csv','manifest_KIRP-323.tsv','TCGA-KIRP.star_counts.tsv'],
    'LGG':['TCGA-LGG.somaticmutation_wxs.tsv','gProfiler_hsapiens_LGG.csv','manifest_LGG_534.tsv','TCGA-LGG.star_counts.tsv'],
    'READ':['TCGA-READ.somaticmutation_wxs.tsv','gProfiler_hsapiens_READ.csv','manifest_READ_177.tsv','TCGA-READ.star_counts.tsv'],
    'UCS':['TCGA-UCS.somaticmutation_wxs.tsv','gProfiler_hsapiens_UCS.csv','manifest_UCS_57.tsv','TCGA-UCS.star_counts.tsv'],
    'SARC':['TCGA-SARC.somaticmutation_wxs.tsv','gProfiler_hsapiens_sarc.csv','manifest_SARC_265.tsv','TCGA-SARC.star_counts.tsv'],
    'TGCT':['TCGA-TGCT.somaticmutation_wxs.tsv','gProfiler_hsapiens_TGCT.csv','manifest_TGCT_156.tsv','TCGA-TGCT.star_counts.tsv'],
    'THYM':['TCGA-THYM.somaticmutation_wxs.tsv','gProfiler_hsapiens_THYM.csv','manifest_THYM_122.tsv','TCGA-THYM.star_counts.tsv'],
    'THCA':['TCGA-THCA.somaticmutation_wxs.tsv','gProfiler_hsapiens_thca.csv','manifest_THCA_572.tsv','TCGA-THCA.star_counts.tsv'],
    'UVM':['TCGA-UVM.somaticmutation_wxs.tsv','gProfiler_hsapiens_uvm.csv','manifest_UVM_80.tsv','TCGA-UVM.star_counts.tsv'],
    'OS':['TARGET-OS.mutect2_snv.tsv','gProfiler_hsapiens_os.csv','manifest_OS_88.tsv','TARGET-OS.htseq_counts.tsv'],
    'PRAD':['TCGA-PRAD.mutect2_snv.tsv','gProfiler_hsapiens_prad.csv','manifest_PRAD_554.tsv','TCGA-PRAD.htseq_counts.tsv'],
    'UCEC':['TCGA-UCEC.somaticmutation_wxs.tsv','gProfiler_hsapiens_ucec.csv','manifest_UCEC_589.tsv','TCGA-UCEC.star_counts.tsv']
}

root_des_dir = '/gpfs/data/yarmarkovichlab/public/ImmunoVerse/molecular_catalogue'
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
    fusion_all = os.path.join(atlas_dir,'fusion.txt')
    fusion_recurrent = os.path.join(atlas_dir,'fusion_recurrent.txt')
    te_all = os.path.join(atlas_dir,'tumor_erv.h5ad')
    te_rec = os.path.join(atlas_dir,'ERV.txt')
    te_good = os.path.join(atlas_dir,'good_erv.txt')
    manifest = os.path.join(atlas_dir,custom_dic[c][2])
    gene_count = os.path.join(gene_dir,custom_dic[c][3])
    real_common = os.path.join(atlas_dir,'real_common.txt')
    real_common_membrane = os.path.join(atlas_dir,'real_common_membrane.txt')
    if not c in ['NBL','OS','RT','CCSK','WT','ALL_BALL_P1']:
        hla_types = os.path.join(atlas_dir,'hla_types.txt')
        manifest_fusion = None
    else:
        hla_types = None
        manifest_fusion = os.path.join(atlas_dir,'manifest_fusion.tsv')
    if not c in ['NBL','OS','RT','CCSK','WT','ALL_BALL_P1','CHOL','MESO','UVM','SARC']:
        deg_file = os.path.join(atlas_dir,'deg.txt')
    else:
        deg_file = None
    if custom_dic[c][0] is None:
        raw_mutation = None
        map_mutation = None
        mutation_rec = None
    else:
        raw_mutation = os.path.join(variant_dir,custom_dic[c][0])
        map_mutation = os.path.join(variant_dir,custom_dic[c][1])
        mutation_rec = os.path.join(atlas_dir,'mutation_rec.txt')

    need_to_cp = [gene_tpm,gene_lfc,pathogen_all,pathogen_rec,intron_all,intron_peptide_all,intron_rec,splicing_all,splicing_rec,
                  fusion_all,fusion_recurrent, manifest,gene_count, real_common, real_common_membrane, hla_types, manifest_fusion, deg_file, mutation_rec, raw_mutation, map_mutation, te_all, te_rec, te_good]

    des_dir = os.path.join(root_des_dir,c)
    if not os.path.exists(des_dir):
        os.mkdir(des_dir)

    for f in need_to_cp:
        if f is not None:
            subprocess.run('cp {} {}'.format(f,des_dir),shell=True)

# search space tesorai
root_des_dir = '/gpfs/data/yarmarkovichlab/public/ImmunoVerse/search_space_tesorai'
for c in cancers:
    db_fasta_dir = os.path.join(root_atlas_dir,c,'db_fasta_tesorai')
    des_dir = os.path.join(root_des_dir,c)
    if not os.path.exists(des_dir):
        os.mkdir(des_dir)
    subprocess.run('cp -R {} {}'.format(db_fasta_dir,des_dir),shell=True)

# search space nt
root_des_dir = '/gpfs/data/yarmarkovichlab/public/ImmunoVerse/search_space_nt'
for c in cancers:
    db_fasta_dir = os.path.join(root_atlas_dir,c,'db_fasta_tesorai_nt')
    des_dir = os.path.join(root_des_dir,c)
    if not os.path.exists(des_dir):
        os.mkdir(des_dir)
    ori_dir = os.getcwd()
    os.chdir(db_fasta_dir)
    all_fs = subprocess.run("for file in *.txt; do echo $file; done",shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
    os.chdir(ori_dir)
    for f in all_fs:
        subprocess.run('cp {} {}'.format(os.path.join(db_fasta_dir,f),des_dir),shell=True)


