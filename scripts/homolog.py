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
    412,
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
    541,
    502,
    35,
    472
]

root_atlas_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas'
bayests_xy_path = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/gene/full_results_XY_essential_tissues.txt'
gtex_median_path = '/gpfs/data/yarmarkovichlab/chordoma/NeoVerse_analysis/bulk-gex_v8_rna-seq_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct'
membrane_path = '/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/NeoVerse/data/human_membrane_proteins_acc2ens.txt'

def process_gtex(gtex_median_path):
    gtex = pd.read_csv(gtex_median_path,sep='\t',skiprows=2,index_col=0)
    gtex = gtex.loc[~gtex.index.str.contains('_PAR'),:]
    gtex.index = [item.split('.')[0] for item in gtex.index]
    gtex.drop(labels=['Cells - Cultured fibroblasts', 'Cells - EBV-transformed lymphocytes'],axis=1,inplace=True)
    gtex.drop(labels=['Ovary','Prostate','Testis','Vagina','Adrenal Gland','Cervix - Endocervix','Cervix - Ectocervix','Fallopian Tube','Pituitary'],axis=1,inplace=True)
    ensg2medians = {}
    for ensg in gtex.index:
        medians = gtex.loc[ensg,:].iloc[1:].values.tolist()
        ensg2medians[ensg] = medians
    all_tissues = gtex.columns[1:].values.tolist()
    ensg2symbol = pd.Series(index=gtex.index.tolist(),data=gtex['Description'].tolist()).to_dict()
    return ensg2medians,all_tissues,ensg2symbol

def process_tumor_gene():
    dic = {}
    for c in cancers:
        gene_lfc_path = os.path.join(root_atlas_dir,c,'gene_lfc.txt')
        gene_lfc = pd.read_csv(gene_lfc_path,sep='\t',index_col=0)
        c_dic = gene_lfc['median_tumor'].to_dict()
        dic[c] = c_dic
    return dic
        
    

data = []
for c in cancers:
    final_path = os.path.join(root_atlas_dir,c,'antigen','fdr','final_enhanced.txt')
    final = pd.read_csv(final_path,sep='\t')
    cond = [False if '[]' in item else True for item in final['presented_by_each_sample_hla']]
    final = final.loc[cond,:]
    final = final.loc[final['typ']=='self_gene',:]
    final = final.loc[final['unique']!=True,:]
    data.append(final)
final = pd.concat(data,axis=0,keys=cancers).reset_index(level=-2).rename(columns={'level_0':'cancer'})

final.to_csv('homolog_self_gene.txt',sep='\t',index=None)


reported_pairs = [
    ('AGR2','AGR3'),
    ('ASCL1','ASCL2'),
    ('ANXA8','ANXA8L1'),
    ('ASF1A','ASF1B'),
    ('AURKA','AURKB'),
    ('ENPP2','ENPP1','ENPP3')
]


ensg2medians,all_tissues,ensg2symbol = process_gtex(gtex_median_path)
all_genes = []
for item in final['ensgs']:
    item = literal_eval(item)
    all_genes.extend(item)
all_genes = list(set(all_genes))


dic = process_tumor_gene()
ensg2tumors = {}
for gene in all_genes:
    data = []
    for k,v in dic.items():
        data.append(v[gene])
    ensg2tumors[gene] = data


df = pd.read_csv(membrane_path,sep='\t')
membrane_ensg = set(df['Ens'].tolist())

all_data = []
for gene in all_genes:
    data = []
    data.extend(ensg2tumors[gene])
    data.extend(ensg2medians[gene])
    all_data.append(data)
df = pd.DataFrame.from_records(all_data,columns=cancers+all_tissues,index=all_genes)

ori_array = [tuple(['cancer']*21+['normal']*43),tuple(df.columns.tolist())]
mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
df.columns = mi

ori_array = [tuple(df.index.tolist()),tuple([ensg2symbol[item] for item in df.index]),tuple([True if item in membrane_ensg else False for item in df.index])]
mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
df.index = mi
df.to_csv('homolog_gene_morpheus.txt',sep='\t')
    
    





    



    