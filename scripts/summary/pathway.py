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
import anndata as ad
from scipy.sparse import csr_matrix

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

root_atlas_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas'
gtex_median_path = '/gpfs/data/yarmarkovichlab/chordoma/NeoVerse_analysis/bulk-gex_v8_rna-seq_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct'

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

# looking at RT and DLBC peptide
final = pd.read_csv(os.path.join(root_atlas_dir,'DLBC','antigen','fdr','final_enhanced.txt'),sep='\t')
cond = [False if ('[]' in item) and ('(\'HLA-' not in item) else True for item in final['presented_by_each_sample_hla']]
df1 = final.loc[cond,:]
final = pd.read_csv(os.path.join(root_atlas_dir,'RT','antigen','fdr','final_enhanced.txt'),sep='\t')
cond = [False if ('[]' in item) and ('(\'HLA-' not in item) else True for item in final['presented_by_each_sample_hla']]
df2 = final.loc[cond,:]

from scipy.stats import chi2_contingency
mat = []
for df in [df1,df2]:
    cond = [True if pep.endswith('D') or pep.endswith('R') else False for pep in df['pep']]
    n = sum(cond)
    t = df.shape[0]
    mat.append((n,t-n))
    p = n / t
    print(p)

mat = np.array(mat)
res = chi2_contingency(mat)
print(res)

mat_df = pd.concat([df1,df2],axis=0,keys=['DLBC','RT']).reset_index(level=-2).rename(columns={'level_0':'cancer'})
mat_df.to_csv('pathway_immunoproteosome_df.txt',sep='\t',index=None)

mat_df = pd.DataFrame(data=mat,index=['DLBC','RT'],columns=['acidic','non_acidic'])
mat_df.to_csv('pathway_immunoproteosome_stat.txt',sep='\t')
sys.exit('stop')

genes = [
    'PSMB8',
    'PSMB9',
    'PSMB10',
    'TAP1',
    'TAP2',
    'ERAP1',
    'ERAP2',
    'TAPBP',
    'TAPBPL',
    'CALR',
    'CANX',
    'PDIA3',
    # 'NLRC5',
    'B2M',
    'HLA-A',
    'HLA-B',
    'HLA-C',
    # 'HLA-E',
    # 'HLA-G',
    # 'CIITA',
    # 'HLA-DPA1',
    # 'HLA-DPB1',
    # 'HLA-DQA1',
    # 'HLA-DQB1',
    # 'HLA-DRA',
    # 'HLA-DRB1',
    'CXCL8',
    'CXCL9',
    'CXCL10',
    'NFKB1',
    'IFNGR1',
    'IFNGR2',
    'JAK1',
    'JAK2',
    'TYK2',
    'STAT1',
    'STAT2',
    'STAT3'
]


ensg2medians,all_tissues,ensg2symbol = process_gtex(gtex_median_path)
symbol2ensg = {v:k for k,v in ensg2symbol.items()}
all_genes = [symbol2ensg[item] for item in genes]
dic = process_tumor_gene()
ensg2tumors = {}
for gene in all_genes:
    data = []
    for k,v in dic.items():
        data.append(v[gene])
    ensg2tumors[gene] = data

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

ori_array = [tuple(df.index.tolist()),
             tuple([ensg2symbol[item] for item in df.index])]
mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
df.index = mi
df.to_csv('gene_pathway.txt',sep='\t')
    