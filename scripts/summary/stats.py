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
import anndata as ad

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


def run_self_gene_de(ensgs,cancer):

    final = pd.read_csv(os.path.join(root_atlas_dir,cancer,'gene_tpm.txt'),sep='\t',index_col=0)
    final = final.loc[~final.index.duplicated(),:]
    tumor_expr = final.loc[ensgs,:].values
    tumor_label = ['tumor{}'.format(i+1) for i in range(tumor_expr.shape[1])]

    gtex = pd.read_csv(os.path.join(database_dir,'bulk-gex_v8_rna-seq_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct'),sep='\t',skiprows=2,index_col=0)
    cond = ~gtex.columns.isin(['Cells - EBV-transformed lymphocytes','Cells - Cultured fibroblasts','Testis'])
    gtex = gtex.loc[:,cond]
    gtex.index = [item.split('.')[0] for item in gtex.index]
    ensg2symbol = pd.Series(index=gtex.index.tolist(),data=gtex['Description'].tolist()).to_dict()
    series = gtex.loc[[ensgs[0]],:].iloc[0,:].iloc[1:]

    adata = ad.read_h5ad(os.path.join(database_dir,'gtex_gene_all.h5ad'))  # 56200 × 17382
    adata.obs_names = [item.split('.')[0] for item in adata.obs_names]
    adata.obs_names_make_unique()
    adata_gene = adata[ensgs,:]
    normal_expr_list = []
    normal_expr_list_label = []
    for t in series.index.tolist():
        values = adata_gene[:,adata_gene.var['tissue']==t].X.toarray()
        normal_label = ['{}{}'.format(t,i+1) for i in range(values.shape[1])]
        normal_expr_list.append(values)
        normal_expr_list_label.extend(normal_label)

    all_array = [tumor_expr] + normal_expr_list
    df_data = np.concatenate(all_array,axis=1)
    df = pd.DataFrame(data=df_data,columns=tumor_label+normal_expr_list_label,index=ensgs)
    df.to_csv('exp.original-steady-state.txt',sep='\t')
    with open('groups.txt','w') as f:
        for item in tumor_label:
            f.write('{}\t1\ttumor\n'.format(item))
        for item in normal_expr_list_label:
            f.write('{}\t2\tnormal\n'.format(item))
    with open('comps.txt','w') as f:
        f.write('1\t2\n')

    # run altanalyze, assuming singularity is enabled
    os.makedirs('./diff_dir/altanalyze_output/ExpressionInput')
    os.rename('exp.original-steady-state.txt','./diff_dir/altanalyze_output/ExpressionInput/exp.original-steady-state.txt')
    os.rename('groups.txt','./diff_dir/groups.txt')
    os.rename('comps.txt','./diff_dir/comps.txt')
    old_dir = os.getcwd()
    os.chdir('./diff_dir')
    writable_path = '/gpfs/data/yarmarkovichlab/Frank/SNAF_ecosystem/altanalyze'
    cmd = 'singularity run -B $PWD:/mnt --writable {} DE altanalyze_output groups.txt'.format(writable_path)
    subprocess.run(cmd,shell=True)
    os.chdir(old_dir)
    os.rename('./diff_dir/altanalyze_output/ExpressionInput/DEGs-LogFold_0.0_adjp/GE.tumor_vs_normal.txt','DE_result_{}.txt'.format(cancer))
    subprocess.run('rm -rf ./diff_dir',shell=True)

    


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
database_dir = '/gpfs/data/yarmarkovichlab/public/ImmunoVerse/database'
bayests_xy_path = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/gene/full_results_XY_essential_tissues.txt'

# self_gene 
final = pd.read_csv('../final_all_ts_antigens.txt',sep='\t')
final = final.loc[final['typ']=='self_gene',:]
total_ensg_intra = list(set(final['ensgs'].values))
mem = pd.read_csv('../gene_morpheus_mem.txt',sep='\t',index_col=[0,1])
total_ensg_mem = mem.index.to_frame(index=False)[0][1:].values.tolist()
total_ensg = list(set(total_ensg_intra + total_ensg_mem))

# for cancer in cancers:
#     run_self_gene_de(total_ensg,cancer)

col1 = []
col2 = []
for row in final.itertuples():
    if row.typ == 'self_gene':
        c = row.cancer
        ensg = row.ensgs
        df = pd.read_csv('DE_result_{}.txt'.format(c),sep='\t',index_col=0)
        rawp = df.loc[ensg,:]['rawp']
        adjp = df.loc[ensg,:]['adjp']
        col1.append(rawp)
        col2.append(adjp)
    else:
        col1.append(None)
        col2.append(None)
final['ts_rawp'] = col1
final['ts_adjp'] = col2
final.to_csv('check.txt',sep='\t',index=None)