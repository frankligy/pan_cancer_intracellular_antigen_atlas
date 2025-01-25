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

data = []
for c in cancers:
    final_path = os.path.join(root_atlas_dir,c,'antigen','fdr','final_enhanced.txt')
    final = pd.read_csv(final_path,sep='\t')
    cond = [False if ('[]' in item) and ('(\'HLA-' not in item) else True for item in final['presented_by_each_sample_hla']]
    final = final.loc[cond,:]
    final = final.loc[final['typ']=='pathogen',:]
    data.append(final)
final = pd.concat(data,axis=0,keys=cancers).reset_index(level=-2).rename(columns={'level_0':'cancer'})

col = []
for item in final['source']:
    if 'Humancytomegalovirus' in item:
        col.append('CMV')
    elif 'Fusobacteriumnucleatum' in item:
        col.append('F.Nucleatum')
    elif 'Epstein-Barrvirus' in item:
        col.append('EBV')
    elif 'Humanpapillomavirus' in item:
        col.append('HPV')
    elif 'HepatitisBvirus' in item:
        col.append('HBV')
    elif 'Helicobacterpylori' in item:
        col.append('H.Pylori')
    elif 'Nialliacirculans' in item:
        col.append('N.Circulans')
    elif 'Clostridiumintestinale' in item:
        col.append('C.Intestinale')
    elif 'Campylobacterureolyticus' in item:
        col.append('C.Ureolyticus')
    else:
        col.append('unknown')
final['strain'] = col
final = final.loc[final['strain']!='unknown',:]
final.to_csv('all_pathogen.txt',sep='\t',index=None)


# # look for N circulans gene in ov
# final_n = final.loc[(final['strain']=='N.Circulans') & (final['cancer']=='OV'),:]
# genes = list(set([item.split('|')[1] for item in final_n['source']]))
# print(len(genes))
# sys.exit('stop')


# # look for cmv
# common_cmv = ['DLLSALQQL',
#               'TLLVYLFSL',
#               'VLEETSVML',
#               'IARLAKIPL',
#               'LLDGVTVSL',
#               'LPVESLPLL',
#               'YTSRGALYLY']
# final_cmv = final.loc[(final['strain']=='CMV') & (final['pep'].isin(common_cmv)),:]
# final_cmv.to_csv('final_cmv.txt',sep='\t')
# cmv_cancers = ['OV','BRCA','LIHC','NBL','GBM','CESC','COAD']

# store_data = []
# for pep in common_cmv:
#     final_p = final_cmv.loc[final_cmv['pep']==pep,:]
#     all_occur = final_p['cancer'].values.tolist()
#     tmp = []
#     for c in cmv_cancers:
#         if c in all_occur:
#             sub = final_p.loc[final_p['cancer']==c,:]
#             all_intensity = []
#             for item in sub['detailed_intensity']:
#                 all_intensity.extend(literal_eval(item))
#             med_intensity = np.median(all_intensity)
#             tmp.append(med_intensity)
#         else:
#             tmp.append(0)
#     store_data.append(tmp)

# cmv_df = pd.DataFrame(data=store_data,index=common_cmv,columns=cmv_cancers)
# cmv_df.to_csv('cmv_df.txt',sep='\t')




data = []
all_strains = []
for s,sub_df in final.groupby(by='strain'):
    all_strains.append(s)
    vc = sub_df['cancer'].value_counts().to_dict()
    row = []
    for c in cancers:
        row.append(vc.get(c,0))
    data.append(row)
df = pd.DataFrame.from_records(data=data,index=all_strains,columns=cancers)
df.to_csv('pathogen_strain_df.txt',sep='\t')
        


