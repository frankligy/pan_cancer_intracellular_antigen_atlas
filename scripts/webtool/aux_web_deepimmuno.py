#!/gpfs/data/yarmarkovichlab/Frank/test_snaf/test_snaf_env/bin/python3.7 

import pandas as pd
import numpy as np
import sys,os
from tqdm import tqdm
import snaf
from snaf.deepimmuno import run_deepimmuno
from ast import literal_eval


meta_dic = {}

root_atlas_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas'
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

for c in cancers:
    final_path = os.path.join(root_atlas_dir,c,'antigen','fdr','final_enhanced.txt')
    final = pd.read_csv(final_path,sep='\t')
    cond = [False if ('[]' in item) and ('(\'HLA-' not in item) else True for item in final['presented_by_each_sample_hla']]
    final = final.loc[cond,:]
    for row in tqdm(final.itertuples()):
        pep = row.pep
        lists = literal_eval(row.additional_query)
        hlas = literal_eval(row.presented_by_each_sample_hla)
        # adding hlas to lists
        for k,vs in hlas.items():
            if len(vs) == 0:
                continue
            else:
                for v in vs:
                    if v[0] is not None:
                        lists.append(v)
        lists = list(set(lists))
        df = pd.DataFrame.from_records(lists,columns=['hla','rank_pert','nM','id'])
        all_hla = [item.replace(':','') for item in df['hla']]
        if pep in meta_dic.keys():
            meta_dic[pep] = list(set(meta_dic[pep]).union(set(all_hla)))
        else:
            meta_dic[pep] = all_hla

data = []
for k,vs in meta_dic.items():
    l = len(k)
    for v in vs:
        data.append((k,v,l))
df = pd.DataFrame.from_records(data,columns=['pep','hla','length'])
df = df.loc[df['length'].isin([9,10]),:]

df_input = pd.DataFrame(data={0:df['pep'].values.tolist(),1:df['hla'].values.tolist()})
df_output = run_deepimmuno(df_input)
df_output.to_csv('all_deepimmuno_immunogenicity.txt',sep='\t',index=None)





