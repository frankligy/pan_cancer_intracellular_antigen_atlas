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
import re
import bisect

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

def process_cancer_splicing(splicing_path):
    splicing = pd.read_csv(splicing_path,sep='\t',index_col=0)
    col1 = []
    col2 = []
    col3 = []
    for item in tqdm(splicing.index):
        chrom,coord = item.split(':')
        start,end = coord.split('-')
        start,end = int(start),int(end)
        col1.append(chrom)
        col2.append(start)
        col3.append(end)
    splicing['chrom'] = col1
    splicing['start'] = col2
    splicing['end'] = col3
    
    splicing_dic = {}
    for chrom,sub_df in tqdm(splicing.groupby(by='chrom')):
        a = sub_df['start'].values.tolist()
        b = sub_df['end'].values.tolist()
        c = a + b
        c = sorted(c)
        splicing_dic[chrom] = c

    return splicing_dic

def process_cancer_intron(intron_path):
    ir = pd.read_csv(intron_path,sep='\t',index_col=0)
    col1 = []
    col2 = []
    col3 = []
    col4 = []
    for item in ir.index:
        _,_,chrom,strand,_,start,end = item.split(',')
        start,end = int(start),int(end)
        col1.append(chrom)
        col2.append(strand)
        col3.append(start)
        col4.append(end)
    ir['chrom'] = col1
    ir['strand'] = col2
    ir['start'] = col3
    ir['end'] = col4

    ir_dic = {}
    for chrom,sub_df in ir.groupby(by='chrom'):
        ir_dic[chrom] = {}
        for strand,sub_df2 in sub_df.groupby(by='strand'):
            tmp = []
            for start,end in zip(sub_df2['start'],sub_df2['end']):
                tmp.append((start,end))
            tmp = list(set(tmp))
            tmp_now = []
            for item in tmp:
                tmp_now.extend(list(item))
            tmp_now = sorted(tmp_now)
            ir_dic[chrom][strand] = tmp_now
    
    return ir_dic


def has_splicing_site(coord,splicing_dic):
    flag = True
    chrom,coord,strand = coord.split(':')
    start,end = coord.split('-')
    start,end = int(start),int(end)
    all_ss = splicing_dic[chrom]
    pos_s = bisect.bisect_left(all_ss,start)
    pos_e = bisect.bisect_left(all_ss,end)
    if pos_s != pos_e:
        pos_s = bisect.bisect_right(all_ss,start)
        pos_e = bisect.bisect_right(all_ss,end)
        if pos_s != pos_e:
            flag = False
    return flag

def has_ir(coord,ir_dic):
    flag = True
    chrom,coord,strand = coord.split(':')
    start,end = coord.split('-')
    start,end = int(start),int(end)
    # plus
    all_bound = ir_dic[chrom]['+'] 
    pos_s = bisect.bisect_left(all_bound,start)
    pos_e = bisect.bisect_left(all_bound,end)
    try:
        assert pos_s == pos_e
        assert pos_s % 2 == 0
    except AssertionError:
        pos_s = bisect.bisect_right(all_bound,start)
        pos_e = bisect.bisect_right(all_bound,end)
        try:
            assert pos_s == pos_e
            assert pos_s % 2 == 0
        except AssertionError:
            flag = False

    if flag:
        # minus
        all_bound = ir_dic[chrom]['-'] 
        pos_s = bisect.bisect_left(all_bound,start)
        pos_e = bisect.bisect_left(all_bound,end)
        try:
            assert pos_s == pos_e
            assert pos_s % 2 == 0
        except AssertionError:
            pos_s = bisect.bisect_right(all_bound,start)
            pos_e = bisect.bisect_right(all_bound,end)
            try:
                assert pos_s == pos_e
                assert pos_s % 2 == 0
            except AssertionError:
                flag = False

    return flag


df = pd.read_csv('ts_te_antigen.txt',sep='\t')

safety_screen_df = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/code/post_safety_screen.txt',sep='\t')
safety_screen_df = safety_screen_df.loc[~safety_screen_df['cond_stringent'],:]
safety_screen_bl = list(set(safety_screen_df['pep'].values.tolist()))

df = df.loc[~df['pep'].isin(safety_screen_bl),:]


for cancer in cancers:
    root_atlas_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas'
    # splicing_dic = process_cancer_splicing(os.path.join(root_atlas_dir,cancer,'splicing_all.txt'))
    # ir_dic = process_cancer_intron(os.path.join(root_atlas_dir,cancer,'intron_all.txt'))

    # with open('splicing_ir_dic/{}_splicing_dic.p'.format(cancer),'wb') as f:
    #     pickle.dump(splicing_dic,f)

    # with open('splicing_ir_dic/{}_ir_dic.p'.format(cancer),'wb') as f:
    #     pickle.dump(ir_dic,f)


    with open('splicing_ir_dic/{}_splicing_dic.p'.format(cancer),'rb') as f:
        splicing_dic = pickle.load(f)

    with open('splicing_ir_dic/{}_ir_dic.p'.format(cancer),'rb') as f:
        ir_dic = pickle.load(f)

    df_now = df.loc[df['cancer']==cancer,:]
    col1 = []
    col2 = []
    for row in df_now.itertuples():
        cancer = row.cancer
        coord = row.source.split('|')[1]
        flag = has_splicing_site(coord,splicing_dic)
        col1.append(flag)
        flag = has_ir(coord,ir_dic)
        col2.append(flag)
    df_now['not_has_ss'] = col1
    df_now['not_in_ir'] = col2
    df_now.to_csv('splicing_ir_dic/autonomy_ts_te_antigens_{}.txt'.format(cancer),sep='\t')


# assemble
df_list = []
for c in cancers:
    df = pd.read_csv('splicing_ir_dic/autonomy_ts_te_antigens_{}.txt'.format(c),sep='\t')
    df_list.append(df)
final = pd.concat(df_list,axis=0,keys=cancers).reset_index(level=-2).rename(columns={'level_0':'cancer'}).drop(columns='Unnamed: 0')
final.to_csv('splicing_ir_dic/final.txt',sep='\t',index=None)


# draw
final = pd.read_csv('splicing_ir_dic/final.txt',sep='\t')
te_all = pd.read_csv('te_all_antigens.txt',sep='\t')
te_all = te_all.loc[~te_all['pep'].isin(safety_screen_bl),:]
te_all = te_all.loc[te_all['typ']=='TE_chimeric_transcript',:]
te_all = te_all.loc[te_all['unique'],:]
te_all['cid'] = [item.split('TE_info:')[1].split(',')[3] for item in te_all['source']]

col1 = []
col2 = []
col3 = []
for cid, sub_df in final.groupby(by='cid'):
    n_auto = sub_df.loc[(sub_df['not_has_ss']) & (sub_df['not_in_ir']),:].shape[0]
    n_chi = sub_df.shape[0] - n_auto + te_all.loc[te_all['cid']==cid,:].shape[0]
    col1.append(n_auto)
    col2.append(n_chi)
    col3.append(cid)
df = pd.DataFrame(index=col3,data={'autonomous':col1,'chimera':col2})
df.plot.bar(rot=90)
plt.savefig('te_auto_chi.pdf',bbox_inches='tight')
plt.close()

df.T.plot.pie(subplots=True,figsize=(20,6))
plt.savefig('te_auto_chi_pie.pdf',bbox_inches='tight')
plt.close()

df.to_csv('te_auto_chi.txt',sep='\t')
    




