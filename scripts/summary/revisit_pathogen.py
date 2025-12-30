#!/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/bin/python3.7

import pandas as pd
import numpy as np
import sys,os
from Bio.SeqIO.FastaIO import SimpleFastaParser
import subprocess
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from tqdm import tqdm
from Bio.Seq import Seq

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

added = {
    'STAD': ['Helicobacter pylori','Campylobacter ureolyticus','Niallia circulans','Clostridium intestinale'],
    'COAD': ['Fusobacterium nucleatum'],
    'SKCM': ['Fusobacterium nucleatum'],
    'HNSC': ['Fusobacterium nucleatum'],
    'ESCA': ['Fusobacterium nucleatum','Campylobacter ureolyticus','Niallia circulans','Clostridium intestinale'],
    'OV':['Campylobacter ureolyticus','Niallia circulans','Clostridium intestinale'],
    'CESC':['Campylobacter ureolyticus','Niallia circulans','Clostridium intestinale']
}

root_atlas_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas'
T2T_NORMAL_META = '/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/NeoVerse/GTEx/selected/t2t_normal_meta.txt'
T2T_NORMAL_DIR = '/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/NeoVerse/GTEx/selected/t2t_pathogen_intron'

# make sure not in genus level in normal
meta = pd.read_csv(T2T_NORMAL_META,sep='\t',index_col=0,header=None)
meta.columns = ['srr']
tissue2srr = meta['srr'].to_dict()
normal_pathogen_list = []
all_srrs = []
all_tissues = []
for t,ss in tissue2srr.items():
    for s in ss.split(','):
        kraken2_result_path = os.path.join(T2T_NORMAL_DIR,s,'test_report.txt')
        kraken2_df = pd.read_csv(kraken2_result_path,sep='\t',header=None)
        kraken2_df = kraken2_df.loc[(kraken2_df[0].str.contains('|g__',regex=False)) & (~kraken2_df[0].str.contains('|s__',regex=False)),:]
        normal_pathogen_list.append(kraken2_df)
        all_srrs.append(s)
        all_tissues.append(t)
all_uids = [s+','+t for s,t in zip(all_srrs,all_tissues)]
normal_pathogen = pd.concat(normal_pathogen_list,axis=0,keys=all_uids).reset_index(level=-2)
normal_pathogen.columns = ['uids','strain','count']
normal_pathogen['sample'] = [item.split(',')[0] for item in normal_pathogen['uids']]
normal_pathogen['tissue'] = [item.split(',')[1] for item in normal_pathogen['uids']]

df_list = []
for c,t in zip(cancers,n_samples):
    print(c)

    # my crude way
    df = pd.read_csv(os.path.join(root_atlas_dir,c,'pathogen_rec.txt'),sep='\t',index_col=0)
    added_list = added.get(c,[])
    df = df.loc[((df['count']>t*0.2) & (df['normal'].isna())) | (df['strain'].isin(added_list)),:]
    # df = df.loc[(df['normal'].isna()) & (df['count']>t*0.2),:]

    # make sure they are abundant in TCMbio (Bracken, remove contaminants)
    dict_mean = {}
    dict_freq = {}

    if c == 'RT' or c == 'NBL':
        df_list.append(df)
        continue
    else:
        ext_bac = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/pathogen/TCMbio/{}_table_gz/{}_otutable_bacteria_species.csv'.format(c,c),sep=',',index_col=0).iloc[:,:-1]
        ext_fun = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/pathogen/TCMbio/{}_table_gz/{}_otutable_fungus_species.csv'.format(c,c),sep=',',index_col=0).iloc[:,:-1]
        ext = pd.concat([ext_bac,ext_fun],axis=1,join='outer')
        ext = ext.rename(columns=lambda x:x.replace('[','').replace(']','')).iloc[:,:-1]
        unknown_sp = []
        for sp in df['strain']:
            uid = 's__{}'.format(sp.replace(' ','_'))
            try:
                s = ext[uid]
            except:  # still need to consider if it is virus
                unknown_sp.append(sp)
            else:
                mean = s.mean()
                freq = np.count_nonzero(s.values) / len(s)
                dict_mean[sp] = mean
                dict_freq[sp] = freq
                if sp in added_list:
                    fig,ax = plt.subplots()
                    sns.histplot(s,ax=ax)
                    ax.set_title('mean:{},freq:{}'.format(round(mean,2),round(freq,2)))
                    plt.savefig(os.path.join(root_atlas_dir,c,'{}_TCMbio.pdf'.format(uid)),bbox_inches='tight')
                    plt.close()
        df['TCMbio_abundance_mean'] = df['strain'].map(dict_mean).values
        df['TCMbio_abundance_freq'] = df['strain'].map(dict_freq).values


    if df.shape[0] > 0:
        genus_list = []
        unique_genus = set()
        for item in df.index:
            genus = item.split('|s__')[0]
            if genus not in unique_genus:
                normal_genus = normal_pathogen.loc[normal_pathogen['strain']==genus,:]
                genus_list.append(normal_genus)
                unique_genus.add(genus)
        normal_genus_final = pd.concat(genus_list,axis=0)
        pathogen_dict = {}
        for strain,sub_df in normal_genus_final.groupby(by='strain'):
            pathogen_dict[strain] = {}
            for tissue,sub_df2 in sub_df.groupby(by='tissue'):
                per_tissue_strain = ','.join([str(item) for item in sub_df2['count'].values])
                pathogen_dict[strain][tissue] = per_tissue_strain
        df['genus'] = [item.split('|s__')[0] for item in df.index]
        df['normal_genus'] = df['genus'].map(pathogen_dict).values
        df.to_csv(os.path.join(root_atlas_dir,c,'pathogen_inspection.txt'),sep='\t')
    
    df_list.append(df)

result = pd.concat(df_list,axis=0,keys=cancers).reset_index(level=-2).rename(columns={'level_0':'cancer'})
result.to_csv('revist_pathogen_cut_add.txt',sep='\t',index=None)


