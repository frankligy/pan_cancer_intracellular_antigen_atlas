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
te_gtf_path = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/splicing/GRCh38_Ensembl_rmsk_TE.gtf'
hg38_normal_dir = '/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/NeoVerse/GTEx/selected/hg38_telocal_intron'
SPLICING_GTEX_HG38='/gpfs/data/yarmarkovichlab/Frank/SNAF_ecosystem/data/controls/GTEx_junction_counts.h5ad'
SPLICING_TCGA_HG38='/gpfs/data/yarmarkovichlab/Frank/SNAF_ecosystem/data/controls/tcga_matched_control_junction_count.h5ad'
SPLICING_GTEX_MAPPING_HG38='/gpfs/data/yarmarkovichlab/chordoma/NeoVerse_analysis/splicing_annotations/gtex_hg38_t2t/tmp_prelift.bed'
SPLICING_TCGA_MAPPING_HG38='/gpfs/data/yarmarkovichlab/chordoma/NeoVerse_analysis/splicing_annotations/tcga_hg38_t2t/tmp_prelift.bed'


def process_splicing_gtex(events):
    adata_gtex = ad.read_h5ad(SPLICING_GTEX_HG38)
    adata_tcga = ad.read_h5ad(SPLICING_TCGA_HG38)
    # add suffix to tcga tissue
    adata_tcga.var['tissue'] = [item + '_paratumor' for item in adata_tcga.var['tissue']]
    # build gtex_map
    gtex_map = pd.read_csv(SPLICING_GTEX_MAPPING_HG38,sep='\t',header=None)
    gtex_map.columns = ['chrom','start_1','end','uid']
    gtex_mapping = {}
    for row in gtex_map.itertuples():
        gtex_mapping['{}:{}-{}'.format(row.chrom,str(int(row.start_1+1)),row.end)] = row.uid
    # build tcga_map
    tcga_map = pd.read_csv(SPLICING_TCGA_MAPPING_HG38,sep='\t',header=None)
    tcga_map.columns = ['chrom','start_1','end','uid']
    tcga_mapping = {}
    for row in tcga_map.itertuples():
        tcga_mapping['{}:{}-{}'.format(row.chrom,str(int(row.start_1+1)),row.end)] = row.uid
    # start to analyze to get adata_single
    uid2medians = {}
    for event in tqdm(events):
        uid = gtex_mapping.get(event,None)
        if uid is not None:
            adata_gtex_single = adata_gtex[uid,:]
        else:
            adata_gtex_single = ad.AnnData(X=csr_matrix(np.full((1,adata_gtex.shape[1]),0)),obs=pd.DataFrame(data={'mean':[0],'std':[0]},index=[uid]),var=adata_gtex.var)
            print('imputing gtex for {}'.format(event))

        uid = tcga_mapping.get(event,None)
        if uid is not None:
            adata_tcga_single = adata_tcga[uid,:]
        else:
            adata_tcga_single = ad.AnnData(X=csr_matrix(np.full((1,adata_tcga.shape[1]),0)),obs=pd.DataFrame(data={'mean':[0],'std':[0]},index=[uid]),var=adata_tcga.var)
            adata_tcga_single.obs_names = [uid]
            print('imputing tcga for {}'.format(event))

        adata_single = ad.concat([adata_gtex_single,adata_tcga_single],axis=1,join='outer',merge='first')
        # get medians for each tissue
        all_tissues = adata_single.var['tissue'].unique().tolist()
        tmp_data = []
        for t in all_tissues:
            values = adata_single[:,adata_single.var['tissue']==t].X.toarray().reshape(-1)
            tmp_data.append(np.mean(values))  # we are still say median but its' mean
        uid2medians[event] = tmp_data
    return uid2medians,all_tissues


def process_tumor_splicing():
    dic = {}
    for c in tqdm(cancers):
        splicing_path = os.path.join(root_atlas_dir,c,'splicing_rec.txt')
        splicing = pd.read_csv(splicing_path,sep='\t',index_col=0)
        dic[c] = splicing['ave_count_tumor'].to_dict()
    return dic



data = []
for c in cancers:
    final_path = os.path.join(root_atlas_dir,c,'antigen','fdr','final_enhanced.txt')
    final = pd.read_csv(final_path,sep='\t')
    cond = [False if ('[]' in item) and ('(\'HLA-' not in item) else True for item in final['presented_by_each_sample_hla']]
    final = final.loc[cond,:]
    final = final.loc[final['typ']=='splicing',:]
    data.append(final)
final = pd.concat(data,axis=0,keys=cancers).reset_index(level=-2).rename(columns={'level_0':'cancer'})
final.to_csv('all_splicing.txt',sep='\t',index=None)
final['pep'].value_counts().to_csv('all_splicing_peptide.txt',sep='\t')
sys.exit('stop')

# peptide view
# peptide = pd.read_csv('splicing_peptides.txt',sep='\t',index_col=0)
# pep2type = peptide['type'].to_dict()
# pep2anno = peptide['annotation'].to_dict()
# store_data = []
# store_type = []
# store_anno = []
# for k,v in pep2type.items():
#     final_p = final.loc[final['pep']==k,:]
#     all_occur = final_p['cancer'].values.tolist()
#     tmp = []
#     for c in cancers:
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
#     store_type.append(v)
#     store_anno.append(pep2anno[k])
    
# df = pd.DataFrame(data=store_data,index=pep2type.keys(),columns=cancers)
# ori_array = [tuple(df.index.tolist()),
#              tuple(store_type),
#              tuple(store_anno)]
# mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
# df.index = mi
# df.to_csv('peptide_view_splicing.txt',sep='\t')




# find tumor specific one
all_uids = []
uid2full = {}
for item in final['source']:
    tmp_list = item.split(';')
    for i in tmp_list:
        if i.startswith('chr'):
            uid = i.split('|')[0]
            all_uids.append(uid)
            uid2full.setdefault(uid,[]).append(i)
            break
all_uids = list(set(all_uids))
tmp = {}
for k,v in uid2full.items():
    if len(v) == 1:
        tmp[k] = v[0]
    else:
        lfcs = []
        for item in v:
            lfcs.append(float(item.split('|')[4]))
        tmp[k] = v[np.argmax(lfcs)]
uid2full = tmp


uid2medians,all_tissues = process_splicing_gtex(all_uids)
dic = process_tumor_splicing()


uid2tumors = {}
for uid in all_uids:
    data = []
    for k,v in dic.items():
        data.append(v.get(uid,0))
    uid2tumors[uid] = data

all_data = []
for uid in all_uids:
    data = []
    data.extend(uid2tumors[uid])
    data.extend(uid2medians[uid])
    all_data.append(data)

df = pd.DataFrame.from_records(all_data,columns=cancers+all_tissues,index=all_uids)

ori_array = [tuple(['cancer']*21+['normal']*51+['paratumor']*21),tuple(df.columns.tolist())]
mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
df.columns = mi

ori_array = [tuple(df.index.tolist()),
             tuple([uid2full[item].split('|')[4] for item in df.index]),
             tuple([uid2full[item].split('|')[-2] for item in df.index])]
mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
df.index = mi

df = df.loc[:,(['cancer','normal'],slice(None))]


df.to_csv('splicing_holy.txt',sep='\t')
sys.exit('stop')



