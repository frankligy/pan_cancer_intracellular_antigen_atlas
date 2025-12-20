#!/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/bin/python3.7

import pandas as pd
import numpy as np
import sys,os
from Bio.SeqIO.FastaIO import SimpleFastaParser
import subprocess
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import rankdata

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'



# screen

def rewrite_msmsScans_new(valids,fold,fdr):

    name = 'msmsScans'
    msmsScans_path = os.path.join(fold,'combined','txt','{}.txt'.format(name))
    final = pd.read_csv(msmsScans_path,sep='\t')
    final['uid'] = [','.join([str(item1),str(item2)]) for item1,item2 in zip(final['Raw file'],final['Scan number'])]
    final['uid_seq'] = [','.join([str(item1),str(item2),str(item3)]) for item1,item2,item3 in zip(final['Raw file'],final['Scan number'],final['Sequence'])]

    cols = []
    for anno,valid in valids.items():
        col = final['uid_seq'].isin(valid).values
        final['Identified_{}'.format(anno)] = col
        cols.append(col)
    
    final_col = ['+' if cond else None for cond in np.any(cols,axis=0)]
    final['Identified'] = final_col

    final.to_csv(os.path.join(fold,'combined','txt','{}_new_{}.txt'.format(name,fdr)),sep='\t',index=None)

def rederive_fdr_vanilla(fold,fdr):
    msms_path = os.path.join(fold,'combined','txt','msms.txt')
    df = pd.read_csv(msms_path,sep='\t')
    df = df.sort_values(by='PEP',ascending=True)
    cond = [1 if isinstance(item,str) else 0 for item in df['Reverse']]
    total_decoy = np.cumsum(cond)
    prop_decoy = total_decoy / (np.arange(len(total_decoy)) + 1)
    col1 = np.where(prop_decoy<fdr,'+',None)
    df['Identified'] = col1

    # now get valid
    df = df.loc[df['Identified']=='+',:]
    valid = []
    for item1,item2,item3 in zip(df['Raw file'],df['Scan number'],df['Sequence']):
        valid.append(','.join([str(item1),str(item2),str(item3)]))  

    return valid


result_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/immuno'
old_dir = os.getcwd()
os.chdir(result_dir)
all_tissues = subprocess.run("for f in *; do echo $f; done",shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
os.chdir(old_dir)

# # derive maxquant 0.05 
# fdr = 0.05
# # rederive
# for t in all_tissues:
#     t_dir = os.path.join(result_dir,t)
#     old_dir = os.getcwd()
#     os.chdir(t_dir)
#     all_batches = subprocess.run("for f in *; do echo $f; done",shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#     os.chdir(old_dir)
#     for b in all_batches:
#         print(fdr,t,b)
#         fold = os.path.join(t_dir,b)
#         valid = rederive_fdr_vanilla(fold,fdr)
#         valids = {'vanilla':valid}
#         rewrite_msmsScans_new(valids,fold,fdr)

# # append tesorai result to each 5% FDR msmsScan
# fdr = 0.05
# tesorai = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/NYU_Tesorai_all_searches/tesorai_peptide_fdr_normal.tsv',sep='\t')
# tesorai = tesorai.loc[tesorai['qval']<0.01,:]
# tesorai['uid'] = [','.join([item1.split('.')[0],str(item2)]) for item1,item2 in zip(tesorai['filename'],tesorai['scan_id'])]
# uid2seq = pd.Series(index=tesorai['uid'].values,data=tesorai['clean_sequence'].values).to_dict()
# uid2protein = pd.Series(index=tesorai['uid'].values,data=tesorai['possible_protein_ids'].values).to_dict()
# uid2score = pd.Series(index=tesorai['uid'].values,data=tesorai['score'].values).to_dict()
# for t in all_tissues:
#     t_dir = os.path.join(result_dir,t)
#     old_dir = os.getcwd()
#     os.chdir(t_dir)
#     all_batches = subprocess.run("for f in *; do echo $f; done",shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#     os.chdir(old_dir)
#     for b in all_batches:
#         print(fdr,t,b)
#         fold = os.path.join(t_dir,b)
#         name = 'msmsScans'
#         final = pd.read_csv(os.path.join(fold,'combined','txt','{}_new_{}.txt'.format(name,fdr)),sep='\t')
#         # remove rev
#         col1 = []
#         col2 = []
#         for i1,i2,i3 in zip(final['Identified'],final['Identified_vanilla'],final['Proteins']):
#             if 'REV' not in i3:
#                 col1.append(i1)
#                 col2.append(i2)
#             else:
#                 col1.append(None)
#                 col2.append(False)
#         final['Identified'] = col1
#         final['Identified_vanilla'] = col2
#         # add tesorai
#         final['tesorai_sequence'] = final['uid'].map(uid2seq).values
#         final['tesorai_proteins'] = final['uid'].map(uid2protein).values
#         final['tesorai_score'] = final['uid'].map(uid2score).values
#         # scrutinize
#         col1 = []
#         col2 = []
#         col3 = []
#         col4 = []
#         col5 = []
#         for i1,i2,i3,i4,i5,i6,i7 in zip(final['Identified'],final['Sequence'],final['Proteins'],final['Score'],final['tesorai_sequence'],final['tesorai_proteins'],final['tesorai_score']):
#             if isinstance(i6,float):  # tesorai empty protein
#                 col1.append(i1)
#                 col2.append(i2)
#                 col3.append(i3)
#                 col4.append(False)
#                 continue
#             if i2 == i5:     # same sequence
#                 if i1 == '+':
#                     col1.append(i1)
#                     col2.append(i2)
#                     col3.append(i3)
#                     col4.append(False)
#                     continue  
#                 else:
#                     col1.append('+')
#                     col2.append(i5)
#                     col3.append(';'.join(i6.split(';;')))
#                     col4.append(True)
#                     continue
#             if i1 == '+':
#                 if (i4 >= 70 and i7 >= 5) or (i4 < 70 and i7 < 5):
#                     col1.append(None)
#                     col2.append(i2)
#                     col3.append(i3)
#                     col4.append(False)
#                 elif i4 >= 70 and i7 < 5:
#                     col1.append(i1)
#                     col2.append(i2)
#                     col3.append(i3)
#                     col4.append(False)
#                 elif i4 < 70 and i7 >= 5:
#                     col1.append('+')
#                     col2.append(i5)
#                     col3.append(';'.join(i6.split(';;')))
#                     col4.append(False)
#             else:
#                 if isinstance(i5,str):
#                     col1.append('+')
#                     col2.append(i5)
#                     col3.append(';'.join(i6.split(';;')))
#                     col4.append(True)
#                 else:
#                     col1.append(i1)
#                     col2.append(i2)
#                     col3.append(i3)
#                     col4.append(False)
#         final['Identified'] = col1
#         final['Sequence'] = col2
#         final['Proteins'] = col3
#         final['Tesorai'] = col4
#         final['Length'] = [len(item) if len(item)>1 else 0 for item in final['Sequence']]
#         final.to_csv(os.path.join(fold,'combined','txt','{}_new_{}_tesorai.txt'.format(name,fdr)),sep='\t',index=None)


# # combining and renaming to resolve small case issue
# fdr = 0.05
# total_dfs = []
# for t in all_tissues:
#     each_tissue_dfs = []
#     each_tissue_raws = []
#     intdir = os.path.join(result_dir,t)
#     old_dir = os.getcwd()
#     os.chdir(intdir)
#     all_batches = subprocess.run("for f in batch*; do echo $f; done",shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#     os.chdir(old_dir)
#     for b in all_batches:
#         msms_path = os.path.join(intdir,b,'combined','txt','msmsScans_new_{}_tesorai.txt'.format(fdr))
#         msms = pd.read_csv(msms_path,sep='\t')
#         msms = msms.loc[msms['Identified']=='+',:]
#         for raw,sub_df in msms.groupby(by='Raw file'):
#             sub_df['Precursor intensity'] = sub_df['Precursor intensity'].fillna(value=1e-5)
#             sub_df['percentile'] = rankdata(sub_df['Precursor intensity'].values,method='min') / sub_df.shape[0]
#             each_raw_data = []
#             for p,sub_df2 in sub_df.groupby(by='Sequence'):
#                 # take the highest
#                 intensity = sub_df2['Precursor intensity'].values.max()
#                 percentile = sub_df2['percentile'].values.max()
#                 each_raw_data.append((p,intensity,percentile))
#             each_raw_df = pd.DataFrame.from_records(data=each_raw_data,columns=['peptide','intensity','percentile'])
#             try:
#                 upper = np.quantile(each_raw_df['intensity'].values,0.75)
#             except:
#                 continue
#             else:
#                 each_raw_df['norm'] = np.log2(each_raw_df['intensity'].values/upper)
#                 each_tissue_dfs.append(each_raw_df)
#                 each_tissue_raws.append(raw)
#     each_tissue_metadf = pd.concat(each_tissue_dfs,axis=0,keys=each_tissue_raws).reset_index(level=-2).rename(columns={'level_0':'raw_file'})
#     total_dfs.append(each_tissue_metadf)

# final = pd.concat(total_dfs,axis=0,keys=all_tissues).reset_index(level=-2).rename(columns={'level_0':'tissue'})
# final['log_intensity'] = np.log2(final['intensity'].values)
# mapping =  {'bladder':'Bladder','brain':'Brain','heart':'Heart','skin':'Skin','smallintestine':'Small intestine','spleen':'Spleen','tongue':'Tongue',
#             'AdrenalGland':'Adrenal gland','BoneMarrow':'Bone marrow','LymphNode':'Lymph node'}
# col = []
# for item in final['tissue']:
#     col.append(mapping.get(item,item))
# final['tissue'] = col
# final.to_csv('hla_ligand_atlas_now_{}_tesorai.txt'.format(fdr),sep='\t',index=None)


# just use this tesorai combined one and original hla ligand atlas
df = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/codes/summary/for_safety_screen.txt',sep='\t')
final = df

fdr = 0.05
df = pd.read_csv('hla_ligand_atlas_now_{}_tesorai.txt'.format(fdr),sep='\t')
dic = {}
for pep,sub_df in df.groupby(by='peptide'):
    ts = ','.join(list(set(sub_df['tissue'].values.tolist())))
    dic[pep] = ts
col = []
for item in final['pep']:
    col.append(dic.get(item,''))
final['normal_{}'.format(fdr)] = col

all_tissues = ['Adrenal gland', 'Aorta', 'Bladder', 'Bone marrow', 'Brain', 'Cerebellum', 'Colon', 'Esophagus', 'Gallbladder', 'Heart', 'Kidney', 'Liver', 
                'Lung', 'Lymph node', 'Mamma', 'Muscle', 'Myelon', 'Ovary', 'Pancreas', 'Prostate', 'Skin', 'Small intestine', 'Spleen', 'Stomach', 'Testis', 
                'Thymus', 'Thyroid', 'Tongue', 'Trachea', 'Uterus']

non_essential = ['Adrenal gland','Ovary','Prostate','Testis','Thymus']
essential = list(set(all_tissues).difference(set(non_essential)))

cond = []
for i1,i2 in zip(final['hla_ligand_atlas'],final['normal_0.05']):
    lis = []
    for i in [i1,i2]:
        if isinstance(i,str):
            lis.extend(i.split(','))
    common = set(lis).intersection(set(essential))
    if len(common) > 0:
        cond.append(False)
    else:
        cond.append(True)
final['cond_stringent'] = cond
final.to_csv('post_safety_screen.txt',sep='\t',index=None)
sys.exit('stop')


# making db
df = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/codes/summary/for_safety_screen.txt',sep='\t')
with open('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/db_fasta/peptides.fasta','w') as f:
    for cancer,pep,typ in zip(df['cancer'],df['pep'],df['typ']):
        f.write('>query|{}|{}|{}\n{}\n'.format(pep,typ,cancer,pep))
sys.exit('stop')

# transfer
# aman_dir = '/gpfs/data/yarmarkovichlab/Aman/hla_ligand_atlas'
# frank_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen'
# df = pd.read_csv(os.path.join(aman_dir,'final.txt'),sep='\t')
# all_tissues = df['tissue'].unique()
# all_tissues = ['brain','bladder','heart','skin','smallintestine','spleen','tongue']  # case issue
# for tissue in all_tissues:
#     print(tissue)
#     d = os.path.join(aman_dir,tissue)
#     all_raws = subprocess.run('find {} -type f -name "*.raw"'.format(d),shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#     # create batch folder if needed
#     batches = set()
#     for raw in all_raws:
#         if len(raw.split('/')) == 9:
#             b = raw.split('/')[-2]
#             batches.add(b)
#     for b in batches:
#         if not os.path.exists(os.path.join(frank_dir,tissue,b)):
#             os.makedirs(os.path.join(frank_dir,tissue,b))
#     if len(batches) == 0:
#         if not os.path.exists(os.path.join(frank_dir,tissue)):
#             os.makedirs(os.path.join(frank_dir,tissue))
#     # copy raw files
#     if len(batches) > 0:
#         for raw in all_raws:
#             b = raw.split('/')[-2]
#             f = raw.split('/')[-1]
#             subprocess.run('cp {} {}'.format(raw,os.path.join(frank_dir,tissue,b)),shell=True)
#     else:
#         for raw in all_raws:
#             f = raw.split('/')[-1]
#             subprocess.run('cp {} {}'.format(raw,os.path.join(frank_dir,tissue)),shell=True)

# # I then remove the case incorrect ones 
# # I gave each tissue a batch to be compatible


