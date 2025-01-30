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
VARIANT_ENSEMBL_GTF = '/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/NeoVerse/data/Homo_sapiens.GRCh38.110.gtf'

all_tissues = ['Adrenal gland', 'Aorta', 'Bladder', 'Bone marrow', 'Brain', 'Cerebellum', 'Colon', 'Esophagus', 'Gallbladder', 'Heart', 'Kidney', 'Liver', 
                'Lung', 'Lymph node', 'Mamma', 'Muscle', 'Myelon', 'Ovary', 'Pancreas', 'Prostate', 'Skin', 'Small intestine', 'Spleen', 'Stomach', 'Testis', 
                'Thymus', 'Thyroid', 'Tongue', 'Trachea', 'Uterus']

def get_enst2gs():
    gtf = pd.read_csv(VARIANT_ENSEMBL_GTF,sep='\t',skiprows=5,header=None)
    gtf.columns = ['chrom','source','feature','start','end','score','strand','phase','attribute']
    gtf_gene = gtf.loc[gtf['feature']=='transcript',:]
    pat1 = re.compile(r'transcript_id "(ENST\d+)";')
    pat2 = re.compile(r'gene_name "(.+?)";')
    dic1 = {}
    for row in gtf_gene.itertuples():  
        enst = re.search(pat1,row.attribute)
        if enst is not None:
            enst = enst.group(1)
            try:
                symbol = re.search(pat2,row.attribute).group(1)
            except:   # lncRNA, no gene_name
                symbol = 'unknown'
            dic1[enst] = symbol
    return dic1


# # nuorf peptide assemble
# mapping = {
#     'BRCA':'breast cancer',
#     'KIRC':'clear cell renal cell carcinoma',
#     'COAD':'colon cancer',
#     'STAD':'gastric cancer',
#     'MESO':'mesothelioma',
#     'LIHC':'liver cancer',
#     'ESCA':'esophageal cancer',
#     'CESC':'cervical cancer',
#     'BLCA':'bladder cancer',	
#     'RT':'rhabdoid tumor',
#     'AML':'acute myeloid leukemia',
#     'DLBC':'diffuse large B cell lymphoma',	
#     'GBM':'glioblastoma',	
#     'NBL':'neuroblastoma',
#     'PAAD':'pancreatic cancer',
#     'HNSC':'head and neck cancer',
#     'OV':'ovarian cancer',
#     'LUSC':'lung squamous cell carcinoma',
#     'LUAD':'lung adenocarcinoma',
#     'CHOL':'bile duct cancer',	
#     'SKCM':'melanoma'
# }

# pd.Series(mapping,name='full_name').to_csv('abbreviation.txt',sep='\t');sys.exit('stop')

# df = pd.read_csv('peptide_view_nuorf.txt',sep='\t',header=None).iloc[2:,:24]
# col = []
# for item1,item2 in zip(df[2],df[1]):
#     col.append(','.join([item1,item2.replace(' ','_')]))
# df['uid'] = col
# col = []
# for i in np.arange(df.shape[0]):
#     s = df.iloc[i]
#     s = s.iloc[3:-1]
#     s = s.astype(float)
#     s = s.where(s>0)
#     s = s.loc[s.notna()]
#     need = ','.join([mapping[cancers[item-3]] for item in s.index])
#     col.append(need)
# df['need'] = col
# df.to_csv('mannual_input_nuorf_in_patent.txt',sep='\t')



# figure out nbl and ribo
# final = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas/NBL/antigen/fdr/final_enhanced.txt',sep='\t')
# cond = [False if '[]' in item else True for item in final['presented_by_each_sample_hla']]
# final = final.loc[cond,:]
# final = final.loc[final['typ']=='nuORF',:]
# final = final.loc[final['unique'],:]
# cell_lines = ['NB_1691','NB_1771','NB-Ebc1','COG-N-415x','COG-N-440x','COG-N-471x','NB-SD','SK-N-AS']
# cond = []
# for item in final['samples']:
#     lis = item.split(',')
#     if len(set(cell_lines).intersection(set(lis))) > 0:
#         cond.append(True)
#     else:
#         cond.append(False)
# final = final.loc[cond,:]
# final.to_csv('nbl_cell_line_nuorf.txt',sep='\t',index=None)



data = []
for c in cancers:
    final_path = os.path.join(root_atlas_dir,c,'antigen','fdr','final_enhanced.txt')
    final = pd.read_csv(final_path,sep='\t')
    cond = [False if ('[]' in item) and ('(\'HLA-' not in item) else True for item in final['presented_by_each_sample_hla']]
    final = final.loc[cond,:]
    final = final.loc[final['typ']=='nuORF',:]
    data.append(final)
final = pd.concat(data,axis=0,keys=cancers).reset_index(level=-2).rename(columns={'level_0':'cancer'})
final.to_csv('all_nuorf.txt',sep='\t',index=None)
final['pep'].value_counts().to_csv('all_nuorf_peptide.txt',sep='\t')

vc = final['pep'].value_counts()
candidates = vc.loc[vc>8].index.tolist()
added = ['STIRVLSGY','RYLPSSVFL','LYLETRSEF','VLMEVTLEGK','AYPASLQTL','KSLAAELLVLK']
candidates = candidates + added

# enst2gs = get_enst2gs()
# dic = {i1:i2 for i1,i2 in zip(final['pep'],final['source'])}
# with open('freq_nuorf_pep2anno.txt','w') as f:
#     f.write('peptide\tsource\tannotation\n')
#     for k,v in dic.items():
#         if k in candidates and ';' not in v:
#             f.write('{}\t{}\t{}\n'.format(k,v,enst2gs.get(v.split('.')[0],'unknown')))


freq_nuorf_pep2anno = pd.read_csv('freq_nuorf_pep2anno.txt',sep='\t',index_col=0)['annotation'].to_dict()
candidates = list(freq_nuorf_pep2anno.keys())

pep2type = {p:sub_df['nuorf_type'].iloc[0] for p,sub_df in final.groupby(by='pep')}
safety_screen_df = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/code/hla_ligand_atlas_now_0.01.txt',sep='\t')

store_data = []
store_type = []
for pep in candidates:
    final_p = final.loc[final['pep']==pep,:]
    all_occur = final_p['cancer'].values.tolist()
    tmp = []
    for c in cancers:
        if c in all_occur:
            sub = final_p.loc[final_p['cancer']==c,:]
            all_intensity = []
            for item in sub['detailed_intensity']:
                all_intensity.extend(literal_eval(item))
            med_intensity = np.median(all_intensity)
            tmp.append(med_intensity)
        else:
            tmp.append(0)
    all_tissues = np.array(all_tissues)
    tmp_normal = np.full(shape=len(all_tissues),fill_value=0.0)
    tmp_normal_df = safety_screen_df.loc[safety_screen_df['peptide']==pep,:]
    for t,sub_df in tmp_normal_df.groupby(by='tissue'):
        med_intensity = np.median(sub_df['percentile'].values)
        indices = np.where(all_tissues == t)[0]
        tmp_normal[indices[0]] = med_intensity
    tmp = tmp + tmp_normal.tolist()
    store_data.append(tmp)
    store_type.append(pep2type[pep])


df = pd.DataFrame(data=store_data,index=candidates,columns=cancers+list(all_tissues))

ori_array = [tuple(['cancer']*21+['normal']*30),tuple(df.columns.tolist())]
mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
df.columns = mi

ori_array = [tuple(df.index.tolist()),
             tuple(store_type),
             tuple([freq_nuorf_pep2anno[item] for item in df.index])]
mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
df.index = mi

df.to_csv('peptide_view_nuorf.txt',sep='\t')


