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
from io import StringIO

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


# nuorf peptide assemble
mapping = {
    'BRCA':'breast cancer',
    'KIRC':'clear cell renal cell carcinoma',
    'COAD':'colon cancer',
    'STAD':'gastric cancer',
    'MESO':'mesothelioma',
    'LIHC':'liver cancer',
    'ESCA':'esophageal cancer',
    'CESC':'cervical cancer',
    'BLCA':'bladder cancer',	
    'RT':'rhabdoid tumor',
    'AML':'acute myeloid leukemia',
    'DLBC':'diffuse large B cell lymphoma',	
    'GBM':'glioblastoma',	
    'NBL':'neuroblastoma',
    'PAAD':'pancreatic cancer',
    'HNSC':'head and neck cancer',
    'OV':'ovarian cancer',
    'LUSC':'lung squamous cell carcinoma',
    'LUAD':'lung adenocarcinoma',
    'CHOL':'bile duct cancer',	
    'SKCM':'melanoma'
}

pd.Series(mapping,name='full_name').to_csv('abbreviation.txt',sep='\t')

# # ribo-seq nbl and immunopeptidome (no safety screen)
# ribo_map = {
#     'OHMX20210030_001':'COG-N-415x',
#     'OHMX20210030_002':'NB-Ebc1',
#     'OHMX20210030_003':'NB_1691',
#     'OHMX20210030_004':'COG-N-471x',
#     'OHMX20210030_005':'COG-N-440x',
#     'OHMX20210030_006':'NB_1771',
#     'OHMX20210030_007':'NB-SD',
#     'OHMX20210030_008':'SK-N-AS'
# }
# reverse_ribo_map = {v:k for k,v in ribo_map.items()}
# ribo_root_dir = '/gpfs/data/yarmarkovichlab/neuroblastoma/riboseq/ribo_real'

# final_path = os.path.join(root_atlas_dir,'NBL','antigen','0.01','final_enhanced.txt')
# final = pd.read_csv(final_path,sep='\t')
# cond = [False if ('[]' in item) and ('(\'HLA-' not in item) else True for item in final['presented_by_each_sample_hla']]
# final = final.loc[cond,:]
# final = final.loc[final['typ']=='nuORF',:]
# pdx = ['SK-N-AS','NB_1691','NB_1771','NB-Ebc1','COG-N-415x','COG-N-440x','COG-N-471x','NB-SD']
# df_list = []
# for p in pdx:
#     print(p)
#     col = []
#     sub = final.loc[final['samples'].str.contains(p),:]
#     table_path = os.path.join(ribo_root_dir,reverse_ribo_map[p],'outputDir','riborf_table.txt')
#     for pep in sub['pep']:
#         cmd = 'grep "{}" {}'.format(pep,table_path)
#         try:
#             df = pd.read_csv(StringIO(subprocess.run(cmd,shell=True,capture_output=True).stdout.decode('utf-8')),sep='\t',header=None)
#         except:
#             col.append(0)
#         else:
#             max_prob = df.iloc[:,3].values.max()
#             col.append(max_prob)
#     sub['max_prob'] = col
#     df_list.append(sub)
# result = pd.concat(df_list,axis=0,keys=pdx).reset_index(level=-2).rename(columns={'level_0':'pdx'})
# count = 0
# total_pep = 0
# prob_list = []
# for pep,sub_df in result.groupby(by='pep'):
#     max_prob_pep = sub_df.sort_values(by='max_prob')['max_prob'].iloc[-1]
#     total_pep += 1
#     prob_list.append(max_prob_pep)
#     if max_prob_pep > 0:
#         count += 1
# prop = (count + 1) / total_pep  # rescue LYLETRSEF
# print(prop)
# result.to_csv('nbl_nuorf_pdx.txt',sep='\t',index=None)
# sns.ecdfplot(data=prob_list)
# plt.savefig('nbl_nuorf_pdx.pdf',bbox_inches='tight')
# plt.close()

df = pd.read_csv('./stats/final_all_ts_antigens.txt',sep='\t')
final = df.loc[df['typ']=='nuORF',:]
final.to_csv('all_nuorf_final.txt',sep='\t',index=None)

ts = ['5\' uORF','5\' Overlap uORF','Out-of-Frame','3\' Overlap dORF','3\' dORF','lncRNA','Pseudogene']
dic = {t:[] for t in ts}
for pep,sub_df in final.groupby(by='pep'):
    nc = sub_df.shape[0]
    nt = sub_df['nuorf_type'].iloc[0]
    if nt in ts:
        dic[nt].append((pep,nc))
new_dic = {}
for k,v in dic.items():
    tmp = sorted(v,key=lambda x:x[1],reverse=True)
    p_sorted, _ = zip(*tmp)
    new_dic[k] = p_sorted[:3]

candidates = []
for k,v in new_dic.items():
    candidates.extend(v)

enst2gs = get_enst2gs()
dic = {i1:i2 for i1,i2 in zip(final['pep'],final['source'])}
with open('freq_nuorf_pep2anno.txt','w') as f:
    f.write('peptide\tsource\tannotation\n')
    for k,v in dic.items():
        if k in candidates and ';' not in v:
            f.write('{}\t{}\t{}\n'.format(k,v,enst2gs.get(v.split('.')[0],'unknown')))

freq_nuorf_pep2anno = pd.read_csv('freq_nuorf_pep2anno.txt',sep='\t',index_col=0)['annotation'].to_dict()
candidates = list(freq_nuorf_pep2anno.keys())

pep2type = {p:sub_df['nuorf_type'].iloc[0] for p,sub_df in final.groupby(by='pep')}
safety_screen_df = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/code/hla_ligand_atlas_now_0.05_tesorai.txt',sep='\t')

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
    tmp_ribo = np.full(shape=8,fill_value=0.0) # healthy ribo
    tmp_normal = np.full(shape=len(all_tissues),fill_value=0.0)
    tmp_normal_df = safety_screen_df.loc[safety_screen_df['peptide']==pep,:]
    for t,sub_df in tmp_normal_df.groupby(by='tissue'):
        med_intensity = np.median(sub_df['percentile'].values)
        indices = np.where(all_tissues == t)[0]
        tmp_normal[indices[0]] = med_intensity
    tmp = tmp + tmp_ribo.tolist() + tmp_normal.tolist()
    store_data.append(tmp)
    store_type.append(pep2type[pep])

all_ribo_tissues = ['Brain','fat','Fibroblast','HAEC','HCAEC','Hepatocytes','VSMC','Kidney']
df = pd.DataFrame(data=store_data,index=candidates,columns=cancers+all_ribo_tissues+list(all_tissues))

ori_array = [tuple(['cancer']*21+['ribo']*8+['normal']*30),tuple(df.columns.tolist())]
mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
df.columns = mi

ori_array = [tuple(df.index.tolist()),
             tuple(store_type),
             tuple([freq_nuorf_pep2anno[item] for item in df.index])]
mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
df.index = mi

df.to_csv('peptide_view_nuorf.txt',sep='\t')


