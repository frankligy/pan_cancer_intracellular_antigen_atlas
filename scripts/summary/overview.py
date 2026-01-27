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
import math
import itertools

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

cancers2immuno = {
    'BRCA':'breast_cancer',
    'KIRC':'kidney_clear_cell',
    'COAD':'colon_cancer',
    'STAD':'stomach_cancer',
    'MESO':'mesothelioma',
    'LIHC':'liver_cancer',
    'ESCA':'esophageal_cancer',
    'CESC':'cervical_cancer',
    'BLCA':'bladder_cancer',
    'RT':'rhabdoid_tumor',
    'AML':'AML',
    'DLBC':'DLBC',
    'GBM':'GBM',
    'NBL':'neuroblastoma',
    'PAAD':'pancreatic_cancer',
    'HNSC':'head_and_neck',
    'OV':'ovarian_cancer',
    'LUSC':'lung_LUSC',
    'LUAD':'lung_LUAD',
    'CHOL':'bile_duct_cancer',
    'SKCM':'melanoma'
}

with open('manual_bl.txt','r') as f:
    manual_bl = [gene.rstrip('\n') for gene in f.readlines()]

root_atlas_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas'
bayests_xy_path = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/gene/full_results_XY_essential_tissues.txt'
protein_path = '/gpfs/data/yarmarkovichlab/public/ImmunoVerse/database/ensembl_protein.fasta'


def chop_normal_pep_db(fasta_path,output_path,mers,allow_duplicates):
    '''
    chop any normal human proteome to certain mers

    :param fasta_path: string, the path to the human protein fasta file
    :param output_path: string, the path to the output fasta file
    :param mers: list, like [9,10] will generate 9mer and 10mer
    :param allow_duplicates: boolean. whether allow duplicate or not

    Example::

        chop_normal_pep_db(fasta_path='human_uniprot_proteome.fasta',output_path='./human_uniprot_proteome_mer9_10.fasta',mers=[9,10],allow_duplicates=False)
    '''
    # for human genome in uniprot, 9-10mer, remove duplicates will decrease from 44,741,578 to 41,638,172
    if allow_duplicates:
        with open(fasta_path,'r') as in_handle, open(output_path,'w') as out_handle:
            for title,seq in tqdm(SimpleFastaParser(in_handle)):
                count = 0
                length = len(seq)
                for i in range(length):
                    for mer in mers:
                        if i+mer <= length:
                            out_handle.write('>{}_{}_{}\n'.format(title,mer,count))
                            out_handle.write('{}\n'.format(seq[i:i+mer]))
                            count += 1
    else:
        with open(fasta_path,'r') as in_handle, open(output_path,'w') as out_handle:
            existing = set()
            for title,seq in tqdm(SimpleFastaParser(in_handle)):
                count = 0
                length = len(seq)
                for i in range(length):
                    for mer in mers:
                        if i+mer <= length:
                            subseq = seq[i:i+mer]
                            if subseq not in existing:                                
                                out_handle.write('>{}_{}_{}\n'.format(title,mer,count))
                                out_handle.write('{}\n'.format(subseq))
                                existing.add(subseq)
                                count += 1  

def get_ts_gene(atlas_dir):

    bayests_cutoff = 0.3
    bayests = pd.read_csv(bayests_xy_path,sep='\t',index_col=0)
    ts_gene = bayests.loc[bayests['BayesTS']<bayests_cutoff,:].index.tolist()
    common = ts_gene

    try:
        deg = pd.read_csv(os.path.join(atlas_dir,'deg.txt'),sep='\t',index_col=1)
    except:
        pass
    else:
        deg = deg.loc[(deg['adjp']<0.05) & (deg['Log2(Fold Change)']>0.58),:].index.tolist()
        common = list(set(ts_gene).intersection(set(deg)))


    gene_lfc = pd.read_csv(os.path.join(atlas_dir,'gene_lfc.txt'),sep='\t',index_col=0)
    ensg2symbol = gene_lfc['gene_symbol'].to_dict()
    symbol2ensg = {v:k for k,v in ensg2symbol.items()}
    gene_lfc = gene_lfc.loc[(gene_lfc['median_tumor']>20) & (gene_lfc['median_tumor'] > gene_lfc['max_median_gtex']),:]
    real_common = list(set(gene_lfc.index).intersection(set(common)))


    # remove the manual bl
    manual_bl_ensg = [symbol2ensg[item] for item in manual_bl]
    real_common = list(set(real_common).difference(set(manual_bl_ensg)))

    # membrane
    mem = pd.read_csv('human_membrane_proteins_acc2ens.txt',sep='\t')
    mem = mem.loc[mem['Ens'].notna(),:]
    membrane_ensg = mem['Ens'].values.tolist()

    mem = pd.read_csv('./compartment/human_cell_membrane_protein_postdoc_final.txt',sep='\t')
    mem = mem.loc[mem['Subcellular location [CC]'].str.contains('anchor'),:]
    mem = mem.loc[mem['ensg'].notna(),:]
    membrane_ensg_add = mem['ensg'].values.tolist()

    membrane_ensg = set(membrane_ensg + membrane_ensg_add)

    real_common_membrane = []
    for item in real_common:
        if item in membrane_ensg:
            real_common_membrane.append(item)


    return real_common,real_common_membrane


def get_ts_splicing(atlas_dir,s):
    splicing_path = os.path.join(atlas_dir,'splicing_rec.txt')
    total_sample = s
    vc = pd.read_csv(splicing_path,sep='\t',index_col=0)
    cond = (vc['n_sample'] > total_sample * 0.2) & (vc['ave_count_normal'] < 1) & (vc['ave_count_tumor'] > 10) & (vc['logFC'] > 4) & (vc['has_known_ss'])
    vc = vc.loc[cond,:]

    return vc

def get_ts_intron_retention(atlas_dir,s):
    intron_path = os.path.join(atlas_dir,'intron_rec.txt')
    normal_recurrency = 0
    total_number = s
    final = pd.read_csv(intron_path,sep='\t',index_col=0)
    final = final.loc[final['cond'],:]
    final['normal'] = [literal_eval(item) if isinstance(item,str) else {} for item in final['normal']]
    final['normal_recurrency'] = [len(item) for item in final['normal']]
    final = final.loc[final['normal_recurrency'] <= normal_recurrency, :]
    final = final.loc[final['n_sample'] > 0.15*total_number,:]

    return final

def get_ts_erv(atlas_dir):
    good_erv = os.path.join(atlas_dir,'good_erv.txt')
    df = pd.read_csv(good_erv,sep='\t')

    return df

def get_ts_variant(atlas_dir,c):

    if c == 'RT':

        n_variant = 544

    else:

        mutation_rec_path = os.path.join(atlas_dir,'mutation_rec.txt')
        variants = pd.read_csv(mutation_rec_path,sep='\t')
        variants = variants.loc[(variants['ensg'].notna()) & (variants['ensg'].str.contains(r'^ENSG')), :]   
        n_variant = variants.shape[0]

    return n_variant

def get_ts_pathogen(c):
    fixed_dict = {
        'BRCA':1,
        'KIRC':0,
        'COAD':2,
        'STAD':5,
        'MESO':0,
        'LIHC':2,
        'ESCA':5,
        'CESC':3,
        'BLCA':0,
        'RT':0,
        'AML':0,
        'DLBC':0,
        'GBM':1,
        'NBL':1,
        'PAAD':0,
        'HNSC':3,
        'OV':4,
        'LUSC':0,
        'LUAD':0,
        'CHOL':0,
        'SKCM':1
    }

    return fixed_dict[c]

def get_ts_fusion(atlas_dir):

    try:
        df = pd.read_csv(os.path.join(atlas_dir,'fusion_recurrent.txt'),sep='\t',index_col=0)
    except:
        df = pd.read_csv(os.path.join(atlas_dir,'fusion_rec.txt'),sep='\t',index_col=0)

    return df


# help with membrane protein filter
df_list = []
for c,s in zip(cancers,n_samples):
    atlas_dir = os.path.join(root_atlas_dir,c)
    real_common, real_common_membrane = get_ts_gene(atlas_dir)
    s = pd.Series(real_common_membrane)
    df_list.append(s)
df = pd.concat(df_list,axis=0,keys=cancers)
df.to_csv('my_filter_membrane.txt',sep='\t')


# just generate for_safety_screen.txt
data = []
for c in cancers:
    final_path = os.path.join(root_atlas_dir,c,'antigen','0.05','final_enhanced.txt')
    final = pd.read_csv(final_path,sep='\t')
    cond = [False if ('[]' in item) and ('(\'HLA-' not in item) else True for item in final['presented_by_each_sample_hla']]
    final = final.loc[cond,:]
    data.append(final)
final = pd.concat(data,axis=0,keys=cancers).reset_index(level=-2).rename(columns={'level_0':'cancer'})
final.to_csv('for_safety_screen.txt',sep='\t',index=None)

# making db
df = pd.read_csv('for_safety_screen.txt',sep='\t')
with open('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/db_fasta/peptides.fasta','w') as f:
    for cancer,pep,typ in zip(df['cancer'],df['pep'],df['typ']):
        f.write('>query|{}|{}|{}\n{}\n'.format(pep,typ,cancer,pep))

# # post safety screen, add normal ribo, conduct I/L
# safety_screen_df = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/code/post_safety_screen.txt',sep='\t')

# # for mer in [8,9,10,11,12,13,14,15]:
# #     print(mer)
# #     chop_normal_pep_db(fasta_path=protein_path,
# #                        output_path='ensembl_protein_{}.fasta'.format(mer),mers=[mer],allow_duplicates=False)


# fasta_dic = {}
# for mer in [8,9,10,11,12,13,14,15]:
#     lis = []
#     with open('ensembl_protein_{}.fasta'.format(mer),'r') as in_handle:
#         for title,seq in SimpleFastaParser(in_handle):
#             lis.append(seq)
#     lis = set(lis)
#     fasta_dic[mer] = lis

# col = []
# for pep,typ in tqdm(zip(safety_screen_df['pep'],safety_screen_df['typ']),total=safety_screen_df.shape[0]):
#     if typ == 'self_gene':
#         col.append(False)
#     else:
#         if 'I' in pep or 'L' in pep:
#             options = []
#             for c in pep:
#                 if c == 'I' or c == 'L':
#                     options.append(['I','L'])
#                 else:
#                     options.append([c])
#             expanded = [''.join(p) for p in itertools.product(*options)]
#             l = len(pep)
#             lis = fasta_dic[l]
#             flag = False
#             for item in expanded:
#                 if item in lis:
#                     flag = True
#                     break
#             col.append(flag)
                
#         else:
#             col.append(False)

# safety_screen_df['is_ambiguous_IL'] = col

# old_dir = os.getcwd()
# os.chdir('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/normal_ribo/result')
# all_fasta = subprocess.run('find . -type f -name "riborf.fasta"',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
# anno_ribo = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/normal_ribo/ely_et_al_anno.txt',sep='\t',index_col=0)
# srr2tissue = anno_ribo['Cell_type'].to_dict()
# meta_dict = {}
# for f in all_fasta:
#     _,srr,_,_ = f.split('/')
#     micro_dict = {}
#     with open(f,'r') as in_handle:
#         for title,seq in SimpleFastaParser(in_handle):
#             micro_dict[title] = seq
#     tissue = srr2tissue[srr]
#     meta_dict['{};{}'.format(srr,tissue)] = micro_dict
# os.chdir(old_dir)

# condensed_dict = {}   # tissue:[seq1,seq2,...]
# for k,v in meta_dict.items():
#     srr,tissue = k.split(';')
#     condensed_dict.setdefault(tissue,[]).extend(list(v.values()))
# condensed_dict = {k:set(v) for k,v in condensed_dict.items()}

# col = []
# not_in_normal_ribo = []
# for row in tqdm(safety_screen_df.itertuples(),total=safety_screen_df.shape[0]):
#     if row.typ == 'nuORF':
#         normal_tissues = []
#         pep = row.pep
#         for t,ss in condensed_dict.items():
#             for s in ss:
#                 if pep in s:
#                     normal_tissues.append(t)
#                     break
#         col.append(','.join(normal_tissues))
#         if len(normal_tissues) > 0:
#             not_in_normal_ribo.append(False)
#         else:
#             not_in_normal_ribo.append(True)
#     else:
#         col.append(None)
#         not_in_normal_ribo.append(True)
# safety_screen_df['normal_ribo'] = col
# safety_screen_df['not_in_normal_ribo'] = not_in_normal_ribo
# safety_screen_df.to_csv('post_safety_screen_add_ribo.txt',sep='\t',index=None)


safety_screen_df = pd.read_csv('post_safety_screen_add_ribo.txt',sep='\t')
safety_screen_bl = set(safety_screen_df.loc[(~safety_screen_df['cond_stringent']) | (~safety_screen_df['not_in_normal_ribo']) | (safety_screen_df['is_ambiguous_IL']),:]['pep'].values)

data = []
for c in cancers:
    final_path = os.path.join(root_atlas_dir,c,'antigen','0.01','final_enhanced.txt')
    final = pd.read_csv(final_path,sep='\t')
    cond = [False if ('[]' in item) and ('(\'HLA-' not in item) else True for item in final['presented_by_each_sample_hla']]
    final = final.loc[cond,:]
    data.append(final)
final = pd.concat(data,axis=0,keys=cancers).reset_index(level=-2).rename(columns={'level_0':'cancer'})
final = final.loc[~final['pep'].isin(safety_screen_bl),:]

self_df = final.loc[(final['typ']=='self_gene') & (final['unique']!=False),:] 
self_df = self_df.loc[~self_df['gene_symbol'].isin(manual_bl),:]

self_translate_te_df = final.loc[final['typ']=='ERV',:]
real_autonomy_check = pd.read_csv('splicing_ir_dic/final.txt',sep='\t')
cond = []
for row in real_autonomy_check.itertuples():
    if row.not_has_ss and row.not_in_ir:
        cond.append(True)
    else:
        cond.append(False)
real_autonomy = set(real_autonomy_check.loc[cond,:]['pep'].values)
self_translate_te_df = self_translate_te_df.loc[self_translate_te_df['pep'].isin(real_autonomy),:]
self_translate_te_df['typ'] = np.full(shape=self_translate_te_df.shape[0],fill_value='self_translate_te')

te_chimeric_df = final.loc[final['typ']=='TE_chimeric_transcript',:]
original_all_self_translate =  final.loc[final['typ']=='ERV',:]
reclassified_te_chimeric = original_all_self_translate.loc[~original_all_self_translate['pep'].isin(real_autonomy),:]
te_chimeric_df = pd.concat([te_chimeric_df,reclassified_te_chimeric],axis=0)
te_chimeric_df['typ'] = np.full(shape=te_chimeric_df.shape[0],fill_value='TE_chimeric_transcript')

splicing_df = final.loc[final['typ']=='splicing',:]
cond = [False if 'nc|ENSG00000100146|P56693|SOX10' in item else True for item in splicing_df['source']]
splicing_df = splicing_df.loc[cond,:]

nuorf_df = final.loc[final['typ']=='nuORF',:]

variant_df = final.loc[final['typ']=='variant',:]

fusion_df = final.loc[final['typ']=='fusion',:]
fusion_df = fusion_df.loc[~fusion_df['source'].str.contains('nc'),:]

ir_df = final.loc[final['typ']=='intron_retention',:]

pathogen_df = final.loc[final['typ']=='pathogen',:]

patent_df = pd.concat([self_df,self_translate_te_df,te_chimeric_df,splicing_df,nuorf_df,variant_df,fusion_df,ir_df,pathogen_df])
patent_df.to_csv('final_all_ts_antigens.txt',sep='\t',index=None)
print(len(patent_df['pep'].unique()))

# compare with iedb
iedb_df = pd.read_csv('all_epitope_no_b_human_linear_mhc_i.tsv',sep='\t')
iedb_pep = set(iedb_df['Epitope - Name'].values)
immunoverse_pep = set(patent_df['pep'].values)
novel = immunoverse_pep.difference(iedb_pep)
novel_rate = len(novel)/len(immunoverse_pep)
print(novel_rate)

data = []
for pep,patent_sub_df in patent_df.groupby(by='pep'):
    patent_sub_df.sort_values(by='highest_score',inplace=True,ascending=False)
    all_c = ','.join(patent_sub_df['cancer'].values.tolist()[:3])
    tmp1 = [item[0].replace('HLA-','').replace('*','').replace(':','') for item in literal_eval(patent_sub_df['additional_query'].iloc[0])]
    tmp2 = [item[2] for item in literal_eval(patent_sub_df['additional_query'].iloc[0])]
    tmp = sorted(zip(tmp1,tmp2),key=lambda x:x[1])[:3]
    all_hla = ','.join([item[0] for item in tmp])
    data.append((pep,all_c,all_hla))
patent_df_final = pd.DataFrame.from_records(data=data,columns=['peptide','indication','HLA'])
patent_df_final.to_csv('patent_df_final.txt',sep='\t',index=None)


# plot peptide overview, order will be gene, splicing, self_TE, chimera_TE, IR, pathogen, fusion, variant, lncRNA, pseudogene, cryptic ORF
data = []
# self_gene
tmp_dic = {}
for c_,sub_df in self_df.groupby(by='cancer'):
    tmp_dic[c_] = sub_df.shape[0]
tmp_data = []
for c in cancers:
    v = tmp_dic.get(c,None)
    if v is None:
        tmp_data.append(0)
    else:
        tmp_data.append(v)
data.append(tmp_data)
# splicing
tmp_dic = {}
for c_,sub_df in splicing_df.groupby(by='cancer'):
    tmp_dic[c_] = sub_df.shape[0]
tmp_data = []
for c in cancers:
    v = tmp_dic.get(c,None)
    if v is None:
        tmp_data.append(0)
    else:
        tmp_data.append(v)
data.append(tmp_data)
# self_TE
tmp_dic = {}
for c_,sub_df in self_translate_te_df.groupby(by='cancer'):
    tmp_dic[c_] = sub_df.shape[0]
tmp_data = []
for c in cancers:
    v = tmp_dic.get(c,None)
    if v is None:
        tmp_data.append(0)
    else:
        tmp_data.append(v)
data.append(tmp_data)
# chimera_TE, need to change
tmp_dic = {}
for c_,sub_df in te_chimeric_df.groupby(by='cancer'):
    tmp_dic[c_] = sub_df.shape[0]
tmp_data = []
for c in cancers:
    v = tmp_dic.get(c,None)
    if v is None:
        tmp_data.append(0)
    else:
        tmp_data.append(v)
data.append(tmp_data)
# IR
tmp_dic = {}
for c_,sub_df in ir_df.groupby(by='cancer'):
    tmp_dic[c_] = sub_df.shape[0]
tmp_data = []
for c in cancers:
    v = tmp_dic.get(c,None)
    if v is None:
        tmp_data.append(0)
    else:
        tmp_data.append(v)
data.append(tmp_data)
# pathogen
tmp_dic = {}
for c_,sub_df in pathogen_df.groupby(by='cancer'):
    tmp_dic[c_] = sub_df.shape[0]
tmp_data = []
for c in cancers:
    v = tmp_dic.get(c,None)
    if v is None:
        tmp_data.append(0)
    else:
        tmp_data.append(v)
data.append(tmp_data)
# fusion
tmp_dic = {}
for c_,sub_df in fusion_df.groupby(by='cancer'):
    tmp_dic[c_] = sub_df.shape[0]
tmp_data = []
for c in cancers:
    v = tmp_dic.get(c,None)
    if v is None:
        tmp_data.append(0)
    else:
        tmp_data.append(v)
data.append(tmp_data)
# variant
tmp_dic = {}
for c_,sub_df in variant_df.groupby(by='cancer'):
    tmp_dic[c_] = sub_df.shape[0]
tmp_data = []
for c in cancers:
    v = tmp_dic.get(c,None)
    if v is None:
        tmp_data.append(0)
    else:
        tmp_data.append(v)
data.append(tmp_data)
# lncRNA, pseudogene, cryptic ORF
lncRNA_df = nuorf_df.loc[nuorf_df['nuorf_type']=='lncRNA',:]
pseudo_df =  nuorf_df.loc[nuorf_df['nuorf_type']=='Pseudogene',:]
cryptic_df = nuorf_df.loc[(nuorf_df['nuorf_type']!='Pseudogene') & (nuorf_df['nuorf_type']!='lncRNA'),:]
prop = 1
# lncRNA
tmp_dic = {}
for c_,sub_df in lncRNA_df.groupby(by='cancer'):
    tmp_dic[c_] = sub_df.shape[0]
tmp_data = []
for c in cancers:
    v = tmp_dic.get(c,None)
    if v is None:
        tmp_data.append(0)
    else:
        tmp_data.append(round(v*prop))
data.append(tmp_data)
# pseudogene
tmp_dic = {}
for c_,sub_df in pseudo_df.groupby(by='cancer'):
    tmp_dic[c_] = sub_df.shape[0]
tmp_data = []
for c in cancers:
    v = tmp_dic.get(c,None)
    if v is None:
        tmp_data.append(0)
    else:
        tmp_data.append(round(v*prop))
data.append(tmp_data)
# cryptic ORF
tmp_dic = {}
for c_,sub_df in cryptic_df.groupby(by='cancer'):
    tmp_dic[c_] = sub_df.shape[0]
tmp_data = []
for c in cancers:
    v = tmp_dic.get(c,None)
    if v is None:
        tmp_data.append(0)
    else:
        tmp_data.append(round(v*prop))
data.append(tmp_data)

plot_df = pd.DataFrame(index=['gene','splicing','self_translate_TE','chimera_TE','intron_retention','pathogen','fusion','variant','lncRNA','pseudogene','cryptic ORF'],
                       columns=cancers,data=np.array(data))

fig,ax = plt.subplots(figsize=(15,15))
sns.heatmap(plot_df,annot=True,linewidth=0.5,fmt='g',annot_kws={"fontsize": 8},cmap='Blues',square=True)
plt.savefig('ts_antigen_overview.pdf',bbox_inches='tight')
plt.close()


hla_dic = {}
dic1 = {}
dic2 = {}
total_immuno = 0  
total_immuno_bio = 0
root_immuno_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/immunopeptidome'
lung = None
with pd.ExcelWriter('all_immuno_meta.xlsx') as writer:
    for c in cancers:
        f = os.path.join(root_immuno_dir,cancers2immuno[c],'metadata.txt')
        df = pd.read_csv(f,sep='\t',index_col=0)
        df = df.loc[~df.index.str.startswith('phase2_'),:]
        df = df.loc[df['special_note'].isna(),:]
        if c == 'LUAD':
            lung = df
        total_immuno += df.shape[0]
        total_immuno_bio += len(df['biology'].unique())
        dic1[c] = df.shape[0]
        dic2[c] = len(df['biology'].unique())
        df.to_excel(writer,sheet_name=c)
        hlas = []
        for item in df['HLA']:
            if isinstance(item,str):
                hlas.extend(list(set(item.split(','))))
        values,counts = np.unique(hlas,return_counts=True)
        hla_dic[c] = {v:c_ for v,c_ in zip(values,counts)}
# do not double count lung cancer
print(total_immuno - lung.shape[0])
print(total_immuno_bio - len(lung['biology'].unique()))


hla_data = []
all_hla = []
for k,v in hla_dic.items():
    all_hla.extend(list(v.keys()))
all_hla = list(set(all_hla))
all_hla = sorted(all_hla)
for c in cancers:
    tmp = []
    for hla in all_hla:
        tmp.append(hla_dic[c].get(hla,0))
    hla_data.append(tmp)

hla_df = pd.DataFrame(index=cancers,columns=all_hla,data=hla_data)
reformatted_hla = ['HLA-' + item.replace(':','').replace('*','') for item in hla_df.columns]
hla_df.columns = reformatted_hla
freq_df = pd.read_csv('/gpfs/data/yarmarkovichlab/medulloblastoma/neoverse_folder/NeoVerse_final_output_new/antigens/US_HLA_frequency.csv',sep=',',index_col=0)
freq_dic = freq_df['Percent US population'].to_dict()
ori_array = [tuple([item[4:] for item in hla_df.columns.tolist()]),
             tuple([freq_dic.get(item,0) for item in hla_df.columns]),
             tuple([item[4] for item in hla_df.columns])]
mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
hla_df.columns = mi

cancer2coverage = {}
for c in cancers:
    s = hla_df.loc[c]
    s = s.loc[s>0]
    freqs = s.index.to_frame(index=False)[1].values.tolist()
    neg = 1
    for f in freqs:
        neg *= (1-f)
    coverage = 1-neg
    cancer2coverage[c] = coverage
ori_array = [tuple(hla_df.index.tolist()),
             tuple([cancer2coverage[item] for item in hla_df.index])]
mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
hla_df.index = mi
hla_df.to_csv('hla_df.txt',sep='\t')

selected_col = [col for col in hla_df.columns if col[1] > 0.01]
hla_df_freq = hla_df.loc[:,selected_col]
hla_df_freq.to_csv('hla_df_freq.txt',sep='\t')

hla_df_freq = hla_df_freq.mask(hla_df_freq>0,1)
hla_df_freq.to_csv('hla_df_freq_binary.txt',sep='\t')


fig,ax = plt.subplots()
bars = ax.bar(x=[i*3 for i in range(len(dic1))],height=list(dic1.values()),width=1)
for bar in bars:
    yval = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2, yval + 0.5, yval, ha='center', va='bottom')
bars = ax.bar(x=[i*3+1 for i in range(len(dic2))],height=list(dic2.values()),width=1)
for bar in bars:
    yval = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2, yval + 0.5, yval, ha='center', va='bottom')
ax.set_xticks([i*3+0.5 for i in range(len(dic2))])
ax.set_xticklabels(list(dic1.keys()),rotation=60)
ax.set_ylabel('Number of sample')
plt.savefig('figs1_stat_immuno.pdf',bbox_inches='tight')
plt.close()


dic = {}
total_rna = 0  
with pd.ExcelWriter('all_rna_manifest.xlsx') as writer:
    for c in cancers:
        print(c)
        cmd = 'find {} -type f -name "manifest_*_*.tsv"'.format(os.path.join(root_atlas_dir,c))
        possible_f = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1][0]
        df = pd.read_csv(possible_f,sep='\t',index_col=0)
        df = df.loc[df['sample_id'].notna(),:]
        tid = [int(item.split('-')[-1][:-1]) for item in df['sample_id']]
        cond = [False if item == 10 or item ==11 else True for item in tid]
        df = df.loc[cond,:]
        total_rna += df.shape[0]
        dic[c] = df.shape[0]
        df.to_excel(writer,sheet_name=c)
print(total_rna)

fig,ax = plt.subplots()
bars = ax.bar(x=np.arange(len(dic)),height=list(dic.values()))
for bar in bars:
    yval = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2, yval + 0.5, yval, ha='center', va='bottom')
ax.set_xticks(np.arange(len(dic)))
ax.set_xticklabels(list(dic.keys()),rotation=60)
ax.set_ylabel('Number of sample')
plt.savefig('figs1_stat_rna.pdf',bbox_inches='tight')
plt.close()

# single-cell stats
scem = pd.read_csv('/gpfs/data/yarmarkovichlab/cancerSCEM2.0/CancerSCEM-Browse-Table.csv',sep=',')
scem = scem.loc[scem['Sample Type']=='Tumour',:]
old_dir = os.getcwd()
os.chdir('/gpfs/data/yarmarkovichlab/Frank/logic_finder_v2')
all_h5ad = subprocess.run('find . -mindepth 3 -type f -name "*.h5ad"',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
os.chdir(old_dir)
valid_cancers = ['LIHC','BLCA','BRCA','SKCM','HNSC','ESCA','KIRC','STAD','LUSC','LUAD','COAD','GBM','OV','PAAD']
sample_list = []
for h5ad in all_h5ad:
    if len(h5ad.split('/')) == 4:
        _,d,_,f = h5ad.split('/')
        dc = d.split('_')[1]
        if dc in valid_cancers:
            sample_list.append(f.split('.h5ad')[0])
print(len(sample_list))
scem = scem.loc[scem['Sample ID'].isin(sample_list),:]
scem.drop_duplicates(subset='Sample ID',inplace=True)
scem.drop(columns=['Transcriptome Analysis','Metabolic Profile'],inplace=True)
scem.to_csv('scem.txt',sep='\t',index=None)


# tumor specific event 
data = []
n_ts_membrane = []
n_ts_total = []
for c,s in zip(cancers,n_samples):
    data_column = []
    atlas_dir = os.path.join(root_atlas_dir,c)
    real_common, real_common_membrane = get_ts_gene(atlas_dir)
    data_column.append(len(real_common))
    n_ts_membrane.append(len(real_common_membrane))

    vc = get_ts_splicing(atlas_dir,s)
    data_column.append(vc.shape[0])

    final = get_ts_intron_retention(atlas_dir,s)
    data_column.append(final.shape[0])

    df = get_ts_erv(atlas_dir)
    data_column.append(df.shape[0])

    n_variant = get_ts_variant(atlas_dir,c)
    data_column.append(n_variant)

    number = get_ts_pathogen(c)
    data_column.append(number)

    df = get_ts_fusion(atlas_dir)
    data_column.append(df.shape[0])

    total = sum(data_column)
    n_ts_total.append(total)
    data.append(data_column)



plot_df = pd.DataFrame(index=['gene','splicing','intron','TE','variant','pathogen','fusion'],columns=cancers,data=np.transpose(np.array(data)))

fig,ax = plt.subplots(figsize=(15,10))
sns.heatmap(plot_df,annot=True,linewidth=0.5,fmt='g',vmax=15000,annot_kws={"fontsize": 5},cmap='Blues',square=True)
plt.savefig('ts_event_overview.pdf',bbox_inches='tight')
plt.close()


# for ts gene, coverage for patients, also consider all antigens
df = pd.read_csv('final_all_ts_antigens.txt',sep='\t')
# so it really does not change that much, this list versus the 66 that I used for enrichment
pan_cancer_ensgs = pd.read_csv('pan_cancer_cluster.txt',sep='\t',header=None)[0].values.tolist()
# remove pan-cancer gene, because they are not ideal as single target
cond = []
for item1,item2 in zip(df['typ'],df['ensgs']):
    if item1 == 'self_gene' and item2 in pan_cancer_ensgs:
        cond.append(False)
    else:
        cond.append(True)
df = df.loc[cond,:]
# add g_prop column
ts_ensg = df.loc[df['typ']=='self_gene',:]['ensgs'].values
self_gene_dic = {}
self_gene_dic_overall = {}
for c in cancers:
    path = os.path.join(root_atlas_dir,c,'gene_tpm.txt')
    exp = pd.read_csv(path,sep='\t',index_col=0)
    exp = exp.loc[ts_ensg,:]
    exp = exp.values
    exp = (exp > 20).astype(int)
    props = np.sum(exp,axis=1) / exp.shape[1]
    mapping = {} # gene specific prop
    for i1,i2 in zip(ts_ensg,props):
        mapping[i1] = i2
    self_gene_dic[c] = mapping
    prop = np.count_nonzero(np.any(exp,axis=0)) / exp.shape[1]
    self_gene_dic_overall[c] = prop

col = []
for row in df.itertuples():
    if row.typ == 'self_gene':
        col.append(self_gene_dic[row.cancer][row.ensgs])
    elif (row.typ == 'splicing' or row.typ == 'TE_chimeric_transcript') and row.unique:
        try:
            n = float(row.source.split('|')[1]) 
            p = n / n_samples[cancers.index(row.cancer)]
        except:
            p = 0
        col.append(p)
    elif row.typ == 'self_translate_te' and row.unique:
        try:
            n = float(row.source.split('|')[5])
            p = n / n_samples[cancers.index(row.cancer)]
        except:
            p = 0
        col.append(p)
    elif row.typ == 'variant' and row.unique:
        try:
            n = float(row.source.split('|')[2])
            p = n / n_samples[cancers.index(row.cancer)]
        except:
            p = 0
        col.append(p)
    elif row.typ == 'intron_retention' and row.unique:
        try:
            n = float(row.source.split('|')[1])
            p = n / n_samples[cancers.index(row.cancer)]
        except:
            p = 0
        col.append(p)
    else:
        col.append(0)

df['g_prop'] = col


data = []
for c in cancers:
    # now consider hla
    if c in ['RT','NBL']:
        dic = pd.read_csv('/gpfs/data/yarmarkovichlab/medulloblastoma/neoverse_folder/NeoVerse_final_output_new/antigens/US_HLA_frequency.csv',sep=',',index_col=0)['Percent US population'].to_dict()
        dic = {k.replace('HLA-',''):v for k,v in dic.items()}
    else:
        dic = {}
        hla_path = os.path.join(root_atlas_dir,c,'hla_types.txt')
        hla = pd.read_csv(hla_path,sep='\t',index_col=0)
        tmp = hla.loc[:,['HLAA1','HLAA2']].values.flatten().tolist()
        values,counts = np.unique(tmp,return_counts=True)
        for v,c_ in zip(values,counts):
            if v.startswith('A'):
                dic[v] = c_/len(tmp)
        tmp = hla.loc[:,['HLAB1','HLAB2']].values.flatten().tolist()
        values,counts = np.unique(tmp,return_counts=True)
        for v,c_ in zip(values,counts):
            if v.startswith('B'):
                dic[v] = c_/len(tmp)
        tmp = hla.loc[:,['HLAC1','HLAC2']].values.flatten().tolist()
        values,counts = np.unique(tmp,return_counts=True)
        for v,c_ in zip(values,counts):
            if v.startswith('C'):
                dic[v] = c_/len(tmp)
    # start to look for antigen 
    df_sub = df.loc[(df['cancer']==c) & (df['g_prop']!=0),:]
    mapping_norm_sb = {}
    mapping_norm_wb = {}
    for source,sub_df in tqdm(df_sub.groupby(by='source')):
        g_prop = sub_df['g_prop'].iloc[0]
        all_query = []
        sb_hla = []
        wb_hla = []
        for item in sub_df['additional_query']:
            item = literal_eval(item)
            all_query.extend(item)
        for item in all_query:
            wb_hla.append(item[0])
            if item[3] == 'SB':
                sb_hla.append(item[0])
        sb_hla = list(set(sb_hla))
        wb_hla = list(set(wb_hla))

        tmp = [dic.get(h.replace('HLA-','').replace('*','').replace(':',''),0) for h in sb_hla]
        h_prop = 1
        for item in tmp:
            h_prop *= (1-item)
        norm_prop = g_prop * (1-h_prop)
        mapping_norm_sb[source] = norm_prop

        tmp = [dic.get(h.replace('HLA-','').replace('*','').replace(':',''),0) for h in wb_hla]
        h_prop = 1
        for item in tmp:
            h_prop *= (1-item)
        norm_prop = g_prop * (1-h_prop)
        mapping_norm_wb[source] = norm_prop

    final_prop = 1
    for k,v in mapping_norm_sb.items():
        if not math.isnan(v):
            final_prop *= (1-v)
    final_sb_prop = 1-final_prop

    final_prop = 1
    for k,v in mapping_norm_wb.items():
        if not math.isnan(v):
            final_prop *= (1-v)
    final_wb_prop = 1-final_prop

    data.append((self_gene_dic_overall[c],'gene_prop',c))
    data.append((final_wb_prop,'wb_prop',c))
    data.append((final_sb_prop,'sb_prop',c))

plot_df = pd.DataFrame.from_records(data,columns=['value','category','cancer'])

plot_df_now = plot_df.loc[plot_df['category']=='sb_prop',:].sort_values(by='value',ascending=False)
cancer2n_sample = pd.Series(index=cancers,data=n_samples).to_dict()
plot_df_now['n_sample'] = plot_df_now['cancer'].map(cancer2n_sample).values
coverage = 0
for i1,i2 in zip(plot_df_now['n_sample'],plot_df_now['value']):
    coverage += round(i1*i2)
print(coverage/sum(n_samples))  













    

