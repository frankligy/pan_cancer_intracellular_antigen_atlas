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
    df = df.loc[df['logfc']>5,:]

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
        'CESC':2,
        'BLCA':0,
        'RT':0,
        'AML':0,
        'DLBC':0,
        'GBM':1,
        'NBL':1,
        'PAAD':0,
        'HNSC':3,
        'OV':3,
        'LUSC':0,
        'LUAD':0,
        'CHOL':0,
        'SKCM':1
    }

    return fixed_dict[c]

def get_ts_fusion(atlas_dir):

    df = pd.read_csv(os.path.join(atlas_dir,'fusion_recurrent.txt'),sep='\t',index_col=0)

    return df


# help with membrane protein filter
# df_list = []
# for c,s in zip(cancers,n_samples):
#     atlas_dir = os.path.join(root_atlas_dir,c)
#     real_common, real_common_membrane = get_ts_gene(atlas_dir)
#     s = pd.Series(real_common_membrane)
#     df_list.append(s)
# df = pd.concat(df_list,axis=0,keys=cancers)
# df.to_csv('my_filter_membrane.txt',sep='\t')


# safety_screen_df = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/code/post_safety_screen.txt',sep='\t')
# safety_screen_df = safety_screen_df.loc[~safety_screen_df['cond_stringent'],:]
# safety_screen_bl = list(set(safety_screen_df['pep'].values.tolist()))

# collage the meta and get number
total_antigen = 0  # 27451
self_df = pd.read_csv('ts_final.txt',sep='\t')
# self_df = self_df.loc[~self_df['pep'].isin(safety_screen_bl),:]
total_antigen += len(set(self_df['pep'].values.tolist()))

self_translate_te_df = pd.read_csv('ts_te_antigen.txt',sep='\t')  # remember, after safety screen, do autonomy check and update 
# real_autonomy_check = pd.read_csv('splicing_ir_dic/final.txt',sep='\t')
# real_autonomy = set(real_autonomy_check.loc[real_autonomy_check['not_has_ss'] & real_autonomy_check['not_in_ir'],:]['pep'].values)
# self_translate_te_df = self_translate_te_df.loc[self_translate_te_df['pep'].isin(real_autonomy),:]
# te_all_df = pd.read_csv('te_all_antigens.txt',sep='\t')
# orf2_taa = te_all_df.loc[te_all_df['source'].str.contains('L1_ORF2'),:]
# self_translate_te_df = pd.concat([self_translate_te_df,orf2_taa],axis=0)
# self_translate_te_df['typ'] = np.full(shape=self_translate_te_df.shape[0],fill_value='self_translate_te')
# self_translate_te_df = self_translate_te_df.loc[~self_translate_te_df['pep'].isin(safety_screen_bl),:]
total_antigen += len(set(self_translate_te_df['pep'].values.tolist()))

te_chimeric_df = pd.read_csv('te_all_antigens.txt',sep='\t')
te_chimeric_df = te_chimeric_df.loc[te_chimeric_df['typ']=='TE_chimeric_transcript',:]
# original_all_self_translate =  pd.read_csv('ts_te_antigen.txt',sep='\t')
# reclassified_te_chimeric = original_all_self_translate.loc[~original_all_self_translate['pep'].isin(real_autonomy),:]
# te_chimeric_df = pd.concat([te_chimeric_df,reclassified_te_chimeric],axis=0)
# te_chimeric_df['typ'] = np.full(shape=te_chimeric_df.shape[0],fill_value='TE_chimeric_transcript')
# te_chimeric_df = te_chimeric_df.loc[~te_chimeric_df['pep'].isin(safety_screen_bl),:]
total_antigen += len(set(te_chimeric_df['pep'].values.tolist()))

splicing_df = pd.read_csv('all_splicing.txt',sep='\t')
# splicing_df = splicing_df.loc[~splicing_df['pep'].isin(safety_screen_bl),:]
total_antigen += len(set(splicing_df['pep'].values.tolist()))

nuorf_df = pd.read_csv('all_nuorf.txt',sep='\t')
# nuorf_df = nuorf_df.loc[~nuorf_df['pep'].isin(safety_screen_bl),:]
total_antigen += len(set(nuorf_df['pep'].values.tolist()))

variant_df = pd.read_csv('all_variants.txt',sep='\t')
# variant_df = variant_df.loc[~variant_df['pep'].isin(safety_screen_bl),:]
total_antigen += len(set(variant_df['pep'].values.tolist()))

fusion_df = pd.read_csv('all_fusion.txt',sep='\t')
fusion_df = fusion_df.loc[~fusion_df['source'].str.contains('nc'),:]
# fusion_df = fusion_df.loc[~fusion_df['pep'].isin(safety_screen_bl),:]
total_antigen += len(set(fusion_df['pep'].values.tolist()))

ir_df = pd.read_csv('all_ir.txt',sep='\t')
# ir_df = ir_df.loc[~ir_df['pep'].isin(safety_screen_bl),:]
total_antigen += len(set(ir_df['pep'].values.tolist()))

pathogen_df = pd.read_csv('all_pathogen.txt',sep='\t')
# pathogen_df = pathogen_df.loc[~pathogen_df['pep'].isin(safety_screen_bl),:]
total_antigen += len(set(pathogen_df['pep'].values.tolist()))

patent_df = pd.concat([self_df,self_translate_te_df,te_chimeric_df,splicing_df,nuorf_df,variant_df,fusion_df,ir_df,pathogen_df])
patent_df.to_csv('for_safety_screen.txt',sep='\t',index=None)   # you can generate for safety screen as well
sys.exit('stop')
# patent_df.to_csv('final_all_ts_antigens.txt',sep='\t',index=None)

# data = []
# for pep,patent_sub_df in patent_df.groupby(by='pep'):
#     patent_sub_df.sort_values(by='highest_score',inplace=True,ascending=False)
#     all_c = ','.join(patent_sub_df['cancer'].values.tolist()[:3])
#     tmp1 = [item[0].replace('HLA-','').replace('*','').replace(':','') for item in literal_eval(patent_sub_df['additional_query'].iloc[0])]
#     tmp2 = [item[2] for item in literal_eval(patent_sub_df['additional_query'].iloc[0])]
#     tmp = sorted(zip(tmp1,tmp2),key=lambda x:x[1])[:3]
#     all_hla = ','.join([item[0] for item in tmp])
#     data.append((pep,all_c,all_hla))
# patent_df_final = pd.DataFrame.from_records(data=data,columns=['peptide','indication','HLA'])
# patent_df_final.to_csv('patent_df_final.txt',sep='\t',index=None)

# print(total_antigen)



# # plot peptide overview, order will be gene, splicing, self_TE, chimera_TE, IR, pathogen, fusion, variant, lncRNA, pseudogene, cryptic ORF
# data = []
# # self_gene
# tmp_dic = {}
# for c_,sub_df in self_df.groupby(by='cancer'):
#     tmp_dic[c_] = sub_df.shape[0]
# tmp_data = []
# for c in cancers:
#     v = tmp_dic.get(c,None)
#     if v is None:
#         tmp_data.append(0)
#     else:
#         tmp_data.append(v)
# data.append(tmp_data)
# # splicing
# tmp_dic = {}
# for c_,sub_df in splicing_df.groupby(by='cancer'):
#     tmp_dic[c_] = sub_df.shape[0]
# tmp_data = []
# for c in cancers:
#     v = tmp_dic.get(c,None)
#     if v is None:
#         tmp_data.append(0)
#     else:
#         tmp_data.append(v)
# data.append(tmp_data)
# # self_TE
# tmp_dic = {}
# for c_,sub_df in self_translate_te_df.groupby(by='cancer'):
#     tmp_dic[c_] = sub_df.shape[0]
# tmp_data = []
# for c in cancers:
#     v = tmp_dic.get(c,None)
#     if v is None:
#         tmp_data.append(0)
#     else:
#         tmp_data.append(v)
# data.append(tmp_data)
# # chimera_TE, need to change
# tmp_dic = {}
# for c_,sub_df in te_chimeric_df.groupby(by='cancer'):
#     tmp_dic[c_] = sub_df.shape[0]
# tmp_data = []
# for c in cancers:
#     v = tmp_dic.get(c,None)
#     if v is None:
#         tmp_data.append(0)
#     else:
#         tmp_data.append(v)
# data.append(tmp_data)
# # IR
# tmp_dic = {}
# for c_,sub_df in ir_df.groupby(by='cancer'):
#     tmp_dic[c_] = sub_df.shape[0]
# tmp_data = []
# for c in cancers:
#     v = tmp_dic.get(c,None)
#     if v is None:
#         tmp_data.append(0)
#     else:
#         tmp_data.append(v)
# data.append(tmp_data)
# # pathogen
# tmp_dic = {}
# for c_,sub_df in pathogen_df.groupby(by='cancer'):
#     tmp_dic[c_] = sub_df.shape[0]
# tmp_data = []
# for c in cancers:
#     v = tmp_dic.get(c,None)
#     if v is None:
#         tmp_data.append(0)
#     else:
#         tmp_data.append(v)
# data.append(tmp_data)
# # fusion
# tmp_dic = {}
# for c_,sub_df in fusion_df.groupby(by='cancer'):
#     tmp_dic[c_] = sub_df.shape[0]
# tmp_data = []
# for c in cancers:
#     v = tmp_dic.get(c,None)
#     if v is None:
#         tmp_data.append(0)
#     else:
#         tmp_data.append(v)
# data.append(tmp_data)
# # variant
# tmp_dic = {}
# for c_,sub_df in variant_df.groupby(by='cancer'):
#     tmp_dic[c_] = sub_df.shape[0]
# tmp_data = []
# for c in cancers:
#     v = tmp_dic.get(c,None)
#     if v is None:
#         tmp_data.append(0)
#     else:
#         tmp_data.append(v)
# data.append(tmp_data)
# # lncRNA, pseudogene, cryptic ORF
# lncRNA_df = nuorf_df.loc[nuorf_df['nuorf_type']=='lncRNA',:]
# pseudo_df =  nuorf_df.loc[nuorf_df['nuorf_type']=='Pseudogene',:]
# cryptic_df = nuorf_df.loc[(nuorf_df['nuorf_type']!='Pseudogene') & (nuorf_df['nuorf_type']!='lncRNA'),:]
# prop = 1
# # lncRNA
# tmp_dic = {}
# for c_,sub_df in lncRNA_df.groupby(by='cancer'):
#     tmp_dic[c_] = sub_df.shape[0]
# tmp_data = []
# for c in cancers:
#     v = tmp_dic.get(c,None)
#     if v is None:
#         tmp_data.append(0)
#     else:
#         tmp_data.append(round(v*prop))
# data.append(tmp_data)
# # pseudogene
# tmp_dic = {}
# for c_,sub_df in pseudo_df.groupby(by='cancer'):
#     tmp_dic[c_] = sub_df.shape[0]
# tmp_data = []
# for c in cancers:
#     v = tmp_dic.get(c,None)
#     if v is None:
#         tmp_data.append(0)
#     else:
#         tmp_data.append(round(v*prop))
# data.append(tmp_data)
# # cryptic ORF
# tmp_dic = {}
# for c_,sub_df in cryptic_df.groupby(by='cancer'):
#     tmp_dic[c_] = sub_df.shape[0]
# tmp_data = []
# for c in cancers:
#     v = tmp_dic.get(c,None)
#     if v is None:
#         tmp_data.append(0)
#     else:
#         tmp_data.append(round(v*prop))
# data.append(tmp_data)

# plot_df = pd.DataFrame(index=['gene','splicing','self_translate_TE','chimera_TE','intron_retention','pathogen','fusion','variant','lncRNA','pseudogene','cryptic ORF'],
#                        columns=cancers,data=np.array(data))

# fig,ax = plt.subplots(figsize=(15,15))
# sns.heatmap(plot_df,annot=True,linewidth=0.5,fmt='g',annot_kws={"fontsize": 8},cmap='Blues',square=True)
# plt.savefig('ts_antigen_overview.pdf',bbox_inches='tight')
# plt.close()



hla_dic = {}
dic = {}
total_immuno = 0  # 1771
total_immuno_bio = 0
root_immuno_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/immunopeptidome'
with pd.ExcelWriter('all_immuno_meta.xlsx') as writer:
    for c in cancers:
        f = os.path.join(root_immuno_dir,cancers2immuno[c],'metadata.txt')
        df = pd.read_csv(f,sep='\t',index_col=0)
        total_immuno += df.shape[0]
        total_immuno_bio += len(df['biology'].unique())
        dic[c] = df.shape[0]
        df.to_excel(writer,sheet_name=c)
        hlas = []
        for item in df['HLA']:
            if isinstance(item,str):
                hlas.extend(list(set(item.split(','))))
        values,counts = np.unique(hlas,return_counts=True)
        hla_dic[c] = {v:c_ for v,c_ in zip(values,counts)}
print(total_immuno)
print(total_immuno_bio)
sys.exit('stop')

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
sys.exit('stop')









fig,ax = plt.subplots()
bars = ax.bar(x=np.arange(len(dic)),height=list(dic.values()))
for bar in bars:
    yval = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2, yval + 0.5, yval, ha='center', va='bottom')
ax.set_xticks(np.arange(len(dic)))
ax.set_xticklabels(list(dic.keys()),rotation=60)
ax.set_ylabel('Number of sample')
plt.savefig('figs1_stat_immuno.pdf',bbox_inches='tight')
plt.close()
    

# dic = {}
# total_rna = 0  # 7473
# with pd.ExcelWriter('all_rna_manifest.xlsx') as writer:
#     for c in cancers:
#         cmd = 'find {} -type f -name "manifest_*_*.tsv"'.format(os.path.join(root_atlas_dir,c))
#         possible_f = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1][0]
#         df = pd.read_csv(possible_f,sep='\t',index_col=0)
#         total_rna += df.shape[0]
#         dic[c] = df.shape[0]
#         df.to_excel(writer,sheet_name=c)

# fig,ax = plt.subplots()
# bars = ax.bar(x=np.arange(len(dic)),height=list(dic.values()))
# for bar in bars:
#     yval = bar.get_height()
#     ax.text(bar.get_x() + bar.get_width() / 2, yval + 0.5, yval, ha='center', va='bottom')
# ax.set_xticks(np.arange(len(dic)))
# ax.set_xticklabels(list(dic.keys()),rotation=60)
# ax.set_ylabel('Number of sample')
# plt.savefig('figs1_stat_rna.pdf',bbox_inches='tight')
# plt.close()


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












    

