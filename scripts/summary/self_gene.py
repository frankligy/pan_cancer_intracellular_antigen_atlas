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
    'SKCM',
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
bayests_xy_path = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/gene/full_results_XY_essential_tissues.txt'
gtex_median_path = '/gpfs/data/yarmarkovichlab/chordoma/NeoVerse_analysis/bulk-gex_v8_rna-seq_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct'
membrane_path = '/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/NeoVerse/data/human_membrane_proteins_acc2ens.txt'

def process_gtex(gtex_median_path):
    gtex = pd.read_csv(gtex_median_path,sep='\t',skiprows=2,index_col=0)
    gtex = gtex.loc[~gtex.index.str.contains('_PAR'),:]
    gtex.index = [item.split('.')[0] for item in gtex.index]
    gtex.drop(labels=['Cells - Cultured fibroblasts', 'Cells - EBV-transformed lymphocytes'],axis=1,inplace=True)
    gtex.drop(labels=['Ovary','Prostate','Testis','Vagina','Adrenal Gland','Cervix - Endocervix','Cervix - Ectocervix','Fallopian Tube','Pituitary'],axis=1,inplace=True)
    ensg2medians = {}
    for ensg in gtex.index:
        medians = gtex.loc[ensg,:].iloc[1:].values.tolist()
        ensg2medians[ensg] = medians
    all_tissues = gtex.columns[1:].values.tolist()
    ensg2symbol = pd.Series(index=gtex.index.tolist(),data=gtex['Description'].tolist()).to_dict()
    return ensg2medians,all_tissues,ensg2symbol

def process_tumor_gene():
    dic = {}
    for c in cancers:
        gene_lfc_path = os.path.join(root_atlas_dir,c,'gene_lfc.txt')
        gene_lfc = pd.read_csv(gene_lfc_path,sep='\t',index_col=0)
        c_dic = gene_lfc['median_tumor'].to_dict()
        dic[c] = c_dic
    return dic
        

# so it really does not change that much, this list versus the 66 that I used for enrichment
pan_cancer_ensgs = pd.read_csv('pan_cancer_cluster.txt',sep='\t',header=None)[0].values.tolist()

# pan-cancer cell cycle plot, let's just use the 68 that I selected and saved in fig2_raw
result = pd.read_csv('Reactome_Pathways_2024_table.txt',sep='\t').iloc[:6,:]
fig,ax = plt.subplots()
ax.barh(y=np.arange(result.shape[0]),width=np.flip(np.negative(np.log10(result['Adjusted P-value'].values))))
ax.set_xlabel('-log10(adjusted p-value)')
ax.set_yticks(np.arange(result.shape[0]))
ax.set_yticklabels(np.flip(result['Term'].values),fontsize=4)
plt.savefig('pan_cancer_cell_cycle_enrichr.pdf',bbox_inches='tight')
plt.close()

# for ts gene, coverage for patients, also consider all antigens
df = pd.read_csv('./stats/final_all_ts_antigens.txt',sep='\t')
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
df = df.loc[df['typ']=='self_gene',:]

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

c2c = plot_df_now.set_index(keys='cancer')['value'].to_dict()
c2c = {k:round(v,2) for k,v in c2c.items()}
custom_order = []
for i in plot_df_now['cancer']:
    if i not in custom_order:
        custom_order.append(i)
plot_df['cancer'] = pd.Categorical(plot_df['cancer'], categories=custom_order, ordered=True)
plot_df.sort_values(by=['cancer','category'],inplace=True)
plot_df.to_csv('self_gene_coverage.txt',sep='\t')
fig,ax = plt.subplots()
sns.barplot(plot_df,x='cancer',y='value',hue='category',ax=ax)
ax.set_xticklabels(['{}({})'.format(item.get_text(),str(c2c[item.get_text()])) for item in ax.get_xticklabels()], rotation=60)
plt.savefig('self_gene_coverage.pdf',bbox_inches='tight')
plt.close()


# get comparments, mainly for membrane protein
ensg2medians,all_tissues,ensg2symbol = process_gtex(gtex_median_path)
symbol2ensg = {v:k for k,v in ensg2symbol.items()}
df = pd.read_csv('uniprotkb_proteome_UP000005640_AND_revi_2024_12_26.tsv',sep='\t')
col = []
for i1,i2 in zip(df['Entry Name'],df['Gene Names']):
    gs_lis = []
    gs_lis.append(i1.split('_')[0])
    if isinstance(i2,str):
        for item in i2.split(' '):
            gs_lis.append(item)
    gs_lis = list(set(gs_lis))
    for gs in gs_lis:
        ensg = symbol2ensg.get(gs,None)
        if not ensg is None:
            break
    col.append(ensg)
df['ensg'] = col

keyword_map = {
    'Cell membrane':'cell_membrane',
    'Secreted':'secreted',
    'Lysosome':'lysosome',
    'Endosome':'endosome',
    'Golgi apparatus': 'golgi',
    'Endoplasmic reticulum':'er',
    'Nucleus':'nucleus',
    'Cytoplasm':'cytoplasma',
    'Mitochondrion': 'mitochondria',
    'Peroxisome':'peroxisome' 
}


for k,v in keyword_map.items():
    col = []
    for item in df['Subcellular location [CC]']:
        if isinstance(item,str) and k in item:  
            col.append(True)
        else:
            col.append(False)
    now_df = df.copy()
    now_df['is_cell_membrane'] = col
    now_df = now_df.loc[:,['ensg','Entry Name','Entry','is_cell_membrane','Subcellular location [CC]']]

    now_df = now_df.loc[now_df['is_cell_membrane']]
    now_df['gs'] = [item.split('_')[0] for item in now_df['Entry Name']]
    now_df.to_csv('./compartment/human_{}_protein_postdoc_final.txt'.format(v),sep='\t',index=None)


all_cancers = ['PAAD','OV','GBM','COAD','LUAD','LUSC','STAD','LIHC','BLCA','BRCA','KIRC','SKCM','HNSC','ESCA']
var_list = []
for c in all_cancers:
    var = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/logic_finder_v2/output_{}/var.csv'.format(c),sep=',',index_col=0)
    var_list.append(var)
d = pd.concat(var_list,axis=1,join='outer')
d.columns = all_cancers
d.to_csv('zach_single_cell_data.txt',sep='\t')
mapping = pd.read_csv('gProfiler_hsapiens_zach_sc.csv',sep=',',index_col=0)['converted_alias'].to_dict()
col = []
for item in d.index:
    col.append(mapping.get(item,None))
d.index = col
d.to_csv('zach_single_cell_data_ensg.txt',sep='\t')

# construct ensg2adjp
de_list = []
for c_ in cancers:
    de = pd.read_csv('./stats/DE_result_{}.txt'.format(c_),sep='\t',index_col=0)
    de_list.append(de['adjp'])
de_df = pd.concat(de_list,axis=1,join='outer')
de_df.columns = cancers
s = de_df.min(axis=1)
ensg2adjp = pd.Series(data=[True if i < 1e-5 else False for i in s],index=s.index)

# membrane protein coverage
mem = pd.read_csv('my_filter_membrane.txt',sep='\t')
membrane_ensg = list(set(mem['0'].values.tolist()))

mem_bl_ensg = [
    'ENSG00000144681',
    'ENSG00000050438',
    'ENSG00000179639',
    'ENSG00000187867',
    'ENSG00000146013',
    'ENSG00000160856',
    'ENSG00000143297',
    'ENSG00000169258',
    'ENSG00000182866',
    'ENSG00000158315',
    'ENSG00000139193',
    'ENSG00000188389',
    'ENSG00000134061',
    'ENSG00000101082',
    'ENSG00000196358',
    'ENSG00000126353',
    'ENSG00000094755',
    'ENSG00000227191',
    'ENSG00000152939',
    'ENSG00000156738',
    'ENSG00000116824',
    'ENSG00000167083',
    'ENSG00000253313',
    'ENSG00000196209',
    'ENSG00000161682',
    'ENSG00000164175'
]

mem_bl_no_extra_ensg = [
    'ENSG00000255587',
    'ENSG00000104537',
    'ENSG00000135605',
    'ENSG00000168421',
    'ENSG00000276231',
    'ENSG00000248905',
    'ENSG00000165304',
    'ENSG00000141293',
    'ENSG00000185686'
]

other_car = {
    'ENSG00000113361':'https://www.dl.begellhouse.com/download/article/2975dc725ae0753a/JEP(T)-40339.pdf',
    'ENSG00000134258':'https://pubmed.ncbi.nlm.nih.gov/27439899/',
    'ENSG00000102524':'https://pubmed.ncbi.nlm.nih.gov/35017485',
    'ENSG00000105369':'https://ashpublications.org/blood/article/140/Supplement%201/12716/493015/CRC-403-A-Phase-1-2-Study-of-bbT369-a-Dual-CD79a',
    'ENSG00000086548':'https://aacrjournals.org/cancerimmunolres/article/5/3_Supplement/A74/468728/Abstract-A74-CAR-T-cells-harboring-camelid-single',
    'ENSG00000066294':'https://ashpublications.org/blood/article/140/Supplement%201/7379/487194/CD84-A-Novel-Target-for-CAR-T-Cell-Therapy-for',
    'ENSG00000154269':'https://patentscope.wipo.int/search/es/detail.jsf;jsessionid=CCC5B7359899ECB0F03B7125D045EDF2.wapp2nC?docId=WO2024226829&_gid=202444',
    'ENSG00000172061':'https://pmc.ncbi.nlm.nih.gov/articles/PMC9604383/',
    'ENSG00000148848':'https://patents.google.com/patent/CA3133633A1/en',
    'ENSG00000114638':'https://patents.google.com/patent/WO2019232503A1/en',
    'ENSG00000137101':'https://ashpublications.org/blood/article/140/Supplement%201/7394/492100/Humanized-Nanobody-Anti-CD72-CAR-T-Cells',
    'ENSG00000137101':'https://pubmed.ncbi.nlm.nih.gov/36739093/',

}


# just bystander t cell
artifact_car = [
    'ENSG00000167286',
    'ENSG00000211753',
    'ENSG00000089692',
    'ENSG00000211697',
    'ENSG00000181847',
    'ENSG00000137078',
    'ENSG00000211689',
    'ENSG00000142484',
    'ENSG00000153283',
    'ENSG00000134460',
    'ENSG00000103522'
]

known_car = pd.read_csv('cart_targets.txt',sep='\t')
known_car = known_car.loc[known_car['Category'] == 'in clinical trials',:]
known_car_ensg = known_car['Ensembl ID'].values.tolist()
known_car_ensg = known_car_ensg + ['ENSG00000184697','ENSG00000213420','ENSG00000079112']  # missed CLDN6, GPC2, CDH17 from the 71 targets list

membrane_ensg = list(set(membrane_ensg).difference(set(mem_bl_ensg)))
membrane_ensg = list(set(membrane_ensg).difference(set(mem_bl_no_extra_ensg)))
membrane_ensg = list(set(membrane_ensg).difference(set(artifact_car)))

ensg2medians,all_tissues,ensg2symbol = process_gtex(gtex_median_path)
symbol2ensg = {v:k for k,v in ensg2symbol.items()}
all_genes = list(set(membrane_ensg).intersection(set(ensg2medians.keys())))
dic = process_tumor_gene()
ensg2tumors = {}
for gene in all_genes:
    data = []
    for k,v in dic.items():
        data.append(v.get(gene,[0]*21))
    ensg2tumors[gene] = data

all_data = []
for gene in all_genes:
    data = []
    data.extend(ensg2tumors[gene])
    data.extend(ensg2medians[gene])
    all_data.append(data)
df = pd.DataFrame.from_records(all_data,columns=cancers+all_tissues,index=all_genes)


ori_array = [tuple(['cancer']*21+['normal']*43),tuple(df.columns.tolist())]
mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
df.columns = mi

cond_col = []
for item in df.index:
    if item in ['ENSG00000048462','ENSG00000177455']:
        cond_col.append('1_FDR_approved')
    elif item in known_car_ensg:
        cond_col.append('2_clinical_trial')
    elif item in other_car.keys():
        cond_col.append('3_preclinical_test')
    else:
        cond_col.append('4_other_viable')

ori_array = [tuple(df.index.tolist()),
             tuple([ensg2symbol[item] for item in df.index]),
             tuple([ensg2adjp[item] for item in df.index]),
             tuple(cond_col)]
mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
df.index = mi
df.to_csv('gene_morpheus_mem.txt',sep='\t')
print(df.shape[0])


# intracellular
ts_final = pd.read_csv('./stats/final_all_ts_antigens.txt',sep='\t')
ts_final = ts_final.loc[ts_final['typ'] == 'self_gene',:]

ensg2dep = {}
for item1,item2 in zip(ts_final['ensgs'],ts_final['depmap_median']):
    ensg2dep[item1] = item2

ts_final.to_csv('ts_final_final.txt',sep='\t',index=None)
all_genes = list(set(ts_final['ensgs'].values))

# most common 
with pd.ExcelWriter('self_gene_common/self_peptide_common.xlsx') as writer:
    for cancer_,sub_df in ts_final.groupby(by='cancer'):
        sub_df.sort_values(by='recurrence',ascending=False).to_excel(writer,sheet_name='{}_recurrence'.format(cancer_),index=None)
    vc = ts_final['pep'].value_counts()
    lis = []
for p in vc.index:
    lis.append(ts_final.loc[ts_final['pep']==p,:])
s = pd.concat(lis,axis=0)
s.to_csv('self_gene_common/common_by_cancer.txt',sep='\t')

mem = pd.read_csv('human_membrane_proteins_acc2ens.txt',sep='\t')
mem = mem.loc[mem['Ens'].notna(),:]
membrane_ensg = mem['Ens'].values.tolist()
mem = pd.read_csv('./compartment/human_cell_membrane_protein_postdoc_final.txt',sep='\t')
mem = mem.loc[mem['Subcellular location [CC]'].str.contains('anchor'),:]
mem = mem.loc[mem['ensg'].notna(),:]
membrane_ensg_add = mem['ensg'].values.tolist()
membrane_ensg = set(membrane_ensg + membrane_ensg_add)

sc = pd.read_csv('zach_single_cell_data_ensg.txt',sep='\t',index_col=0)
sc = sc.loc[[True if isinstance(item,str) and item.startswith('ENSG') else False for item in sc.index],:]
s = sc.max(axis=1).to_frame()
s.columns = ['value']
s.sort_values(by='value',ascending=False,inplace=True)
s['percentile'] = 1 - np.arange(len(s)) / len(s)
s_map = s['percentile'].to_dict()

dic = process_tumor_gene()
ensg2tumors = {}
for gene in all_genes:
    data = []
    for k,v in dic.items():
        data.append(v[gene])
    ensg2tumors[gene] = data

all_data = []
for gene in all_genes:
    data = []
    data.extend(ensg2tumors[gene])
    data.extend(ensg2medians[gene])
    all_data.append(data)
df = pd.DataFrame.from_records(all_data,columns=cancers+all_tissues,index=all_genes)


ori_array = [tuple(['cancer']*21+['normal']*43),tuple(df.columns.tolist())]
mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
df.columns = mi

ori_array = [tuple(df.index.tolist()),
             tuple([ensg2symbol[item] for item in df.index]),
             tuple([True if item in membrane_ensg else False for item in df.index]),
             tuple([ensg2dep[item] for item in df.index]),
             tuple([ensg2adjp[item] for item in df.index]),
             tuple([s_map.get(item,None) for item in df.index])]
mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
df.index = mi
df.to_csv('gene_morpheus.txt',sep='\t')
print(df.shape[0])
print(len(ts_final['pep'].unique()))
    
    





    



    