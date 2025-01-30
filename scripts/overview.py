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

manual_bl = [
    'FCGR1A',
    'TDGF1',
    'EPS8L3',
    'PCSK9',
    'HOXA9',
    'MUC13',
    'CTSW',
    'ADGRE2',
    'IRF4',
    'PRR15L',
    'MYEOV',
    'LYPD6B',
    'ATP6V1B1',
    'LGALS4',
    'GPA33',
    'DUOXA2',
    'FOXA3',
    'FOXA2',
    'SLC44A4',
    'SLAMF7',
    'ECEL1',
    'DPEP2',
    'UGT2A3',
    'PPP1R14D',
    'ATP2C2',
    'CBLC',
    'RASEF',
    'REG4',
    'OLIG2',
    'OVOL2',
    'B3GNT3',
    'TENM2',
    'CCNF',
    'PLEK2',
    'LIX1',
    'ESRP1',
    'ETV7',
    'CA14',
    'ITGB6',
    'NTS',
    'PRR15',
    'IGSF9',
    'MAP3K7CL',
    'LRRTM1',
    'TJP3',
    'LIPH',
    'HKDC1',
    'SCN9A',
    'GUCY2C',
    'TP63',
    'CXCL17',
    'CDX2',
    'LGR5',
    'CNTNAP2',
    'PLS3',
    'GPR87',
    'LARGE2',
    'ELF3',
    'KLHL14',
    'HLA-DQB2',
    'SLC28A1',
    'TRAF3IP3',
    'SLC17A4',
    'AP1M2',
    'KCTD14',
    'HLA-G',
    'GCNT1',
    'KLK8',
    'SLC6A15',
    'CXCL8',
    'IHH',
    'RNF175',
    'MCTP2',
    'CELF5',
    'INHBA',
    'TMPRSS4',
    'NLRP3',
    'GFRA2',
    'MAB21L1',
    'CEMIP',
    'RNASE2',
    'SNX20',
    'NKAIN4',
    'ARHGAP36',
    'ESPN',
    'ATP10B',
    'RNASE7',
    'GCSAM',
    'POU3F2',
    'CTSE',
    'SEZ6L',
    'IDO1',
    'CELSR1',
    'WNT5A',
    'ANKS4B',
    'SCN3A',
    'TMEM200A',
    'TMEM45B',
    'SUCNR1',
    'MISP',
    'TMEM150B',
    'HGF',
    'SATB2',
    'SLC5A10',
    'MIA',
    'CDX1',
    'BLK',
    'MSR1',
    'ELAVL2',
    'FAM129C',
    'TMEM139',
    'CLRN3',
    'KCNG1',
    'PMAIP1',
    'TMEM178B',
    'BIRC7',
    'KCNK15',
    'FBLL1',
    'C15orf48',
    'C9orf152',
    'C16orf54',
    'CSF2RA',
    'PLEKHG6'
]



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
    mem = pd.read_csv('human_membrane_protein_postdoc_final_no_edit.txt',sep='\t')
    mem = mem.loc[mem['ensg'].notna(),:]
    membrane_ensg = mem['ensg'].values.tolist()
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
        variants = variants.loc[(variants['ensg'].str.contains(r'^ENSG')) &                  
                                (~variants['mutation'].str.contains(r'^HLA')) &
                                (~variants['mutation'].str.contains(r'^IG')) &
                                (~variants['mutation'].str.contains(r'^TRAV')) & 
                                (~variants['mutation'].str.contains(r'^TRAJ')) &
                                (~variants['mutation'].str.contains(r'^TRBV')) & 
                                (~variants['mutation'].str.contains(r'^TRBD')) & 
                                (~variants['mutation'].str.contains(r'^TRBJ')) &
                                ((variants['n_samples']>1) | (variants['is_driver'])), :]   

        n_variant = variants.shape[0]

        
        
    return n_variant

def get_ts_pathogen(c):
    fixed_dict = {
        'BRCA':1,
        'KIRC':0,
        'COAD':2,
        'STAD':4,
        'MESO':0,
        'LIHC':2,
        'ESCA':3,
        'CESC':2,
        'BLCA':0,
        'RT':0,
        'AML':0,
        'DLBC':0,
        'GBM':1,
        'NBL':1,
        'PAAD':0,
        'HNSC':1,
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

safety_screen_df = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/code/post_safety_screen.txt',sep='\t')
safety_screen_df = safety_screen_df.loc[~safety_screen_df['cond'],:]
safety_screen_bl = list(set(safety_screen_df['pep'].values.tolist()))

# collage the meta and get number
total_antigen = 0  # 16077
self_df = pd.read_csv('ts_final.txt',sep='\t')
self_df = self_df.loc[self_df['highest_abundance'].notna(),:]
self_df = self_df.loc[~self_df['pep'].isin(safety_screen_bl),:]
total_antigen += len(set(self_df['pep'].values.tolist()))

self_translate_te_df = pd.read_csv('ts_te_antigen.txt',sep='\t')
self_translate_te_df = self_translate_te_df.loc[self_translate_te_df['highest_abundance'].notna(),:]
self_translate_te_df = self_translate_te_df.loc[~self_translate_te_df['pep'].isin(safety_screen_bl),:]
total_antigen += len(set(self_translate_te_df['pep'].values.tolist()))

te_chimeric_df = pd.read_csv('te_all_antigens.txt',sep='\t')
te_chimeric_df = te_chimeric_df.loc[te_chimeric_df['typ']=='TE_chimeric_transcript',:]
te_chimeric_df = te_chimeric_df.loc[te_chimeric_df['highest_abundance'].notna(),:]
te_chimeric_df = te_chimeric_df.loc[~te_chimeric_df['pep'].isin(safety_screen_bl),:]
total_antigen += len(set(te_chimeric_df['pep'].values.tolist()))

splicing_df = pd.read_csv('all_splicing.txt',sep='\t')
splicing_df = splicing_df.loc[splicing_df['highest_abundance'].notna(),:]
splicing_df = splicing_df.loc[~splicing_df['pep'].isin(safety_screen_bl),:]
total_antigen += len(set(splicing_df['pep'].values.tolist()))

nuorf_df = pd.read_csv('all_nuorf.txt',sep='\t')
nuorf_df = nuorf_df.loc[nuorf_df['highest_abundance'].notna(),:]
nuorf_df = nuorf_df.loc[nuorf_df['highest_score']>40,:]
nuorf_df = nuorf_df.loc[~nuorf_df['pep'].isin(safety_screen_bl),:]
total_antigen += len(set(nuorf_df['pep'].values.tolist()))

variant_df = pd.read_csv('all_variants.txt',sep='\t')
variant_df = variant_df.loc[variant_df['highest_abundance'].notna(),:]
variant_df = variant_df.loc[~variant_df['pep'].isin(safety_screen_bl),:]
total_antigen += len(set(variant_df['pep'].values.tolist()))

fusion_df = pd.read_csv('all_fusion.txt',sep='\t')
fusion_df = fusion_df.loc[fusion_df['highest_abundance'].notna(),:]
fusion_df = fusion_df.loc[~fusion_df['source'].str.contains('nc'),:]
fusion_df = fusion_df.loc[~fusion_df['pep'].isin(safety_screen_bl),:]
total_antigen += len(set(fusion_df['pep'].values.tolist()))

ir_df = pd.read_csv('all_ir.txt',sep='\t')
ir_df = ir_df.loc[ir_df['highest_abundance'].notna(),:]
ir_df = ir_df.loc[~ir_df['pep'].isin(safety_screen_bl),:]
total_antigen += len(set(ir_df['pep'].values.tolist()))

pathogen_df = pd.read_csv('all_pathogen.txt',sep='\t')
pathogen_df = pathogen_df.loc[pathogen_df['highest_abundance'].notna(),:]
pathogen_df = pathogen_df.loc[pathogen_df['unique'],:]
pathogen_df = pathogen_df.loc[pathogen_df['highest_score']>40,:]
pathogen_df = pathogen_df.loc[~pathogen_df['pep'].isin(safety_screen_bl),:]
total_antigen += len(set(pathogen_df['pep'].values.tolist()))

patent_df = pd.concat([self_df,self_translate_te_df,te_chimeric_df,splicing_df,nuorf_df,variant_df,fusion_df,ir_df,pathogen_df])
patent_df.to_csv('submission_all_antigens.txt',sep='\t',index=None)   # you can generate for safety screen as well
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

print(total_antigen)


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
# chimera_TE
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

sys.exit('stop')

hla_dic = {}
dic = {}
total_immuno = 0  # 1564
root_immuno_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/immunopeptidome'
with pd.ExcelWriter('all_immuno_meta.xlsx') as writer:
    for c in cancers:
        f = os.path.join(root_immuno_dir,cancers2immuno[c],'metadata.txt')
        df = pd.read_csv(f,sep='\t',index_col=0)
        total_immuno += df.shape[0]
        dic[c] = df.shape[0]
        df.to_excel(writer,sheet_name=c)
        hlas = []
        for item in df['HLA']:
            if isinstance(item,str):
                hlas.extend(list(set(item.split(','))))
        values,counts = np.unique(hlas,return_counts=True)
        hla_dic[c] = {v:c_ for v,c_ in zip(values,counts)}

hla_data = []
all_hla = []
for k,v in hla_dic.items():
    all_hla.extend(list(v.keys()))
all_hla = list(set(all_hla))
for c in cancers:
    tmp = []
    for hla in all_hla:
        tmp.append(hla_dic[c].get(hla,0))
    hla_data.append(tmp)

hla_df = pd.DataFrame(index=cancers,columns=all_hla,data=hla_data)
freq_dic = pd.read_csv('/gpfs/data/yarmarkovichlab/medulloblastoma/neoverse_folder/NeoVerse_final_output_new/antigens/US_HLA_frequency.csv',sep=',',index_col=0)['Percent US population'].to_dict()

reformatted_hla = ['HLA-' + item.replace(':','').replace('*','') for item in hla_df.columns]
freq_df = pd.read_csv('/gpfs/data/yarmarkovichlab/medulloblastoma/neoverse_folder/NeoVerse_final_output_new/antigens/US_HLA_frequency.csv',sep=',',index_col=0)
# HLA-A
freq_a = freq_df.loc[[True if item.startswith('HLA-A') else False for item in freq_df.index],:]
freq_a = freq_a.loc[freq_a['Percent US population']>0.01,:]
have_a = [item for item in reformatted_hla if item.startswith('HLA-A')]
prop = len(set(freq_a.index).intersection(set(have_a))) / len(freq_a)
print(prop)
# HLA-B
freq_a = freq_df.loc[[True if item.startswith('HLA-B') else False for item in freq_df.index],:]
freq_a = freq_a.loc[freq_a['Percent US population']>0.01,:]
have_a = [item for item in reformatted_hla if item.startswith('HLA-B')]
prop = len(set(freq_a.index).intersection(set(have_a))) / len(freq_a)
print(prop)
# HLA-C
freq_a = freq_df.loc[[True if item.startswith('HLA-C') else False for item in freq_df.index],:]
freq_a = freq_a.loc[freq_a['Percent US population']>0.01,:]
have_a = [item for item in reformatted_hla if item.startswith('HLA-C')]
prop = len(set(freq_a.index).intersection(set(have_a))) / len(freq_a)
print(prop)
sys.exit('stop')

ori_array = [tuple([freq_dic.get('HLA-' + item.replace(':','').replace('*',''),0) for item in hla_df.columns]),tuple(hla_df.columns.tolist())]
mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
hla_df.columns = mi
hla_df.to_csv('hla_df.txt',sep='\t')




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
    

dic = {}
total_rna = 0  # 7473
with pd.ExcelWriter('all_rna_manifest.xlsx') as writer:
    for c in cancers:
        cmd = 'find {} -type f -name "manifest_*_*.tsv"'.format(os.path.join(root_atlas_dir,c))
        possible_f = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1][0]
        df = pd.read_csv(possible_f,sep='\t',index_col=0)
        total_rna += df.shape[0]
        dic[c] = df.shape[0]
        df.to_excel(writer,sheet_name=c)

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

sys.exit('stop')

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

print(n_ts_total,n_ts_membrane)

fig,ax = plt.subplots(figsize=(15,10))
sns.heatmap(plot_df,annot=True,linewidth=0.5,fmt='g',annot_kws={"fontsize": 5},cmap='Blues',square=True)
plt.savefig('ts_event_overview.pdf',bbox_inches='tight')
plt.close()












    

