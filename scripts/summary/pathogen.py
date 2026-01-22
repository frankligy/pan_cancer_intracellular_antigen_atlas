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
import anndata as ad
from matplotlib_venn import venn2,venn3

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
gtex_median_path = '/gpfs/data/yarmarkovichlab/chordoma/NeoVerse_analysis/bulk-gex_v8_rna-seq_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct'
GTEX_GENE_ALL_H5AD = '/gpfs/data/yarmarkovichlab/chordoma/NeoVerse_analysis/gtex_gene_all.h5ad'


def categorize_sets(set1, set2, set3):
    return [
        set1 - set2 - set3,  # Unique to set1
        set2 - set1 - set3,  # Unique to set2
        (set1 & set2) - set3,  # Overlap between set1 and set2 but not in set3
        set3 - set1 - set2,  # Unique to set3
        (set3 & set1) - set2,  # Overlap between set3 and set1 but not in set2
        (set3 & set2) - set1,  # Overlap between set3 and set2 but not in set1
        set1 & set2 & set3  # Overlap among all three
    ]


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

def plot_gtex_all(ensg):

    # process gtex
    gtex = pd.read_csv(gtex_median_path,sep='\t',skiprows=2,index_col=0)
    cond = ~gtex.columns.isin(['Cells - EBV-transformed lymphocytes','Cells - Cultured fibroblasts','Testis'])
    gtex = gtex.loc[:,cond]
    gtex.index = [item.split('.')[0] for item in gtex.index]
    ensg2symbol = pd.Series(index=gtex.index.tolist(),data=gtex['Description'].tolist()).to_dict()
    series = gtex.loc[ensg,:].iloc[1:]
    gs = ensg2symbol[ensg]

    # all
    adata = ad.read_h5ad(GTEX_GENE_ALL_H5AD)  # 56200 × 17382
    adata.obs_names = [item.split('.')[0] for item in adata.obs_names]
    adata.obs_names_make_unique()
    adata_gene = adata[[ensg],:]
    normal_expr_list = []
    for t in series.index.tolist():
        values = adata_gene[:,adata_gene.var['tissue']==t].X.toarray().reshape(-1)
        normal_expr_list.append(values)
    
    tumor_expr_list = []
    for c in cancers:
        gene_path = os.path.join(root_atlas_dir,c,'gene_tpm.txt')
        final = pd.read_csv(gene_path,sep='\t',index_col=0)
        tumor_expr = final.loc[ensg,:].values
        tumor_expr_list.append(tumor_expr)

    all_list = tumor_expr_list + normal_expr_list

    fig,ax = plt.subplots()
    bp = ax.boxplot(x=all_list,positions=np.arange(len(all_list)),patch_artist=True)
    for flier in bp['fliers']:
        flier.set_markersize(1)
        flier.set_marker('o')
    for box in bp['boxes']:
        box.set_facecolor('green')
        box.set_edgecolor('black')
        box.set_linewidth(1)
    ax.set_xticks(np.arange(len(all_list)))
    ax.set_xticklabels(cancers + series.index.tolist(),fontsize=5,rotation=90)

    fig.suptitle(gs)
    plt.savefig('apobec_{}_{}.pdf'.format(ensg,gs),bbox_inches='tight')
    plt.close()

def process_tumor_gene():
    dic = {}
    for c in cancers:
        gene_lfc_path = os.path.join(root_atlas_dir,c,'gene_lfc.txt')
        gene_lfc = pd.read_csv(gene_lfc_path,sep='\t',index_col=0)
        c_dic = gene_lfc['median_tumor'].to_dict()
        dic[c] = c_dic
    return dic

# # APOBEC3 and ADAR analysis
# genes = {
#     'APOBEC3A':'ENSG00000128383',
#     'APOBEC3B':'ENSG00000179750',
#     'APOBEC3C':'ENSG00000244509',
#     'APOBEC3D':'ENSG00000243811',
#     'APOBEC3F':'ENSG00000128394',
#     'APOBEC3G':'ENSG00000239713',
#     'APOBEC3H':'ENSG00000100298',
#     'APOBEC1':'ENSG00000111701',
#     'APOBEC2':'ENSG00000124701',
#     'APOBEC4':'ENSG00000173627',
#     'ADAR':'ENSG00000160710',
#     'AICDA':'ENSG00000111732'
# }

# for ensg in genes.values():
#     plot_gtex_all(ensg)


# ensg2medians,all_tissues,ensg2symbol = process_gtex(gtex_median_path)
# symbol2ensg = {v:k for k,v in ensg2symbol.items()}
# all_genes = [symbol2ensg[item] for item in genes.keys()]
# dic = process_tumor_gene()
# ensg2tumors = {}
# for gene in all_genes:
#     data = []
#     for k,v in dic.items():
#         data.append(v[gene])
#     ensg2tumors[gene] = data

# all_data = []
# for gene in all_genes:
#     data = []
#     data.extend(ensg2tumors[gene])
#     data.extend(ensg2medians[gene])
#     all_data.append(data)
# df = pd.DataFrame.from_records(all_data,columns=cancers+all_tissues,index=all_genes)

# ori_array = [tuple(['cancer']*21+['normal']*43),tuple(df.columns.tolist())]
# mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
# df.columns = mi

# ori_array = [tuple(df.index.tolist()),
#              tuple([ensg2symbol[item] for item in df.index])]
# mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
# df.index = mi
# df.to_csv('gene_apobec.txt',sep='\t')
    

df = pd.read_csv('./stats/final_all_ts_antigens.txt',sep='\t')
final = df.loc[df['typ']=='pathogen',:]
final.to_csv('all_pathogen_final.txt',sep='\t',index=None)

# hbv
common_hbv = {
    'FLLTRILTI':True,
    'TRILTIPQSL':False,
    'LTIPQSLDSW':False,
    'FVGLSPTVWL':False,
    'STLPETTVVRR':True,
    'GTLPQEHIVQK':False,
    'GLSPTVWLSV':True,
    'GSTHVSWPK':True,
    'LPSDFFPSI':False,
    'IPIPSSWAF':True,
    'GVWIRTPPAYR':False,
    'LTIPQSLDSWW':False,
    'FVGLSPTVW':False,
    'HLYSHPIIL':True,
    'YPALMPLYA':False,
    'KYTSFPWLL':True,
    'LPFRPTTGR':False,
    'ASRELVVSY':True,
    'ILSPFLPLL':True,
    'GMLPVCPLI':False,
    'GSLPQEHIIQK':False,
    'TASPLSSIF':False,
    'TCIPIPSSW':False,
    'VGLSPTVWL':True
}    # intersect hbv iedb and hbv immunoverse
ratio_hbv = len([k for k,v in common_hbv.items() if v])/len(common_hbv)
print(ratio_hbv)
iedb_hbv = pd.read_csv('hbv_mhc_i_no_b_host_human.tsv',sep='\t')
iedb_hbv = set([item for item in iedb_hbv['Epitope - Name'] if len(item) >=8 and len(item) <=15])
my_hbv = set(final.loc[final['strain']=='HBV','pep'].values.tolist())
print(len(iedb_hbv - my_hbv),len(iedb_hbv & my_hbv), len(my_hbv - iedb_hbv))
tmp = final.loc[final['strain']=='HBV',:]
tmp['documented_immunogenicity'] = [common_hbv.get(item,None) for item in tmp['pep']]
tmp.to_csv('annotated_hbv_table.txt',sep='\t',index=None)

# hpv, none

# cmv
common_cmv = {
    'LLDGVTVSL':False,
    'LDFGDLLKY':False,
    'LPVESLPLL':False,
    'RLQPNVPLV':True,
}    # generated by intersecpting cmv iedb and cmv immunoverse
my_cmv = set(final.loc[final['strain']=='CMV','pep'].values.tolist())
print(len(my_cmv))
iedb_cmv = pd.read_csv('cmv_mhc_i_no_b_host_human.tsv',sep='\t')
iedb_cmv = set([item for item in iedb_cmv['Epitope - Name'] if len(item) >=8 and len(item) <=15])
print(len(iedb_cmv - my_cmv),len(iedb_cmv & my_cmv), len(my_cmv - iedb_cmv))
tmp = final.loc[final['strain']=='CMV',:]
tmp['documented_immunogenicity'] = [common_cmv.get(item,None) for item in tmp['pep']]
tmp.to_csv('annotated_cmv_table.txt',sep='\t',index=None)



# select representative
candidates = []

final_sub = final.loc[final['strain']=='HBV',:]
for c_,final_sub2 in final_sub.groupby(by='cancer'):
    tmp = final_sub2.sort_values(by='n_psm')['pep'].iloc[-2:]
    candidates.extend(list(tmp))

final_sub = final.loc[final['strain']=='HPV',:]
for c_,final_sub2 in final_sub.groupby(by='cancer'):
    tmp = final_sub2.sort_values(by='n_psm')['pep'].iloc[-2:]
    candidates.extend(list(tmp))

final_sub = final.loc[final['strain']=='CMV',:]
for c_,final_sub2 in final_sub.groupby(by='cancer'):
    tmp = final_sub2.sort_values(by='n_psm')['pep'].iloc[-2:]
    candidates.extend(list(tmp))

final_sub = final.loc[final['strain']=='EBV',:]
for c_,final_sub2 in final_sub.groupby(by='cancer'):
    tmp = final_sub2.sort_values(by='n_psm')['pep'].iloc[-2:]
    candidates.extend(list(tmp))

final_sub = final.loc[final['strain']=='F.Nucleatum',:]
for c_,final_sub2 in final_sub.groupby(by='cancer'):
    tmp = final_sub2.sort_values(by='n_psm')['pep'].iloc[-2:]
    candidates.extend(list(tmp))

final_sub = final.loc[final['strain']=='H.Pylori',:]
for c_,final_sub2 in final_sub.groupby(by='cancer'):
    tmp = final_sub2.sort_values(by='n_psm')['pep'].iloc[-2:]
    candidates.extend(list(tmp))

final_sub = final.loc[final['strain']=='C.Intestinale',:]
for c_,final_sub2 in final_sub.groupby(by='cancer'):
    tmp = final_sub2.sort_values(by='n_psm')['pep'].iloc[-2:]
    candidates.extend(list(tmp))

final_sub = final.loc[final['strain']=='C.Ureolyticus',:]
for c_,final_sub2 in final_sub.groupby(by='cancer'):
    tmp = final_sub2.sort_values(by='n_psm')['pep'].iloc[-2:]
    candidates.extend(list(tmp))

final_sub = final.loc[final['strain']=='N.Circulans',:]
for c_,final_sub2 in final_sub.groupby(by='cancer'):
    tmp = final_sub2.sort_values(by='n_psm')['pep'].iloc[-2:]
    candidates.extend(list(tmp))

# start to draw
candidates = list(set(candidates))
pep2patho = {p:sub_df['strain'].iloc[0] for p,sub_df in final.groupby(by='pep')}
safety_screen_df = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/code/hla_ligand_atlas_now_0.05_tesorai.txt',sep='\t')

all_tissues = ['Adrenal gland', 'Aorta', 'Bladder', 'Bone marrow', 'Brain', 'Cerebellum', 'Colon', 'Esophagus', 'Gallbladder', 'Heart', 'Kidney', 'Liver', 
                'Lung', 'Lymph node', 'Mamma', 'Muscle', 'Myelon', 'Ovary', 'Pancreas', 'Prostate', 'Skin', 'Small intestine', 'Spleen', 'Stomach', 'Testis', 
                'Thymus', 'Thyroid', 'Tongue', 'Trachea', 'Uterus']

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
    store_type.append(pep2patho[pep])


df = pd.DataFrame(data=store_data,index=candidates,columns=cancers+list(all_tissues))

ori_array = [tuple(['cancer']*21+['normal']*30),tuple(df.columns.tolist())]
mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
df.columns = mi

ori_array = [tuple(df.index.tolist()),
             tuple(store_type)]
mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
df.index = mi

df.to_csv('peptide_view_pathogen.txt',sep='\t')
sys.exit('stop')

# n.circulans
tmp = final.loc[final['strain']=='N.Circulans',:]
tmp_ov = tmp.loc[tmp['cancer']=='OV',:]
print(len(tmp_ov['pep'].unique()))
unique_genes = set()
for item in tmp_ov['source']:
    unique_genes.add(item.split('|')[1])
print(len(unique_genes))

# assemble to tesorai nbl (only peptide, or whole)
nia_peptide = set(tmp['pep'].unique())
dic = {}
with open('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/immunopeptidome/neuroblastoma/combined_NBL_pan_cancer.fasta','r') as in_handle:
    for title,seq in SimpleFastaParser(in_handle):
        dic[title] = seq
for pep in tmp['pep']:
    dic['test|{}'.format(pep)] = pep
with open('nbl_neg_circulans_only_pep.fasta','w') as f:
    for k,v in dic.items():
        f.write('>{}\n{}\n'.format(k,v))

dic = {}
with open('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/immunopeptidome/neuroblastoma/combined_NBL_pan_cancer.fasta','r') as in_handle:
    for title,seq in SimpleFastaParser(in_handle):
        dic[title] = seq
with open('/gpfs/data/yarmarkovichlab/public/ImmunoVerse/search_space_tesorai/OV/db_fasta_tesorai/Niallia_circulans_UP000319837.fasta','r') as in_handle:
    for title,seq in SimpleFastaParser(in_handle):
        dic[title] = seq
with open('nbl_neg_circulans_whole_proteome.fasta','w') as f:
    for k,v in dic.items():
        f.write('>{}\n{}\n'.format(k,v))






















        


