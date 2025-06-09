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
    
data = []
for c in cancers:
    final_path = os.path.join(root_atlas_dir,c,'antigen','fdr','final_enhanced.txt')
    final = pd.read_csv(final_path,sep='\t')
    cond = [False if ('[]' in item) and ('(\'HLA-' not in item) else True for item in final['presented_by_each_sample_hla']]
    final = final.loc[cond,:]
    final = final.loc[final['typ']=='pathogen',:]
    data.append(final)
final = pd.concat(data,axis=0,keys=cancers).reset_index(level=-2).rename(columns={'level_0':'cancer'})

col = []
for item in final['source']:
    if '_HCMV' in item:
        col.append('CMV')
    elif '_FUSNU' in item:
        col.append('F.Nucleatum')
    elif '_EBV' in item:
        col.append('EBV')
    elif '_HPV' in item:
        col.append('HPV')
    elif '_HBV' in item:
        col.append('HBV')
    elif '_HELPJ' in item:
        col.append('H.Pylori')
    elif '_NIACI' in item:
        col.append('N.Circulans')
    elif '_9CLOT' in item:
        col.append('C.Intestinale')
    elif '_9BACT' in item:
        col.append('C.Ureolyticus')
    else:
        col.append('unknown')
final['strain'] = col
final = final.loc[final['strain']!='unknown',:]
final.to_csv('all_pathogen.txt',sep='\t',index=None)

sys.exit('stop')

safety_screen_df = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/code/post_safety_screen.txt',sep='\t')
safety_screen_df = safety_screen_df.loc[~safety_screen_df['cond_stringent'],:]
safety_screen_bl = list(set(safety_screen_df['pep'].values.tolist()))

final = final.loc[~final['pep'].isin(safety_screen_bl),:]
final.to_csv('all_pathogen_ts.txt',sep='\t',index=None)


# select representative
final_sub = final.loc[final['strain']=='F.Nucleatum',:].sort_values(by='n_psm',ascending=False)
fuso = final_sub.iloc[:5,:]['pep'].values.tolist()
fuso.extend(['GLDLGTLTY','TSAIGGAIAAV'])

final_sub = final.loc[final['strain']=='N.Circulans',:].sort_values(by='n_psm',ascending=False)
circ = final_sub.iloc[:5,:]['pep'].values.tolist()

final_sub = final.loc[final['strain']=='C.Intestinale',:].sort_values(by='n_psm',ascending=False)
intes = final_sub.iloc[:5,:]['pep'].values.tolist()

final_sub = final.loc[final['strain']=='C.Ureolyticus',:].sort_values(by='n_psm',ascending=False)
ureo = final_sub.iloc[:5,:]['pep'].values.tolist()

pylori = ['LSAKQNLEI','VRGKDIFII']

common_cmv = {
    'LLDGVTVSL':False,
    'TLFDEPPPL':False,
    'LDFGDLLKY':False,
    'LPVESLPLL':False,
    'RLQPNVPLV':True,
}

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
    'FPHCLAFSY':False,

}

common_hpv = {
    'YMLDLQPET':True
}

synthesized_bac = [
    'TIDELQKI',
    'LSDLGSGIYR',
    'IEKEVISKY',
    'QPKTKLLLL',
    'LDIHTFGLYY',
    'KLKPGILKK',
    'KMAEVIGLSK',
    'LVGPNGVGK'
]

ebv_peptides = ['LVEMALAGQR','GLKDAVYFL']
hpv_peptides = ['LESRITVFE']




raw_total = list(common_cmv.keys()) + list(common_hbv.keys()) + list(common_hpv.keys()) + synthesized_bac + fuso + circ + intes + ureo + pylori + ebv_peptides + hpv_peptides 
raw_total = list(set(raw_total))

final_now = final.loc[final['pep'].isin(raw_total),:]
candidates = list(set(final_now['pep'].values.tolist()))
final = final_now


# start to draw
pep2patho = {p:sub_df['strain'].iloc[0] for p,sub_df in final_now.groupby(by='pep')}
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

# # look for N circulans gene in ov
# final_n = final.loc[(final['strain']=='N.Circulans') & (final['cancer']=='OV'),:]
# genes = list(set([item.split('|')[1] for item in final_n['source']]))
# print(len(genes))
# sys.exit('stop')


# # look for cmv
# common_cmv = ['DLLSALQQL',
#               'TLLVYLFSL',
#               'VLEETSVML',
#               'IARLAKIPL',
#               'LLDGVTVSL',
#               'LPVESLPLL',
#               'YTSRGALYLY']
# final_cmv = final.loc[(final['strain']=='CMV') & (final['pep'].isin(common_cmv)),:]
# final_cmv.to_csv('final_cmv.txt',sep='\t')
# cmv_cancers = ['OV','BRCA','LIHC','NBL','GBM','CESC','COAD']

# store_data = []
# for pep in common_cmv:
#     final_p = final_cmv.loc[final_cmv['pep']==pep,:]
#     all_occur = final_p['cancer'].values.tolist()
#     tmp = []
#     for c in cmv_cancers:
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

# cmv_df = pd.DataFrame(data=store_data,index=common_cmv,columns=cmv_cancers)
# cmv_df.to_csv('cmv_df.txt',sep='\t')





        


