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

prop = [
    0.76,
    0.99,
    1.0,
    0.56,
    0.70,
    0.9,
    0.90,
    0.97,
    0.0,
    0.16,
    1.0,
    0.99,
    0.98,
    1.0,
    0.99,
    0.99,
    1.0,
    0.99,
    1.0,
    0.64,
    1.0
]

e_t = 0
for i1,i2 in zip(prop,n_samples):
    e_t += i1*i2
print(e_t/sum(n_samples))


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
        

pan_cancer_ensgs = [
    'ENSG00000100162',
    'ENSG00000101057',
    'ENSG00000161888',
    'ENSG00000089685',
    'ENSG00000167513',
    'ENSG00000093009',
    'ENSG00000091651',
    'ENSG00000080986',
    'ENSG00000142945',
    'ENSG00000131153',
    'ENSG00000146670',
    'ENSG00000237649',
    'ENSG00000175643',
    'ENSG00000131351',
    'ENSG00000127564',
    'ENSG00000162062',
    'ENSG00000137807',
    'ENSG00000131747',
    'ENSG00000117724',
    'ENSG00000112984',
    'ENSG00000138160',
    'ENSG00000126787',
    'ENSG00000072571',
    'ENSG00000071539',
    'ENSG00000101003',
    'ENSG00000171848',
    'ENSG00000094804',
    'ENSG00000090889',
    'ENSG00000138180',
    'ENSG00000112742',
    'ENSG00000088325',
    'ENSG00000165304',
    'ENSG00000170312',
    'ENSG00000117650',
    'ENSG00000087586',
    'ENSG00000175063',
    'ENSG00000117399',
    'ENSG00000134690',
    'ENSG00000115163',
    'ENSG00000157456',
    'ENSG00000100526',
    'ENSG00000164087',
    'ENSG00000105011',
    'ENSG00000169679',
    'ENSG00000148773',
    'ENSG00000143228',
    'ENSG00000186185',
    'ENSG00000135451',
    'ENSG00000111206',
    'ENSG00000151725',
    'ENSG00000143476',
    'ENSG00000156970',
    'ENSG00000111247',
    'ENSG00000183856',
    'ENSG00000179750',
    'ENSG00000189057',
    'ENSG00000123219',
    'ENSG00000276043',
    'ENSG00000119969',
    'ENSG00000166803',
    'ENSG00000169245',
    'ENSG00000155893',
    'ENSG00000129195',
    'ENSG00000168078',
    'ENSG00000121211',
    'ENSG00000144354',
    'ENSG00000162390',
    'ENSG00000196083',
    'ENSG00000160886',
    'ENSG00000130487',
    'ENSG00000117407',
    'ENSG00000176597',
]

# pan-cancer cell cycle plot, let's just use the 68 that I selected and saved in fig2_raw
result = pd.read_csv('Reactome_Pathways_2024_table.txt',sep='\t').iloc[:6,:]
fig,ax = plt.subplots()
ax.barh(y=np.arange(result.shape[0]),width=np.flip(np.negative(np.log10(result['Adjusted P-value'].values))))
ax.set_xlabel('-log10(adjusted p-value)')
ax.set_yticks(np.arange(result.shape[0]))
ax.set_yticklabels(np.flip(result['Term'].values),fontsize=4)
plt.savefig('pan_cancer_cell_cycle_enrichr.pdf',bbox_inches='tight')
plt.close()
sys.exit('stop')

# for ts gene, coverage for patients
df = pd.read_csv('ts_final.txt',sep='\t')
ts_ensg = df['ensgs'].unique().tolist()
ts_ensg = list(set(ts_ensg).difference(set(pan_cancer_ensgs)))
data = []
for c in cancers:
    print(c)
    path = os.path.join(root_atlas_dir,c,'gene_tpm.txt')
    exp = pd.read_csv(path,sep='\t',index_col=0)
    exp = exp.loc[ts_ensg,:]
    exp = exp.values
    exp = (exp > 20).astype(int)
    props = np.sum(exp,axis=1) / exp.shape[1]
    mapping = {} # gene specific prop
    for i1,i2 in zip(ts_ensg,props):
        mapping[i1] = i2
    prop = np.count_nonzero(np.any(exp,axis=0)) / exp.shape[1]
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
    df_sub = df.loc[(df['cancer']==c) & (df['ensgs'].isin(ts_ensg)),:]
    mapping_norm_sb = {}
    mapping_norm_wb = {}
    for ensg,sub_df in tqdm(df_sub.groupby(by='ensgs')):
        g_prop = mapping[ensg]
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
        mapping_norm_sb[ensg] = norm_prop

        tmp = [dic.get(h.replace('HLA-','').replace('*','').replace(':',''),0) for h in wb_hla]
        h_prop = 1
        for item in tmp:
            h_prop *= (1-item)
        norm_prop = g_prop * (1-h_prop)
        mapping_norm_wb[ensg] = norm_prop

    final_prop = 1
    for k,v in mapping_norm_sb.items():
        final_prop *= (1-v)
    final_sb_prop = 1-final_prop

    final_prop = 1
    for k,v in mapping_norm_wb.items():
        final_prop *= (1-v)
    final_wb_prop = 1-final_prop

    data.append((prop,'gene_prop',c))
    data.append((final_wb_prop,'wb_prop',c))
    data.append((final_sb_prop,'sb_prop',c))

plot_df = pd.DataFrame.from_records(data,columns=['value','category','cancer'])
rename = {
    'BRCA':'16-BRCA(0.76)',
    'KIRC':'10-KIRC(0.99)',
    'COAD':'01-COAD(1.0)',
    'STAD':'19-STAD(0.56)',
    'MESO':'17-MESO(0.70)',
    'LIHC':'08-LIHC(0.9)',
    'ESCA':'15-ESCA(0.90)',
    'CESC':'14-CESC(0.97)',
    'BLCA':'21-BLCA(0.0)',
    'RT':'20-RT(0.16)',
    'AML':'02-AML(1.0)',
    'DLBC':'12-DLBC(0.99)',
    'GBM':'13-GBM(0.98)',
    'NBL':'04-NBL(1.0)',
    'PAAD':'07-PAAD(0.99)',
    'HNSC':'11-HNSC(0.99)',
    'OV':'05-OV(1.0)',
    'LUSC':'08-LUSC(0.99)',
    'LUAD':'06-LUAD(1.0)',
    'CHOL':'18-CHOL(0.64)',
    'SKCM':'03-SKCM(1.0)'
}
plot_df['cancer'] = [rename[item] for item in plot_df['cancer']]
plot_df.sort_values(by=['cancer','category'],inplace=True)
plot_df.to_csv('self_gene_coverage.txt',sep='\t')
fig,ax = plt.subplots()
sns.barplot(plot_df,x='cancer',y='value',hue='category',ax=ax)
ax.set_xticklabels(ax.get_xticklabels(), rotation=60)
plt.savefig('self_gene_coverage.pdf',bbox_inches='tight')
plt.close()

fig,ax = plt.subplots()
plot_df = plot_df.loc[plot_df['category']=='sb_prop',:].sort_values(by='value',ascending=False)
bars = ax.bar(x=np.arange(plot_df.shape[0]),height=plot_df['value'])
for bar in bars:
    yval = round(bar.get_height(),3)
    ax.text(bar.get_x() + bar.get_width() / 2, yval + 0.5, yval, ha='center', va='bottom')
ax.set_xticks(np.arange(plot_df.shape[0]))
ax.set_xticklabels(plot_df['cancer'],rotation=60)
ax.set_ylabel('HLA normalized patient coverage')
plt.savefig('self_gene_coverage_sb.pdf',bbox_inches='tight')
plt.close()

sys.exit('stop')


        





data = []
for c in cancers:
    final_path = os.path.join(root_atlas_dir,c,'antigen','fdr','final_enhanced.txt')
    final = pd.read_csv(final_path,sep='\t')
    cond = [False if '[]' in item else True for item in final['presented_by_each_sample_hla']]
    final = final.loc[cond,:]
    final = final.loc[final['typ']=='self_gene',:]
    final = final.loc[final['unique']!=False,:]
    data.append(final)
final = pd.concat(data,axis=0,keys=cancers).reset_index(level=-2).rename(columns={'level_0':'cancer'})

final.to_csv('all_unique_self_gene.txt',sep='\t',index=None)
final['gene_symbol'].value_counts().to_csv('check_gene.txt',sep='\t')
final['pep'].value_counts().to_csv('check_pep.txt',sep='\t')

ensg2dep = {}
for item1,item2 in zip(final['ensgs'],final['depmap_median']):
    ensg2dep[item1] = item2

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



ensg2medians,all_tissues,ensg2symbol = process_gtex(gtex_median_path)
symbol2ensg = {v:k for k,v in ensg2symbol.items()}
manual_bl_ensg = [symbol2ensg[item] for item in manual_bl]
all_genes = list(set(final['ensgs'].values))
all_genes = list(set(all_genes).difference(set(manual_bl_ensg)))

ts_final = final.loc[~final['ensgs'].isin(manual_bl_ensg),:]
ts_final.to_csv('ts_final.txt',sep='\t',index=None)


dic = process_tumor_gene()
ensg2tumors = {}
for gene in all_genes:
    data = []
    for k,v in dic.items():
        data.append(v[gene])
    ensg2tumors[gene] = data


df = pd.read_csv(membrane_path,sep='\t')
membrane_ensg = set(df['Ens'].tolist())

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
             tuple([ensg2dep[item] for item in df.index])]
mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
df.index = mi
df.to_csv('gene_morpheus.txt',sep='\t')
    
    





    



    