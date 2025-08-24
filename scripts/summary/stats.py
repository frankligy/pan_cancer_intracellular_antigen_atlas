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
from scipy.sparse import csr_matrix
from tqdm import tqdm
from scipy.stats import mannwhitneyu,chi2
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
from statsmodels.othermod.betareg import BetaModel

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


def run_self_gene_de(ensgs,cancer):

    final = pd.read_csv(os.path.join(root_atlas_dir,cancer,'gene_tpm.txt'),sep='\t',index_col=0)
    final = final.loc[~final.index.duplicated(),:]
    tumor_expr = final.loc[ensgs,:].values
    tumor_label = ['tumor{}'.format(i+1) for i in range(tumor_expr.shape[1])]

    gtex = pd.read_csv(os.path.join(database_dir,'bulk-gex_v8_rna-seq_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct'),sep='\t',skiprows=2,index_col=0)
    cond = ~gtex.columns.isin(['Cells - EBV-transformed lymphocytes','Cells - Cultured fibroblasts','Testis'])
    gtex = gtex.loc[:,cond]
    gtex.index = [item.split('.')[0] for item in gtex.index]
    ensg2symbol = pd.Series(index=gtex.index.tolist(),data=gtex['Description'].tolist()).to_dict()
    series = gtex.loc[[ensgs[0]],:].iloc[0,:].iloc[1:]

    adata = ad.read_h5ad(os.path.join(database_dir,'gtex_gene_all.h5ad'))  # 56200 × 17382
    adata.obs_names = [item.split('.')[0] for item in adata.obs_names]
    adata.obs_names_make_unique()
    adata_gene = adata[ensgs,:]
    normal_expr_list = []
    normal_expr_list_label = []
    for t in series.index.tolist():
        values = adata_gene[:,adata_gene.var['tissue']==t].X.toarray()
        normal_label = ['{}{}'.format(t,i+1) for i in range(values.shape[1])]
        normal_expr_list.append(values)
        normal_expr_list_label.extend(normal_label)

    all_array = [tumor_expr] + normal_expr_list
    df_data = np.concatenate(all_array,axis=1)
    df = pd.DataFrame(data=df_data,columns=tumor_label+normal_expr_list_label,index=ensgs)
    df.to_csv('exp.original-steady-state.txt',sep='\t')
    with open('groups.txt','w') as f:
        for item in tumor_label:
            f.write('{}\t1\ttumor\n'.format(item))
        for item in normal_expr_list_label:
            f.write('{}\t2\tnormal\n'.format(item))
    with open('comps.txt','w') as f:
        f.write('1\t2\n')

    # run altanalyze, assuming singularity is enabled
    os.makedirs('./diff_dir/altanalyze_output/ExpressionInput')
    os.rename('exp.original-steady-state.txt','./diff_dir/altanalyze_output/ExpressionInput/exp.original-steady-state.txt')
    os.rename('groups.txt','./diff_dir/groups.txt')
    os.rename('comps.txt','./diff_dir/comps.txt')
    old_dir = os.getcwd()
    os.chdir('./diff_dir')
    writable_path = '/gpfs/data/yarmarkovichlab/Frank/SNAF_ecosystem/altanalyze'
    cmd = 'singularity run -B $PWD:/mnt --writable {} DE altanalyze_output groups.txt'.format(writable_path)
    subprocess.run(cmd,shell=True)
    os.chdir(old_dir)
    os.rename('./diff_dir/altanalyze_output/ExpressionInput/DEGs-LogFold_0.0_adjp/GE.tumor_vs_normal.txt','DE_result_{}.txt'.format(cancer))
    subprocess.run('rm -rf ./diff_dir',shell=True)

def run_te_de(tes,cancer,is_te):
    tumor_erv = ad.read_h5ad(os.path.join(root_atlas_dir,cancer,'tumor_erv.h5ad'))
    tumor_erv_aux = pd.read_csv(os.path.join(root_atlas_dir,cancer,'tumor_erv_aux_df.txt'),sep='\t',index_col=0)
    normal_erv = pd.read_csv(os.path.join(HG38_NORMAL_DIR,'normal_erv.txt'),sep='\t',index_col=0)
    normal_erv_aux = pd.read_csv(os.path.join(HG38_NORMAL_DIR,'normal_erv_aux_df.txt'),sep='\t',index_col=0)

    n_data = normal_erv.values / normal_erv_aux['total_count'].values.reshape(1,-1) * 1e6

    # reorg
    tumor_erv.X = tumor_erv.layers['cpm']
    del tumor_erv.layers['cpm']
    normal_erv = pd.DataFrame(data=n_data,index=normal_erv.index,columns=normal_erv.columns)

    # use simple name
    tumor_erv.obs_names = [index.split(':')[0] for index in tumor_erv.obs_names]
    normal_erv.index = [index.split(':')[0] for index in normal_erv.index]
    
    # extract tumor 
    common = list(set(tes).intersection(set(tumor_erv.obs_names)))
    tumor_expr = tumor_erv[common,:].X.toarray()
    tumor_df = pd.DataFrame(data=tumor_expr,index=common,columns=['tumor{}'.format(i+1) for i in range(tumor_expr.shape[1])])

    # extract normal
    common = list(set(tes).intersection(set(normal_erv.index)))
    normal_df = normal_erv.loc[common,:]

    # combine
    final_df = pd.concat([tumor_df,normal_df],axis=1,join='outer').fillna(value=0)

    # LIMMA
    final_df.to_csv('exp.original-steady-state.txt',sep='\t')
    with open('groups.txt','w') as f:
        for item in tumor_df.columns:
            f.write('{}\t1\ttumor\n'.format(item))
        for item in normal_df.columns:
            f.write('{}\t2\tnormal\n'.format(item))
    with open('comps.txt','w') as f:
        f.write('1\t2\n')

    # run altanalyze, assuming singularity is enabled
    os.makedirs('./diff_dir/altanalyze_output/ExpressionInput')
    os.rename('exp.original-steady-state.txt','./diff_dir/altanalyze_output/ExpressionInput/exp.original-steady-state.txt')
    os.rename('groups.txt','./diff_dir/groups.txt')
    os.rename('comps.txt','./diff_dir/comps.txt')
    old_dir = os.getcwd()
    os.chdir('./diff_dir')
    writable_path = '/gpfs/data/yarmarkovichlab/Frank/SNAF_ecosystem/altanalyze'
    cmd = 'singularity run -B $PWD:/mnt --writable {} DE altanalyze_output groups.txt'.format(writable_path)
    subprocess.run(cmd,shell=True)
    os.chdir(old_dir)
    if is_te:
        os.rename('./diff_dir/altanalyze_output/ExpressionInput/DEGs-LogFold_0.0_adjp/GE.tumor_vs_normal.txt','DE_result_TE_chi_{}.txt'.format(cancer))
    else:
        os.rename('./diff_dir/altanalyze_output/ExpressionInput/DEGs-LogFold_0.0_adjp/GE.tumor_vs_normal.txt','DE_result_TE_{}.txt'.format(cancer))
    subprocess.run('rm -rf ./diff_dir',shell=True)

def run_splicing_de(events,cancer,te):
    adata_gtex = ad.read_h5ad('/gpfs/data/yarmarkovichlab/Frank/SNAF_ecosystem/data/controls/GTEx_junction_counts.h5ad')
    adata_tcga = ad.read_h5ad('/gpfs/data/yarmarkovichlab/Frank/SNAF_ecosystem/data/controls/tcga_matched_control_junction_count.h5ad')
    # build gtex_map
    gtex_map = pd.read_csv('/gpfs/data/yarmarkovichlab/chordoma/NeoVerse_analysis/splicing_annotations/gtex_hg38_t2t/tmp_prelift.bed',sep='\t',header=None)
    gtex_map.columns = ['chrom','start_1','end','uid']
    gtex_mapping = {}
    for row in gtex_map.itertuples():
        gtex_mapping['{}:{}-{}'.format(row.chrom,str(int(row.start_1+1)),row.end)] = row.uid
    # build tcga_map
    tcga_map = pd.read_csv('/gpfs/data/yarmarkovichlab/chordoma/NeoVerse_analysis/splicing_annotations/tcga_hg38_t2t/tmp_prelift.bed',sep='\t',header=None)
    tcga_map.columns = ['chrom','start_1','end','uid']
    tcga_mapping = {}
    for row in tcga_map.itertuples():
        tcga_mapping['{}:{}-{}'.format(row.chrom,str(int(row.start_1+1)),row.end)] = row.uid

    # extract normal
    final_data = np.empty((len(events),adata_gtex.shape[1]+adata_tcga.shape[1]),dtype=np.float32)
    for i,event in tqdm(enumerate(events)):
        store = []
        # gtex
        uid = gtex_mapping.get(event,None)
        if uid is not None:
            tmp = adata_gtex[uid,:].X.toarray().reshape(-1)
            store.append(tmp)
        else:
            tmp = np.zeros(adata_gtex.shape[1])
            store.append(tmp)
        # tcga
        uid = tcga_mapping.get(event,None)
        if uid is not None:
            tmp = adata_tcga[uid,:].X.toarray().reshape(-1)
            store.append(tmp)
        else:
            tmp = np.zeros(adata_tcga.shape[1])
            store.append(tmp)
        store_data = np.concatenate(store)
        final_data[i,:] = store_data
    normal_df = pd.DataFrame(data=final_data,index=events,columns=['normal{}'.format(i+1) for i in range(final_data.shape[1])])

    # extract tumor
    all_splicing = pd.read_csv(os.path.join(root_atlas_dir,cancer,'splicing_all.txt'),sep='\t')
    tmp = all_splicing.loc[all_splicing['Unnamed: 0'].isin(events),:]
    tmp.columns = ['event','count','sample']
    tumor_df = tmp.groupby(by=['event','sample'])['count'].sum().unstack(fill_value=0)

    # run mannwhitneyU
    rawps = []
    for event in tqdm(events):
        if event in tumor_df.index:
            t = tumor_df.loc[event,:].values
        else:
            t = np.zeros(tumor_df.shape[1])
        if event in normal_df.index:
            n = normal_df.loc[event,:].values
        else:
            n = np.zeros(normal_df.shape[1])
        u1,p = mannwhitneyu(t,n,alternative='greater')
        rawps.append(p)
    _,adjps,_,_ = multipletests(rawps,alpha=0.05,method='fdr_bh',is_sorted=False,returnsorted=False)
    de_df = pd.DataFrame(data={'rawp':rawps,'adjp':adjps},index=events)
    if te:
        de_df.to_csv('DE_splicing_te_{}.txt'.format(cancer),sep='\t')
    else:
        de_df.to_csv('DE_splicing_{}.txt'.format(cancer),sep='\t')

    
def run_nuorf_de(nuorfs,cancer,method):

    # extract normal df, percentile
    normal = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/code/hla_ligand_atlas_now_0.05_tesorai.txt',sep='\t')
    n_normal = len(normal['raw_file'].unique())
    normal = normal.loc[normal['peptide'].isin(nuorfs),:]

    # extract tumor df
    tumor_dfs = []
    tumor_raws = []
    giant_msms = pd.read_csv(os.path.join(root_atlas_dir,cancer,'antigen','fdr','msmsScans_all_add_tesorai.txt'),sep='\t')
    for raw,msms in giant_msms.groupby(by='Raw file'):
        msms = msms.loc[msms['Identified']=='+',:]
        each_raw_data = []
        msms['Precursor intensity'] = msms['Precursor intensity'].fillna(value=0)
        msms = msms.sort_values(by='Precursor intensity',ascending=True)
        msms['percentile'] = [(i+1)/msms.shape[0] for i in range(msms.shape[0])]
        for p,sub_df2 in msms.groupby(by='Sequence'):
            intensity = sub_df2['Precursor intensity'].values.max()
            percentile = sub_df2['percentile'].values.max()
            each_raw_data.append((p,intensity,percentile))
        each_raw_df = pd.DataFrame.from_records(data=each_raw_data,columns=['peptide','intensity','percentile'])
        each_raw_df = each_raw_df.loc[each_raw_df['intensity']>0,:]
        if each_raw_df.shape[0] > 0:
            upper = np.quantile(each_raw_df['intensity'].values,0.75)
            each_raw_df['norm'] = np.log2(each_raw_df['intensity'].values/upper)
        else:   # if no peptide identified in that raw file
            each_raw_df['norm'] = []
        tumor_dfs.append(each_raw_df)
        tumor_raws.append(raw)
    final = pd.concat(tumor_dfs,axis=0,keys=tumor_raws).reset_index(level=-2).rename(columns={'level_0':'raw_file'})
    final['log_intensity'] = np.log2(final['intensity'].values)
    tumor = final.loc[final['peptide'].isin(nuorfs),:]

    # GLM fit using beta family, w or w/o group
    # you need data.endog a series, data.exog a dataframe (have const column)
    rawps = []
    for nuorf in tqdm(nuorfs):
        t = tumor.loc[tumor['peptide']==nuorf,'percentile'].values
        if len(t) > 0:
            n = normal.loc[normal['peptide']==nuorf,'percentile'].values
            if len(n) == 0:
                n = np.full(n_normal,fill_value=1e-5)
            endog = pd.Series(data=np.concatenate([t,n]),name='y')
            total = len(t) + n_normal
            exog = pd.DataFrame(data={'const':[1]*total,'group':[0]*len(t)+[1]*n_normal})

            def lrt(endog,exog):
                # alt
                model = BetaModel(endog.values,exog.values)
                results = model.fit()
                alt_llf = results.llf

                # null
                model = BetaModel(endog.values,exog.values[:,0])
                results = model.fit()
                null_llf = results.llf

                # lrt
                lr = 2 * (alt_llf - null_llf)
                p = chi2.sf(lr,df=1)
                
                return p

            while True:
                try:
                    p = lrt(endog,exog)
                    break
                except:
                    continue

        else:
            p = 1
        
        rawps.append(p)
    
    # FDR
    _,adjps,_,_ = multipletests(rawps,alpha=0.05,method='fdr_bh',is_sorted=False,returnsorted=False)
    de_df = pd.DataFrame(data={'rawp':rawps,'adjp':adjps},index=nuorfs)
    de_df.to_csv('DE_nuorf_{}.txt'.format(cancer),sep='\t')
    sys.exit('stop')









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
database_dir = '/gpfs/data/yarmarkovichlab/public/ImmunoVerse/database'
bayests_xy_path = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/gene/full_results_XY_essential_tissues.txt'
HG38_NORMAL_DIR = '/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/NeoVerse/GTEx/selected/hg38_telocal_intron'

# self_gene 
final = pd.read_csv('../final_all_ts_antigens.txt',sep='\t')
final = final.loc[final['typ']=='self_gene',:]
total_ensg_intra = list(set(final['ensgs'].values))
mem = pd.read_csv('../gene_morpheus_mem.txt',sep='\t',index_col=[0,1])
total_ensg_mem = mem.index.to_frame(index=False)[0][1:].values.tolist()
total_ensg = list(set(total_ensg_intra + total_ensg_mem))

# for cancer in cancers:
#     run_self_gene_de(total_ensg,cancer)

# TE
final = pd.read_csv('../final_all_ts_antigens.txt',sep='\t')
final = final.loc[final['typ']=='self_translate_te',:]
col = []
for items in final['source']:
    valid_items = []
    for item in items.split(';'):
        if (not item.startswith('nc|')) and (not 'nuORF' in item) and (not 'TE_info' in item):
            valid_items.append(item)
    designated_te = None
    min_lfc = 1e5
    for item in valid_items:
        lfc = float(item.split('|')[4])
        if lfc < min_lfc:
            designated_te = item
            min_lfc = lfc
    col.append(designated_te)
final['designated_te'] = col
final = final.loc[final['designated_te'].notna(),:]
total_tes = final['designated_te'].values.tolist()
total_tes = list(set([item.split('|')[0] for item in total_tes]))

# for cancer in cancers:
#     run_te_de(total_tes,cancer,False)


# splicing
final = pd.read_csv('../final_all_ts_antigens.txt',sep='\t')
final = final.loc[final['typ']=='splicing',:]
col = []
for items in final['source']:
    valid_items = []
    for item in items.split(';'):
        if (not item.startswith('nc|')) and (not 'nuORF' in item) and (not 'TE_info' in item):
            valid_items.append(item)
    designated_event = None
    min_lfc = 1e5
    for item in valid_items:
        lfc = float(item.split('|')[4])
        if lfc < min_lfc:
            designated_event = item
            min_lfc = lfc
    col.append(designated_event)
final['designated_event'] = col
final = final.loc[final['designated_event'].notna(),:]
total_events = final['designated_event'].values.tolist()
total_events = list(set([item.split('|')[0] for item in total_events]))

# for cancer in cancers:
#     run_splicing_de(total_events,cancer,False)


# TE chimeric part I 
final = pd.read_csv('../final_all_ts_antigens.txt',sep='\t')
final = final.loc[final['typ']=='TE_chimeric_transcript',:]
col = []
for items in final['source']:
    valid_items = []
    for item in items.split(';'):
        if (not item.startswith('nc|')) and (not 'nuORF' in item) and ('TE_info' in item):
            valid_items.append(item)
    designated_event = None
    min_lfc = 1e5
    for item in valid_items:
        lfc = float(item.split('|')[4])
        if lfc < min_lfc:
            designated_event = item
            min_lfc = lfc
    col.append(designated_event)
final['designated_event'] = col
final = final.loc[final['designated_event'].notna(),:]
total_events = final['designated_event'].values.tolist()
total_events = list(set([item.split('|')[0] for item in total_events]))

# for cancer in cancers:
#     run_splicing_de(total_events,cancer,True)

# TE chimeric part II
final = pd.read_csv('../final_all_ts_antigens.txt',sep='\t')
final = final.loc[final['typ']=='TE_chimeric_transcript',:]
col = []
for items in final['source']:
    valid_items = []
    for item in items.split(';'):
        if (not item.startswith('nc|')) and (not 'nuORF' in item) and (not 'TE_info' in item):
            valid_items.append(item)
    designated_te = None
    min_lfc = 1e5
    for item in valid_items:
        lfc = float(item.split('|')[4])
        if lfc < min_lfc:
            designated_te = item
            min_lfc = lfc
    col.append(designated_te)
final['designated_te'] = col
final = final.loc[final['designated_te'].notna(),:]
total_tes = final['designated_te'].values.tolist()
total_tes = list(set([item.split('|')[0] for item in total_tes]))

# for cancer in cancers:
#     run_te_de(total_tes,cancer,True)

# nuorf
final = pd.read_csv('../final_all_ts_antigens.txt',sep='\t')
final = final.loc[final['typ']=='nuORF',:]
nuorfs = list(set(final['pep'].values))


for cancer in cancers:
    run_nuorf_de(nuorfs,cancer)

sys.exit('stop')

col1 = []
col2 = []
for row in final.itertuples():
    if row.typ == 'self_translate_te':
        c = row.cancer
        source = row.designated_te
        source = source.split('|')[0]
        df = pd.read_csv('DE_result_TE_{}.txt'.format(c),sep='\t',index_col=0)
        rawp = df.loc[source,:]['rawp']
        adjp = df.loc[source,:]['adjp']
        col1.append(rawp)
        col2.append(adjp)
    else:
        col1.append(None)
        col2.append(None)
final['ts_rawp'] = col1
final['ts_adjp'] = col2
final.to_csv('check.txt',sep='\t',index=None)