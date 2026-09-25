#!/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/bin/python3.7

import pandas as pd
import numpy as np
import sys,os
import subprocess
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import anndata as ad

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

# compare the expression amongst HNC, OV, NBL, CESC, BRCA

ensg = 'ENSG00000185686'
t_hnc = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas/HNSC/gene_tpm.txt',sep='\t',index_col=0).loc[ensg,:].values
t_ov = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas/OV/gene_tpm.txt',sep='\t',index_col=0).loc[ensg,:].values
t_nbl = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas/NBL/gene_tpm.txt',sep='\t',index_col=0).loc[ensg,:].values
t_cesc = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas/CESC/gene_tpm.txt',sep='\t',index_col=0).loc[ensg,:].values
t_brca = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas/BRCA/gene_tpm.txt',sep='\t',index_col=0).loc[ensg,:].values

# for t in [t_ov,t_brca,t_cesc]:
#     prop = np.count_nonzero(t > 20) / len(t)
#     print(prop)


# fig,ax = plt.subplots()
# sns.swarmplot(data={'Ovarian':t_ov,'Breast':t_brca,'Cervical':t_cesc},size=1,ax=ax)
# ax.set_ylim([-5,500])
# ax.set_ylabel('PRAME TPM')
# ax.set_xlabel('cancer')
# plt.savefig('PRAME_stripplot.pdf',bbox_inches='tight')
# plt.close()


GTEX_GENE = '/gpfs/data/yarmarkovichlab/chordoma/NeoVerse_analysis/bulk-gex_v8_rna-seq_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct'
GTEX_GENE_ALL_H5AD = '/gpfs/data/yarmarkovichlab/chordoma/NeoVerse_analysis/gtex_gene_all.h5ad'
plot_type = 'boxplot+boxplot'
image_format = 'pdf'
ax1_label = 'OV'


def run(ensg,symbol):
    gtex = pd.read_csv(GTEX_GENE,sep='\t',skiprows=2,index_col=0)
    cond = ~gtex.columns.isin(['Cells - EBV-transformed lymphocytes','Cells - Cultured fibroblasts','Testis'])
    gtex = gtex.loc[:,cond]
    gtex.index = [item.split('.')[0] for item in gtex.index]
    ensg2symbol = pd.Series(index=gtex.index.tolist(),data=gtex['Description'].tolist()).to_dict()
    series = gtex.loc[ensg,:].iloc[1:]


    # plot
    t_pt, n_pt = plot_type.split('+')
    fig = plt.figure(figsize=(15,6))
    gs = mpl.gridspec.GridSpec(nrows=1,ncols=2,width_ratios=(0.1,0.9),wspace=0.2)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1],sharey=ax1)
    if t_pt == 'swarmplot':
        sns.swarmplot(data=tumor_expr,color='red',ax=ax1)
    elif t_pt == 'boxplot':
        sns.boxplot(data=t_ov,color='red',ax=ax1)
    elif t_pt == 'swarmplot_cat':
        sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
        from colors import pick_n_colors
        select_colors = pick_n_colors(len(cat_dic))
        for i,(cat,samples) in enumerate(cat_dic.items()):
            samples = [item + ',TPM' for item in samples]
            tumor_expr = final.loc[ensg,samples].values
            sns.swarmplot(data=tumor_expr,color=select_colors[i],ax=ax1)
        import matplotlib.lines as mlines
        ax1.legend(handles=[mlines.Line2D([],[],marker='o',linestyle='',color=i) for i in select_colors],labels=list(cat_dic.keys()),bbox_to_anchor=(0,1),loc='upper right') 
    elif t_pt == 'swarmplot_cat_box':
        sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
        from colors import pick_n_colors
        select_colors = pick_n_colors(len(cat_dic))
        for i,(cat,samples) in enumerate(cat_dic.items()):
            samples = [item + ',TPM' for item in samples]
            tumor_expr = final.loc[ensg,samples].values
            sns.swarmplot(data=tumor_expr,color=select_colors[i],ax=ax1)
        import matplotlib.lines as mlines
        ax1.legend(handles=[mlines.Line2D([],[],marker='o',linestyle='',color=i) for i in select_colors],labels=list(cat_dic.keys()),bbox_to_anchor=(0,1),loc='upper right')
        tumor_expr = final.loc[ensg,:].values
        bp = ax1.boxplot(x=[tumor_expr],positions=[0],patch_artist=False)
        

    ax1.set_ylabel('TPM')
    ax1.set_xlabel(ax1_label)

    if n_pt == 'barplot':
        ax2.bar(x=np.arange(len(series)),height=series.values,width=0.8,color='green')
    elif n_pt == 'boxplot':
        # # build a h5ad
        # tmp_df = pd.read_csv(GTEX_GENE_ALL,sep='\t',index_col=0,skiprows=2).iloc[:,1:]
        # adata = ad.AnnData(X=csr_matrix(tmp_df.values),obs=pd.DataFrame(index=tmp_df.index),var=pd.DataFrame(index=tmp_df.columns))
        # meta = pd.read_csv(GTEX_META,sep='\t',index_col=0)
        # common_samples = list(set(adata.var_names).intersection(set(meta.index)))
        # meta = meta.loc[common_samples,:]
        # mapping = meta['SMTSD'].to_dict()
        # adata.var['tissue'] = [mapping.get(item,'unknown') for item in adata.var_names]
        # adata.write('gtex_gene_all.h5ad')

        # use the h5ad
        adata = ad.read_h5ad(GTEX_GENE_ALL_H5AD)  # 56200 × 17382
        adata.obs_names = [item.split('.')[0] for item in adata.obs_names]
        adata.obs_names_make_unique()
        adata_gene = adata[[ensg],:]
        normal_expr_list = []
        for t in series.index.tolist():
            values = adata_gene[:,adata_gene.var['tissue']==t].X.toarray().reshape(-1)
            normal_expr_list.append(values)
        bp = ax2.boxplot(x=normal_expr_list,positions=np.arange(len(normal_expr_list)),patch_artist=True)
        for flier in bp['fliers']:
            flier.set_markersize(1)
            flier.set_marker('o')
        for box in bp['boxes']:
            box.set_facecolor('green')
            box.set_edgecolor('black')
            box.set_linewidth(1)

    ax2.set_xticks(np.arange(len(series)))
    ax2.set_xticklabels(series.index.tolist(),fontsize=10,rotation=90)
    ax2.set_xlabel('GTEx Normal')


    fig.suptitle('{},{}'.format(ensg,symbol))
    plt.savefig(os.path.join('.','{}_{}_ov_expr_{}.{}'.format(ensg,symbol,plot_type,image_format)),bbox_inches='tight')
    plt.close()


markers = {
    'ENSG00000185686':'PRAME'
}

for ensg,symbol in markers.items():
    run(ensg,symbol)



