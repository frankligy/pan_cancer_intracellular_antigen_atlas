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
from scipy.stats import ttest_ind

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

root_atlas_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas'
human_canonical_proteome = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/gene/ensembl_protein.fasta'


# compartment analysis
normal = pd.read_csv('/gpfs/data/yarmarkovichlab/Aman/hla_ligand_atlas.txt',sep='\t')
cancer = pd.read_csv('all_unique_self_gene.txt',sep='\t')

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


data = []
tissues = []

total_map = {}

for t,sub_df in cancer.groupby(by='cancer'):
    tissues.append(t)
    tmp = []
    for v in keyword_map.values():
        v_ensg = pd.read_csv('./compartment/human_{}_protein_postdoc_final.txt'.format(v),sep='\t')
        v_ensg = v_ensg.loc[v_ensg['ensg'].notna(),:]['ensg'].values.tolist()

        total_map[v] = len(v_ensg)

        total_pep = set()

        for item1, item2 in zip(sub_df['pep'],sub_df['ensgs']):
            if item2 in v_ensg and item1 not in total_pep:
                total_pep.add(item1)
        tmp.append(len(total_pep))
    data.append(tmp)

total = sum(list(total_map.values()))
sf = np.array([i/total for i in total_map.values()]).reshape(1,-1)

data = np.array(data)
norm_data = data / sf

data = data / data.sum(axis=1).reshape(-1,1)
norm_data = norm_data / norm_data.sum(axis=1).reshape(-1,1)

pd.DataFrame(data=data,index=tissues,columns=list(keyword_map.values())).to_csv('./compartment/cancer_data_prop.txt',sep='\t')
pd.DataFrame(data=norm_data,index=tissues,columns=list(keyword_map.values())).to_csv('./compartment/cancer_norm_data_prop.txt',sep='\t')


for n,d in enumerate([data,norm_data]):
    fig,ax = plt.subplots()
    ax.bar(tissues,d[:,0],label=list(keyword_map.values())[0])
    ax.bar(tissues,d[:,1],bottom=d[:,0],label=list(keyword_map.values())[1])
    for i in range(2,10):
        ax.bar(tissues,d[:,i],bottom=d[:,:i].sum(axis=1),label=list(keyword_map.values())[i])
    ax.set_ylabel('Proportion')
    ax.set_ylim([0,1.05])
    ax.set_xticklabels(tissues,rotation=60)
    ax.legend()
    plt.savefig('./compartment/cancer_{}_prop.pdf'.format(n),bbox_inches='tight')
    plt.close()

sys.exit('stop')

# derive a p-value for 85 with at least one signal derived
# metric_df = pd.read_csv('signal_potential_useful.txt',sep='\t',index_col=0)
# metric_df = metric_df.sort_values(by='potential_in_signal',ascending=False)
# metric_df = metric_df.iloc[:85]
# print(ttest_ind(metric_df['potential_in_signal'],metric_df['potential_not_in_signal']))

mem_path_dic = {
    'cell_membrane': 'human_membrane_protein_postdoc_final_no_edit.txt',
    'ER':'human_er_protein_postdoc_final.txt',
    'golgi':'human_golgi_protein_postdoc_final.txt',
    'lysosome':'human_lysosome_protein_postdoc_final.txt',
    'secreted':'human_secreted_protein_postdoc_final.txt',
    'endosome':'human_endosome_protein_postdoc_final.txt'
}

# for mem_id, mem_path in mem_path_dic.items():
#     mem = pd.read_csv(mem_path,sep='\t')
#     mem = mem.loc[mem['ensg'].notna(),:]
#     membrane_ensg = mem['ensg'].values.tolist()

#     with open(human_canonical_proteome,'r') as in_handle, open('./signal/{}_protein.fasta'.format(mem_id),'w') as f:
#         for title,seq in SimpleFastaParser(in_handle):
#             ensg,enst,gs = title.split('|')
#             if ensg in membrane_ensg:
#                 f.write('>{}\n{}\n'.format(title,seq))


normal = pd.read_csv('/gpfs/data/yarmarkovichlab/Aman/hla_ligand_atlas.txt',sep='\t')
cancer = pd.read_csv('all_unique_self_gene.txt',sep='\t')

for mem_id, mem_path in mem_path_dic.items():
    dic = {}
    with open('./signal/{}_protein.fasta'.format(mem_id),'r') as in_handle:
        for title,seq in SimpleFastaParser(in_handle):
            dic[title] = seq


    final = {}
    ensg2symbol = {}
    final2prop = {}
    p = './signal/{}.signalp5'.format(mem_id+'_protein')
    with open(p,'r') as f:
        for line in f:
            if line.startswith('ENSG') and 'SP(Sec/SPI)' in line:
                line = line.rstrip('\n')
                k = line.split('\t')[0].split('|')[0]
                gs = line.split('\t')[0].split('|')[2]
                uid = line.split('\t')[0]
                tmp = line.split('CS pos: ')[1].split('. ')[0]
                p1,p2 = tmp.split('-')
                p1,p2 = int(p1),int(p2)
                signal_seq = dic[uid][:p1]
                prop = (p1,len(dic[uid])-p1)
                ensg2symbol[k] = gs
                final[k] = signal_seq
                final2prop[k] = prop


    # in signal
    data = []
    for item1,item2 in zip(normal['peptide'],normal['ensg']):
        signal_seq = final.get(item2,None)
        if signal_seq is not None:
            if item1 in signal_seq:
                data.append([item1,item2,ensg2symbol[item2],'normal'])

    for item1,item2,item3 in zip(cancer['pep'],cancer['ensgs'],cancer['cancer']):
        signal_seq = final.get(item2,None)
        if signal_seq is not None:
            if item1 in signal_seq:
                data.append([item1,item2,ensg2symbol[item2],item3])

    df = pd.DataFrame(data=data,index=np.arange(len(data)),columns=['pep','ensg','gs','source'])
    df = df.drop_duplicates()
    df.to_csv('./signal/all_signal_peptides_{}.txt'.format(mem_id),sep='\t',index=None)
    df_in_signal = df

    # not in signal
    data = []
    for item1,item2 in zip(normal['peptide'],normal['ensg']):
        signal_seq = final.get(item2,None)
        if signal_seq is not None:
            if item1 not in signal_seq:
                data.append([item1,item2,ensg2symbol[item2],'normal'])

    for item1,item2,item3 in zip(cancer['pep'],cancer['ensgs'],cancer['cancer']):
        signal_seq = final.get(item2,None)
        if signal_seq is not None:
            if item1 not in signal_seq:
                data.append([item1,item2,ensg2symbol[item2],item3])

    df = pd.DataFrame(data=data,index=np.arange(len(data)),columns=['pep','ensg','gs','source'])
    df = df.drop_duplicates()
    df.to_csv('./signal/all_signal_peptides_not_{}.txt'.format(mem_id),sep='\t',index=None)
    df_not_in_signal = df


    # calculate potential
    metrics = []
    for k,v in final2prop.items():
        n_in_signal = len(df_in_signal.loc[df_in_signal['ensg']==k,:]['pep'].unique())
        n_not_in_signal = len(df_not_in_signal.loc[df_not_in_signal['ensg']==k,:]['pep'].unique())
        pot_in_singal = n_in_signal / v[0]
        pot_not_in_signal = n_not_in_signal / v[1]
        metrics.append([n_in_signal,n_not_in_signal,pot_in_singal,pot_not_in_signal])
    metrics_df = pd.DataFrame(data=metrics,index=[ensg2symbol[item] for item in final2prop.keys()],columns=['n_in_signal','n_not_in_signal','potential_in_signal','potential_not_in_signal'])
    metrics_df.to_csv('./signal/signal_potential_{}.txt'.format(mem_id),sep='\t')
    metrics_df = metrics_df.loc[(metrics_df['n_in_signal']>0) | (metrics_df['n_not_in_signal']>0),:]
    metrics_df.to_csv('./signal/signal_potential_useful_{}.txt'.format(mem_id),sep='\t')


    metrics_df.sort_values(by='potential_in_signal',ascending=False,inplace=True)
    # making the line plots
    y1 = []
    y2 = []
    for i1,i2 in zip(metrics_df['potential_in_signal'],metrics_df['potential_not_in_signal']):
        y1.append(i1)
        y2.append(i2)
    x = np.arange(len(y1))
    fig,ax = plt.subplots()
    ax.plot(x,y1,lw=1,marker='o',markersize=1,label='signal_region')
    ax.plot(x,y2,lw=1,marker='o',markersize=1,label='non_signal_region')
    ax.set_ylabel('presentation_potential')
    ax.set_xlabel('membrane_proteins')
    ax.set_xticks([])
    ax.legend()
    plt.savefig('./signal/signal_peptide_potential_{}.pdf'.format(mem_id),bbox_inches='tight')
    plt.close()
    # making the boxplot
    # metrics_df = metrics_df.loc[metrics_df['n_in_signal']>0,:]
    s,p = ttest_ind(metrics_df['potential_in_signal'],metrics_df['potential_not_in_signal'])
    fig,ax = plt.subplots()
    bp = ax.boxplot(x=[metrics_df['potential_in_signal'].values,metrics_df['potential_not_in_signal'].values],positions=[0,1],patch_artist=True)
    for flier in bp['fliers']:
        flier.set_markersize(1)
        flier.set_marker('o')
    for box in bp['boxes']:
        box.set_facecolor('green')
        box.set_edgecolor('black')
        box.set_linewidth(1)
    ax.set_title(mem_id)
    fig.suptitle('{}_{}'.format(p,metrics_df.shape[0]))
    plt.savefig('./signal/signal_peptide_potential_boxplot_{}.pdf'.format(mem_id),bbox_inches='tight')
    plt.close()

# assemble plots
data = []
for mem_id in mem_path_dic.keys():
    metrics_df = pd.read_csv('./signal/signal_potential_useful_{}.txt'.format(mem_id),sep='\t')
    # metrics_df = metrics_df.loc[metrics_df['n_in_signal']>0,:]
    for item in metrics_df['potential_in_signal']:
        data.append((mem_id,item,'signal'))
    for item in metrics_df['potential_not_in_signal']:
        data.append((mem_id,item,'not_signal'))  
final = pd.DataFrame(data=data,columns=['location','potential','identity'])
fig,ax = plt.subplots()
sns.boxplot(final,x='location',y='potential',hue='identity')
plt.savefig('./signal/assembled_plot.pdf',bbox_inches='tight')
plt.close()








    



