#!/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/bin/python3.7

import pandas as pd
import numpy as np
import sys,os
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import SeqIO
import subprocess
import matplotlib.pyplot as plt
import matplotlib as mpl
from tqdm import tqdm
import pycircos
from Bio.Seq import Seq
import pickle
import re

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
root_atlas_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas'


# get ncbi derived proteome
genome_dict = {}
with open('GCA_013267435.1/GCA_013267435.1_ASM1326743v1_genomic.fna','r') as in_handle:
    for title,seq in SimpleFastaParser(in_handle):
        genome_dict[title] = seq
ref = genome_dict['CP053989.1 Niallia circulans strain FDAARGOS_783 chromosome, complete genome']

# df = pd.read_csv('GCA_013267435.1/genomic.gtf',sep='\t',skiprows=5,header=None)
# df.columns = ['chrom','source','feature','start','end','score','strand','phase','attribute']
# df = df.loc[(df['source']=='Protein Homology') & (df['feature']=='CDS'),:]
# pat1 = re.compile(r'protein_id "(.+?)"')
# pat2 = re.compile(r'product "(.+?)"')
# ncbi_pro = {}
# for item1,item2,item3,item4 in zip(df['start'],df['end'],df['strand'],df['attribute']):
#     start = int(item1)
#     end = int(item2)
#     strand = item3
#     try:
#         pid = re.search(pat1,item4).group(1)
#     except:
#         pid = 'unknown'
#     try:
#         pname = re.search(pat2,item4).group(1)
#     except:
#         pname = 'unknown'

#     title = '{}|{}|{}|{}|{}'.format(pid,str(start),str(end),strand,pname.replace(' ','_'))
    
#     cds = ref[start-1:end]
#     if len(cds) % 3 == 0:
#         if strand == '+':
#             aa = str(Seq(cds).translate(table=11,to_stop=False))
#         else:
#             aa = str(Seq(cds).reverse_complement().translate(table=11,to_stop=False))
#     else:
#         continue
#     if '*' not in aa and 'unknown' not in title:
#         ncbi_pro[title] = aa

# uniprot_pro = {}
# with open('Niallia_circulans_UP000319837.fasta','r') as in_handle:
#     for title,seq in SimpleFastaParser(in_handle):
#         uniprot_pro[title] = seq

# pd.Series(ncbi_pro).to_csv('ncbi_pro.txt',sep='\t')
# pd.Series(uniprot_pro).to_csv('uniprot_pro.txt',sep='\t')

# with open('Niallia_circulans_ncbi.fasta','w') as f:
#     for k,v in ncbi_pro.items():
#         f.write('>{}\n{}\n'.format(k,v))

# draw circos plot for peptide
result = pd.read_csv('tesorai_peptide_fdr_OV_ncbi_pathogen.tsv',sep='\t')
cond = []
for item in result['possible_protein_ids']:
    
    if isinstance(item,str):
        sources = item.split(';;')
        identities = []
        for source in sources:
            if source.startswith('QKH'):
                identities.append('bingo')
            else:
                identities.append('nah')
        identities = list(set(identities))

        if len(identities) == 1 and identities[0] == 'bingo':
            cond.append(True)
        else:
            cond.append(False)
    else:
        cond.append(False)

result = result.loc[cond,:]
proteins = []
for item1,item2 in zip(result['clean_sequence'],result['possible_protein_ids']):
    for source in item2.split(';;'):
        proteins.append(source)
col1_p = []
col2_p = []
for pro in proteins:
    start = int(pro.split('|')[1]) - 1
    end = int(pro.split('|')[2])
    width = end - start + 1
    col1_p.append(start)
    col2_p.append(width)


# # get manifest files to programmatically download
# tcga_ids = [
#     'TCGA-13-0899-01A',
#     'TCGA-23-1021-01B',
#     'TCGA-10-0931-01A',
#     'TCGA-LN-A49R-01A',
#     'TCGA-Q9-A6FW-01A',
#     'TCGA-IG-A4QS-01A',
#     'TCGA-BR-A4QM-01A',
#     'TCGA-CG-5721-01A',
#     'TCGA-HU-A4GC-01A'
# ]

# df_list = []
# for c in cancers:
#     cmd = 'find {} -type f -name "manifest_*_*.tsv"'.format(os.path.join(root_atlas_dir,c))
#     possible_f = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1][0]
#     df = pd.read_csv(possible_f,sep='\t')
#     df_list.append(df)
# manifest = pd.concat(df_list,axis=0,keys=cancers).reset_index(level=-2).rename(columns={'level_0':'cancer'})
# manifest = manifest.loc[manifest['sample_id'].isin(tcga_ids),:]
# manifest.to_csv('manifest.txt',sep='\t',index=None)


'''
Go to NCBI assembly to download the species reference assembly, including
genome fasta
annotation gtf
cds fasta
protein fasta

But there are 35 assembly for different strains

color code the cds

only look for ones that mapped to reference assembly
'''

genome_dict = {}
with open('GCA_013267435.1/GCA_013267435.1_ASM1326743v1_genomic.fna','r') as in_handle:
    for title,seq in SimpleFastaParser(in_handle):
        genome_dict[title] = seq
ref = genome_dict['CP053989.1 Niallia circulans strain FDAARGOS_783 chromosome, complete genome']

sample_dict = {
    'OV':['TCGA-10-0931-01A','TCGA-23-1021-01B','TCGA-13-0899-01A'],
    'ESCA':['TCGA-Q9-A6FW-01A','TCGA-IG-A4QS-01A','TCGA-LN-A49R-01A'],
    'STAD':['TCGA-BR-A4QM-01A','TCGA-HU-A4GC-01A','TCGA-CG-5721-01A']
}


# sample_all_true_seqs = {}
# for k,v in sample_dict.items():
#     all_true_seqs = []
#     for s in v:
#         print(s)
#         tmp = subprocess.run('grep -A1 \'kraken:taxid|1397\' {}/test_cseqs_1.fq'.format(s),shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#         seqs1 = [item for item in tmp if item != '--' and 'kraken' not in item]
#         tmp = subprocess.run('grep -A1 \'kraken:taxid|1397\' {}/test_cseqs_2.fq'.format(s),shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#         seqs2 = [item for item in tmp if item != '--' and 'kraken' not in item]
#         seqs = seqs1 + seqs2
#         true_seqs_f = [item for item in tqdm(seqs) if item in ref]
#         true_seqs_r = [str(Seq(item).reverse_complement()) for item in tqdm(seqs) if str(Seq(item).reverse_complement()) in ref]
#         true_seqs = true_seqs_f + true_seqs_r
#         all_true_seqs.extend(true_seqs)
#     sample_all_true_seqs[k] = all_true_seqs

# with open('circulans_sample_all_true_seq.p','wb') as f:
#     pickle.dump(sample_all_true_seqs,f)

with open('circulans_sample_all_true_seq.p','rb') as f:
    sample_all_true_seqs = pickle.load(f)

# draw cds as cytoband
df = pd.read_csv('GCA_013267435.1/genomic.gtf',sep='\t',skiprows=5,header=None)
df.columns = ['chrom','source','feature','start','end','score','strand','phase','attribute']
df_cds = df.loc[(df['attribute'].str.contains('16S ribosomal RNA')) & (df['feature']=='transcript'),:]
col1 = []
col2 = []
for start,end in zip(df_cds['start'],df_cds['end']):
    col1.append(int(start)-1)
    col2.append(int(end)-int(start)+1)

sample_all_true_seqs.pop('ESCA')
sample_all_true_seqs.pop('STAD')

Garc = pycircos.Garc
Gcircle = pycircos.Gcircle

circle = Gcircle(figsize=(6,6))
arc = Garc(arc_id='n_circulans',size=len(ref),raxis_range=(935,985),labelposition=150,label_visible=True)
circle.add_garc(arc)
circle.set_garcs(-90,245)
circle.tickplot('n_circulans',raxis_range=(985,1000),tickinterval=100000,ticklabels=None)
circle.barplot('n_circulans',data=[1]*len(col1),positions=col1,width=col2,raxis_range=(935,985),facecolor='blue')
for i,(k,v) in enumerate(sample_all_true_seqs.items()):
    pos = {}
    for item in v:
        i_ = ref.find(item)
        pos.setdefault(i_,0)
        pos[i_] += 1
    vmin,vmax = min(list(pos.values())),max(list(pos.values()))
    circle.barplot('n_circulans',data=list(pos.values()),positions=list(pos.keys()),width=75,rlim=[vmin-0.05,vmax+0.05],raxis_range=[935-80*(i+1),935-80*i],facecolor='black',spine=True)

circle.barplot('n_circulans',data=[1]*len(col1_p),positions=col1_p,width=col2_p,rlim=[-0.05,1.05],raxis_range=[775,855],facecolor='red',spine=True)

circle.save('n_circulans_plot_ov')

        