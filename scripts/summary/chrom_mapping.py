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
import re
import bisect
import karyopype.karyopype as kp

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


# create cnv (amplification and deletion) region
pat = re.compile('(.+)([pq])(.+)')
ideogram = pd.read_csv('ideogram_ncbi_ucsc.txt',sep='\t')
regions = pd.read_csv('cnv_regions.txt',sep='\t',encoding_errors='replace')
cna = regions.loc[regions['type']=='amplification',:]
with open('./chrom_mapping/cna_region.bed','w') as f:
    for region in cna['notation']:
        match = re.search(pat,region)
        chrom,arm,band = match.groups()
        s = ideogram.loc[(ideogram['#chromosome']==chrom) & (ideogram['arm']==arm) & (ideogram['band']==float(band)),:].iloc[0]
        f.write('chr{} {} {}\n'.format(chrom,s['bp_start'],s['bp_stop']))

cnd = regions.loc[regions['type']=='deletion',:]
with open('./chrom_mapping/cnd_region.bed','w') as f:
    for region in cnd['notation']:
        match = re.search(pat,region)
        if not match:
            print(region);sys.exit('stop')
        chrom,arm,band = match.groups()
        s = ideogram.loc[(ideogram['#chromosome']==chrom) & (ideogram['arm']==arm) & (ideogram['band']==float(band)),:].iloc[0]
        f.write('chr{} {} {}\n'.format(chrom,s['bp_start'],s['bp_stop']))



root_atlas_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas'

prot_fa = {}
with open('/gpfs/data/yarmarkovichlab/public/ImmunoVerse/database/ensembl_protein.fasta','r') as in_handle:
    for title,seq in SimpleFastaParser(in_handle):
        ensg,enst,gs = title.split('|')
        prot_fa[ensg] = seq

cds_df = pd.read_csv('hg38_cds_final.txt',sep='\t',index_col=0)
cds_df = cds_df.loc[cds_df['is_canonical'],:]
cds_df = cds_df.loc[cds_df['dic'].notna(),:]
cds_df['dic'] = [literal_eval(item) for item in cds_df['dic']]
col = []
for dic in cds_df['dic']:
    col.append([exon[2] for e,exon in dic.items()])
cds_df['residue'] = col
ensg2strand = pd.Series(index=cds_df['ensg'].values,data=cds_df['strand'].values).to_dict()
ensg2dic = pd.Series(index=cds_df['ensg'].values,data=cds_df['dic'].values).to_dict()
ensg2residue = pd.Series(index=cds_df['ensg'].values,data=cds_df['residue'].values).to_dict()
ensg2chrom = pd.Series(index=cds_df['ensg'].values,data=cds_df['chrom'].values).to_dict()
valid_ensg = set(ensg2chrom.keys())
# cds_df.to_csv('check.txt',sep='\t')


for c in cancers:
    print(c)
    final_path = os.path.join(root_atlas_dir,c,'antigen','fdr','final_enhanced_all.txt')
    final = pd.read_csv(final_path,sep='\t')
    cond = [False if ('[]' in item) and ('(\'HLA-' not in item) else True for item in final['presented_by_each_sample_hla']]
    final = final.loc[cond,:]
    final = final.loc[final['typ']=='self_gene',:]
    final = final.loc[final['unique'],:]
    
    # create region bed
    with open('./chrom_mapping/{}_region.bed'.format(c),'w') as f:
        for row in final.itertuples():
            ensg = row.ensgs
            prot_seq = prot_fa[ensg]
            pos = prot_seq.index(row.pep) + 1
            if ensg in valid_ensg:
                dic = ensg2dic[ensg]
                strand = ensg2strand[ensg]
                chrom = ensg2chrom[ensg]
                residue = ensg2residue[ensg] 
                n_exon = bisect.bisect(residue,pos)   # residing on nth exon
                interval = pos - residue[n_exon-1] - 1   # from last completed codon to this codon, interval codons in between
                stop_phase = dic[n_exon][3]


                if strand == '+':
                    start_pos = dic[n_exon][0] - stop_phase + interval * 3
                    end_pos = start_pos + len(row.pep) * 3
                else:
                    start_pos = dic[n_exon][0] + stop_phase - interval * 3
                    end_pos = start_pos - len(row.pep) * 3
                    start_pos, end_pos = end_pos, start_pos

                '''
                141 stop_phase=2 start=47377014, interval=13 (155-141-1) retract 2, and then count 13 intervals  
                ***  **        * *** *** *** *** *** *** *** *** *** *** *** *** *** 
                '''

                f.write('{} {} {}\n'.format(chrom,start_pos-1,end_pos))

    # draw karyopype
    kp.plot_karyopype('hg38',['./chrom_mapping/cnd_region.bed','./chrom_mapping/{}_region.bed'.format(c)]) 
    os.rename('hg38_karyopype.pdf','./chrom_mapping/{}_karyopype_cnd.pdf'.format(c))

# all karyopype or comparision
# lis = ['./chrom_mapping/{}_region.bed'.format(c) for c in cancers]
# kp.plot_karyopype('hg38',lis)
# lis = ['./chrom_mapping/{}_region.bed'.format(c) for c in ['NBL','SKCM','AML']]
# kp.plot_karyopype('hg38',lis)

    

        



        
        

            
        