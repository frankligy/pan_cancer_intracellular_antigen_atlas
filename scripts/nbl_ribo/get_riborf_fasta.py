#!/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/bin/python3.7

import pandas as pd
import numpy as np
import sys,os
from tqdm import tqdm
import mygene
import subprocess
import re
import argparse
import pickle
import pysam
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import bisect
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from collections import Counter

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

VALID_CHROM = set(['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12',
                   'chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM'])

def process_gtf(gtf_path):
    gtf = pd.read_csv(gtf_path,sep='\t',header=None)
    gtf.columns = ['chrom', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attrs']
    dic = {}  # {tid:[+,(),(),()]}
    pat = re.compile(r'transcript_id "(ENST\d+)";')
    col = []
    for item in gtf['attrs']:
        tid = re.search(pat,item).group(1)
        col.append(tid)
    gtf['tid'] = col
    for tid,sub_df in tqdm(gtf.groupby(by='tid')):
        sub_df = sub_df.loc[sub_df['type']=='exon',:]
        strand = sub_df['strand'].iloc[0]
        chrom = sub_df['chrom'].iloc[0]
        tup = [strand,chrom]
        for row in sub_df.itertuples():   # no need to process the first line which is transcript
            tup.append((row.start,row.end))
        dic[tid] = tup
    return dic

def process_fasta(fa_path):
    dict_fa = {}
    with open(fa_path,'r') as in_handle:
        for title,seq in SimpleFastaParser(in_handle):
            if title in VALID_CHROM:
                dict_fa[title] = seq
    return dict_fa

def row2orf(row):
    uid = row.Index
    item_list = uid.split(':')
    enst,chrom,strand = item_list[0],item_list[1],item_list[2]
    if chrom not in VALID_CHROM:
        orf_seq = 'NOT VALID CHROM'
        return orf_seq
    strand = strand.split('|')[0]
    start,end = row.codon5,row.codon3
    exons = gtf_dic.get(enst,None)
    if exons is None:  # can not recover
        orf_seq = 'ENST NOT IN HG19 ANNOTATION'
        return orf_seq
    else:
        assert exons[0] == strand
        assert exons[1] == chrom
        exons = exons[2:]
        pos = []
        for exon in exons:
            pos.extend(list(exon))
        start = int(start)
        end = int(end)
        s_pos = bisect.bisect(pos,start)
        e_pos = bisect.bisect(pos,end)
        try:
            assert s_pos % 2 == 1
            assert e_pos % 2 == 1
        except AssertionError:
            if s_pos == e_pos:
                if strand == '-':
                    orf_seq = fa_dic[chrom][start-1:end]
                    orf_seq = str(Seq(orf_seq).reverse_complement())
                elif strand == '+':
                    orf_seq = fa_dic[chrom][start-1:end]
                return orf_seq
            else:    # no way to recover
                orf_seq = 'FALL INTO INTRON'
                return orf_seq
        else:
            s_n_exon = (s_pos + 1) // 2
            e_n_exon = (e_pos + 1) // 2
            effect_exons = exons[s_n_exon-1:e_n_exon]
            total = len(effect_exons)
            if strand == '-':
                if total == 1:  # reside on one exon
                    orf_seq = fa_dic[chrom][start-1:end]
                    orf_seq = str(Seq(orf_seq).reverse_complement())
                else:   # spread over more than one exon
                    orf_seq = ''
                    for i,exon in enumerate(effect_exons[::-1]):
                        if i == 0:
                            s = exon[0]
                            e = end
                        elif i == total - 1:
                            s = start
                            e = exon[1]
                        else:
                            s,e = exon
                        seq = fa_dic[chrom][s-1:e]
                        seq = str(Seq(seq).reverse_complement())
                        orf_seq += seq

            elif strand == '+':
                if total == 1:
                    orf_seq = fa_dic[chrom][start-1:end]
                else:
                    orf_seq = ''
                    for i,exon in enumerate(effect_exons):
                        if i == 0:
                            s = start
                            e = exon[1]
                        elif i == total - 1:
                            s = exon[0]
                            e = end
                        else:
                            s,e = exon
                        seq = fa_dic[chrom][s-1:e]
                        orf_seq += seq
            return orf_seq


        
        


HG19_SEQ = "/gpfs/data/yarmarkovichlab/neuroblastoma/riboseq/riborf_test/hg19.fa"
HG19_GTF = "/gpfs/data/yarmarkovichlab/neuroblastoma/riboseq/riborf_test/hg19.ensGene.gtf"
HG19_DIC = "/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/NeoVerse/nuORF/gtf_dic.p"
RIBORF_RESULT = "/gpfs/data/yarmarkovichlab/neuroblastoma/riboseq/ribo_real/OHMX20210030_001/outputDir/pred.pvalue.parameters.txt"
FASTA_OUTDIR = "/gpfs/data/yarmarkovichlab/neuroblastoma/riboseq/ribo_real/OHMX20210030_001/outputDir/"


# gtf_dic = process_gtf('hg19.ensGene.gtf')  # download from ucsc
# with open(HG19_DIC,'wb') as f:
#     pickle.dump(gtf_dic,f)

fa_dic = process_fasta(HG19_SEQ)   # download from ucsc
with open(HG19_DIC,'rb') as f:
    gtf_dic = pickle.load(f)

df = pd.read_csv(RIBORF_RESULT,sep='\t',index_col=0)   # from nuORF paper
df.rename(columns={'pred.pvalue':'pvalue'},inplace=True)
df['plotType'] = [item.split('|')[-2] for item in df.index]
df = df.loc[df['plotType'].isin(['dORF','odORF','ouORF','uORF','iORF','noncoding']),:]
df = df.loc[[item.startswith('ENST') for item in df.index],:]

data = []
for row in tqdm(df.itertuples(),total=df.shape[0]):
    orf_seq = row2orf(row)
    ori_orf_length = row.length
    data.append((row.Index,orf_seq,ori_orf_length,row.pvalue))
result = pd.DataFrame.from_records(data=data,columns=['uid','orf_seq','ori_orf_length','translatability'])
result = result.loc[(result['orf_seq']!='ENST NOT IN HG19 ANNOTATION') & (result['orf_seq']!='FALL INTO INTRON') & (result['orf_seq']!='NOT VALID CHROM'),:]
result['now_orf_length'] = [len(item)-1 for item in result['orf_seq']]  # we know this is +1 issue
result['cond'] = result['ori_orf_length'] == result['now_orf_length']
result = result.loc[result['cond'],:]

actual_orf_seq = []
for item in result['orf_seq']:
    if item[:3] in set(['ATG','TTG','CTG','GTG']):
        actual_orf_seq.append(item[:-1])
    elif item[1:4] in set(['ATG','TTG','CTG','GTG']):
        actual_orf_seq.append(item[1:])
    else:
        actual_orf_seq.append('NO CLEAR NTG')   # PRICE pipeline allows hamming distance = 1
result['actual_orf_seq'] = actual_orf_seq
result = result.loc[result['actual_orf_seq'] != 'NO CLEAR NTG',:]

col = []
cond = []
for item in result['actual_orf_seq']:
    pep = str(Seq(item).translate(to_stop=False))
    col.append(pep)
    if ('*' in pep) and (not pep.endswith('*')):  # it might be ok that * is the last one
        cond.append(False)
    else:
        cond.append(True)
result['pep'] = col
result['no_stop'] = cond

result = result.loc[result['no_stop'],:]
result['actual_pep'] = [item.rstrip('*') if item.endswith('*') else item for item in result['pep']]
result.to_csv(os.path.join(FASTA_OUTDIR,'riborf_table.txt'),sep='\t',index=None)

with open(os.path.join(FASTA_OUTDIR,'riborf.fasta'),'w') as f:
    for uid,tran,pep in zip(result['uid'],result['translatability'],result['actual_pep']):
        if len(pep) >= 8:
            f.write('>{}|{}|{}\n{}\n'.format(uid,str(tran),'nuORF',pep))











