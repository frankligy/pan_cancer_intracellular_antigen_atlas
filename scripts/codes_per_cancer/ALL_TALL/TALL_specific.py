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
import anndata as ad
from scipy.sparse import csr_matrix
from collections import Counter
from ast import literal_eval
import argparse
import json

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

def get_mut_seq(effect,et,aa_seq):
    if et == 'missense_variant':
        pat = r'([ARNDCQEGHILKMFPSTWYV])(\d+)([ARNDCQEGHILKMFPSTWYV])'
        match = re.search(pat,effect)
        ref_aa = match.group(1)
        pos_aa = int(match.group(2))
        alt_aa = match.group(3)
        documented_ref_aa = aa_seq[pos_aa-1]
        try:
            assert documented_ref_aa == ref_aa
        except:
            variant_seq = 'unknown'
        else:
            if MAX_PEP_LEN > pos_aa:
                actual_needed = pos_aa
            else:
                actual_needed = MAX_PEP_LEN
            variant_seq = aa_seq[pos_aa-1-(actual_needed-1):pos_aa-1] + alt_aa + aa_seq[pos_aa:pos_aa+(MAX_PEP_LEN-1)]

    elif et == 'inframe_insertion':
        if effect.endswith('dup'):
            effect = effect.split('dup')[0]
            if '_' in effect:
                first,second = effect.split('_')
                first = int(first[1:])
                second = int(second[1:])
            else:
                first = int(effect[1:])
                second = first
            preceding = aa_seq[:first-1]
            following = aa_seq[second:]
            impacted = aa_seq[first-1:second]
            updated = impacted + impacted
            n_aa_needed = MAX_PEP_LEN - 1
            variant_seq = preceding[-n_aa_needed:] + updated + following[:n_aa_needed]
        elif 'ins' in effect:
            effect,insert = effect.split('ins')
            if '_' in effect:  # found one in KIRC UHMK1@p.*420delinsYHL*YI, not considered this edge case for now
                first,second = effect.split('_')
                first = int(first[1:])
                second = int(second[1:])
                preceding = aa_seq[:first-1]
                following = aa_seq[second:]
                if len(aa_seq) >= second: # ov PCDH15@p.T1794_F1795insDS, I think due to reporting to another isoform
                    updated = aa_seq[first-1] + insert + aa_seq[second-1]
                    n_aa_needed = MAX_PEP_LEN - 1
                    variant_seq = preceding[-n_aa_needed:] + updated + following[:n_aa_needed]

    elif et == 'inframe_deletion':
        if 'delins' in effect:
            effect,delins = effect.split('delins')
            if '_' in effect:
                first,second = effect.split('_')
                first = int(first[1:])
                second = int(second[1:])
            else:
                first = int(effect[1:])
                second = int(effect[1:])
            preceding = aa_seq[:first-1]
            following = aa_seq[second:]
            updated = delins
            n_aa_needed = MAX_PEP_LEN - 1
            variant_seq = preceding[-n_aa_needed:] + updated + following[:n_aa_needed]
        elif effect.endswith('del'):
            effect = effect.split('del')[0]
            if '_' in effect:
                first,second = effect.split('_')
                first = int(first[1:])
                second = int(second[1:])
            else:      
                first = int(effect[1:])
                second = int(effect[1:])
            preceding = aa_seq[:first-1]
            following = aa_seq[second:]
            updated = ''
            n_aa_needed = MAX_PEP_LEN - 1
            variant_seq = preceding[-n_aa_needed:] + updated + following[:n_aa_needed]

    elif et == 'frameshift_variant':

        dict_fa = seq_dict['hg38']
        coord,replace,ensg = effect.split(';') 
        full_coord = coord
        before, after = replace.split('/')
        chrom,coord = coord.split(':')
        start,end = coord.split('-')
        start,end = int(start),int(end)

        # get mode from before after
        previous_length = 0 if before == '-' else len(before)
        after_length =  0 if after == '-' else len(after)
        if previous_length < after_length:  # insertion
            mode = 'mode1'
        elif previous_length > after_length:  # deletion
            mode = 'mode2'

        # get from cds file
        cds = cds_all.loc[(cds_all['ensg']==ensg) & (cds_all['is_canonical']),:]
        s = cds.iloc[0]
        strand = s['strand']
        cds_seq = s['cds_seq']
        coord = literal_eval(s['coord'])
        dic_info = literal_eval(s['dic'])

        if mode == 'mode2':
            start = start - 1
            end = end + 1
            after = ''

        # now let's go
        exons = []
        total_exon = len(coord)
        for exon in coord:
            exons.extend(list(exon))

        if strand == '+':
            pos = bisect.bisect_left(exons,start)
            if pos % 2 == 1:
                n_exon = pos // 2 + 1
            else:
                pos = bisect.bisect_right(exons,start)
                n_exon = pos // 2 + 1

            i1,i2,i3,i4 = dic_info[n_exon]
            i1,i2,i3,i4 = int(i1),int(i2),int(i3),int(i4)
            start_cds_pos = i2 + (start - i1)   # the start is this nt of the total cds
            end_cds_pos = i2 + (end - i1)     # the end is this nt of the total cds
            span = start - i1 + 1 + i4
            added_residues = span // 3
            now_phase = span % 3
            n_codon = i3 + added_residues 
            if now_phase != 0:
                n_codon = n_codon + 1
                p_codon = now_phase
            else:
                n_codon = n_codon
                p_codon = 3
            


        if strand == '-':
            exons = sorted(exons)
            pos = bisect.bisect_left(exons,end)
            if pos % 2 == 1:
                n_exon = pos // 2 + 1
            else:
                pos = bisect.bisect_right(exons,start)
                n_exon = pos // 2 + 1
            n_exon = total_exon - n_exon + 1
            i1,i2,i3,i4 = dic_info[n_exon]
            i1,i2,i3,i4 = int(i1),int(i2),int(i3),int(i4)
            start_cds_pos = i2 + (i1 - end)   # reverse
            end_cds_pos = i2 + (i1 - start)     # reverse
            span = i1 - end + 1 + i4
            added_residues = span // 3
            now_phase = span % 3
            n_codon = i3 + added_residues 
            if now_phase != 0:
                n_codon = n_codon + 1
                p_codon = now_phase
            else:
                n_codon = n_codon
                p_codon = 3        

            

        '''

        using ACVR2A@p.K437Rfs*5, chr2:147926117-147926117;A/-, ENSG00000121989 as example
        put the coordinate in ucsc, realize it is mode2 (16 and 18 as start and end, A become ''), start 16 occur in codon 434 as the 3rd one
        all the things are in the sense of translation

        so 'ACVR2A':{'chr2:147926117-147926117':['+','mode2',434,3]},

        we will take the 433 aa, and then start-offset2 to point to the actual pos that needs to retrive dna seq

        also Rfs*5 means counting from the K437R---*

        '''

        

        # now you either use new method of cds sequence or bounce back to the old method, we need n_exon either way
        if strand == '+':
            new_cds = cds_seq[:start_cds_pos] + after + cds_seq[end_cds_pos-1:]
        else:
            new_cds = cds_seq[:start_cds_pos] + str(Seq(after).reverse_complement()) + cds_seq[end_cds_pos-1:]
        new_cds_pep = str(Seq(new_cds).translate(to_stop=False))

        if '*' in new_cds_pep:
            variant_seq = new_cds_pep.split('*')[0]
        else:
            preceding = aa_seq[:n_codon-1]
            n_aa_needed = MAX_PEP_LEN - 1
            first_part = preceding[-n_aa_needed:]
            if strand == '+':
                stretch = dict_fa[chrom][start-1-(p_codon-1):start] + after + dict_fa[chrom][end-1:end-1+1000]
            elif strand == '-':
                stretch = dict_fa[chrom][start-1-1000:start] + after + dict_fa[chrom][end-1:end+(p_codon-1)]
                stretch = str(Seq(stretch).reverse_complement())
            second_part = str(Seq(stretch).translate(to_stop=False)).split('*')[0]
            variant_seq = first_part + second_part


    return variant_seq


## main program
MAX_PEP_LEN = 15
MIN_PEP_LEN = 8
cds_all = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/codes/summary/hg38_cds_final.txt',sep='\t',index_col=0)
protein_fasta = '/gpfs/data/yarmarkovichlab/public/reference/ImmunoVerse_data/ensembl_protein.fasta'
dic = {}
gs2ensg = {}
with open(protein_fasta,'r') as in_handle:
    for title,seq in SimpleFastaParser(in_handle):
        ensg,enst,gs = title.split('|')
        gs2ensg[gs] = ensg
        dic[gs] = seq 


hg38_fasta = '/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/NeoVerse/GTEx/circRNA/hg38.fa'
hg38_dict = {}
VALID_CHROM = set(['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12',
                   'chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM'])
with open(hg38_fasta,'r') as in_handle:
    for title,seq in SimpleFastaParser(in_handle):
        if title in VALID_CHROM:
            hg38_dict[title] = seq
seq_dict = {}
seq_dict['hg38'] = hg38_dict

mutation = pd.read_csv('../../variants/ALL_TALL/ALL_TALL_mutations.txt',sep='\t')
mutation = mutation.loc[mutation['anno1']!='stop',:]
mutation1 = mutation.loc[mutation['anno1']!='frameshift',:]
mutation2 = mutation.loc[mutation['anno1']=='frameshift',:]

# get gs,effect,type tuple
shelf = []
for row in mutation1.itertuples():
    gs = row.gs
    effect = row.effect_new
    typ = row.anno2
    if typ == 'SNV':
        typ = 'missense_variant'
    elif typ == 'deletion':
        typ = 'inframe_deletion'
    elif typ == 'insertion':
        if 'delins' in effect:
            typ = 'inframe_deletion'
        else:
            typ = 'inframe_insertion'
    shelf.append((gs,effect.split('.')[1],typ))

# # use liftover ucsc web to map hg19 to hg38
# for row in mutation2.itertuples():
#     effect = row.effect_new
#     print(effect)


mapping = pd.read_csv('../../variants/ALL_TALL/mapping.txt',sep='\t',index_col=0)['hg38'].to_dict()
for row in mutation2.itertuples():
    effect = row.effect_new
    effect = mapping[effect]
    if effect != 'unknown':
        effect = '{};{}/{};{}'.format(effect,row.ref,row.alt,gs2ensg[row.gs])
        shelf.append((row.gs,effect,'frameshift_variant'))

'''
here, chr1:100 -/AC, means AC will be inserted after 100, instead of at 100
some missense are not canonical
some gs are not canonical alias
delis have single first
'''

with open('../../atlas/ALL_TALL/db_fasta/mutation.fasta','w') as f:
    for tup in shelf:
        k,v,typ = tup
        aa = dic.get(k,None)
        if aa is not None:
            variant_seq = get_mut_seq(v,typ,aa)
        else:
            variant_seq = 'unknown'
        if variant_seq == 'unknown':
            print(tup)
        else:
            f.write('>{}|{}|na|0.5|{}|chr0:1-999|A/G|{}\n{}\n'.format(k,v,gs2ensg[k],typ,variant_seq))