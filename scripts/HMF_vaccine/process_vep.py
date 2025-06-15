#!/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/bin/python3.7

import pandas as pd
import numpy as np
import sys,os
import subprocess
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import anndata as ad
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
import re
from ast import literal_eval
import bisect

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


# compare
# data = []
# with open('wgs_variant.fasta','r') as in_handle:
#     for title,seq in SimpleFastaParser(in_handle):
#         gs = title.split('|')[2]
#         effect = title.split('|')[6]
#         et = title.split('|')[7]
#         if 'missense_variant' in et:
#             et = 'missense_variant'
#         elif 'inframe_insertion' in et:
#             et = 'inframe_insertion'
#         elif 'inframe_deletion' in et:
#             et = 'inframe_deletion'
#         elif 'frameshift_variant' in et:
#             et = 'frameshift_variant'
#         data.append((et,gs,effect,gs+'_'+effect))
# df = pd.DataFrame.from_records(data,columns=['et','gs','effect','uid'])
# df.to_csv('dna.txt',sep='\t',index=None)

# data = []
# with open('/gpfs/data/yarmarkovichlab/HMF_vaccine/immunoverse_result/db_fasta/CPCT02030461T_HC3NNBGX9_S7/tmp_variants.fasta','r') as in_handle:
#     for title,seq in SimpleFastaParser(in_handle):
#         gs = title.split('|')[2]
#         effect = title.split('|')[4]
#         et = title.split('|')[5]
#         if 'missense_variant' in et:
#             et = 'missense_variant'
#         elif 'inframe_insertion' in et:
#             et = 'inframe_insertion'
#         elif 'inframe_deletion' in et:
#             et = 'inframe_deletion'
#         elif 'frameshift_variant' in et:
#             et = 'frameshift_variant'
#         data.append((et,gs,effect,gs+'_'+effect))
# df = pd.DataFrame.from_records(data,columns=['et','gs','effect','uid'])
# df.to_csv('rna.txt',sep='\t',index=None)

dna = pd.read_csv('dna.txt',sep='\t')['uid'].unique()
rna = pd.read_csv('rna.txt',sep='\t')['uid'].unique()
common = set(dna).intersection(set(rna))
somatic = pd.read_csv('/gpfs/data/yarmarkovichlab/HMF_vaccine/immunoverse_result/variants_possible_somatic.txt',sep='\t')
col = []
for info,gs in zip(somatic['INFO'],somatic['symbol']):
    dic = dict(pair.split('=') for pair in info.split(';'))
    effect = dic['CE'].split('_')[2]
    col.append(gs+'_'+effect)
somatic['uid'] = col
print(somatic)
somatic_sub = somatic.loc[somatic['uid'].isin(common),:]
print(somatic_sub)
sys.exit('sotp')







def variant_effect_analyzer(pos,change,enst,cat,impact,genetic,chrom,strand,coord,enst2protein,can_enst2protein,dict_fa,cds_all):

    protein = can_enst2protein.get(enst,None) if cat == 'canonical' else enst2protein.get(enst,None)

    if protein is None or '*' in change:
        return None
    else:
        if 'missense_variant' in impact:  # 45 E/G 


            pos = int(pos)
            ref,alt = change.split('/')

            previous = protein[:pos-1]
            after = protein[pos:]
            seq = previous[-(MAX_PEP_LEN-1):] + alt + after[:(MAX_PEP_LEN-1)]


        elif 'inframe_deletion' in impact:

            if '-' in pos: # 67-68|-/S
                first, second = pos.split('-')
                first, second = int(first), int(second)
                ref,alt = change.split('/')

            else:  # 446 G/EG
                first = int(pos)
                second = int(pos)
                ref,alt = change.split('/')
            
            previous = protein[:first-1]
            after = protein[second:]
            impacted = protein[first-1:second]

            if ref == '-':
                updated = impacted[0] + alt + impacted[1]
            else:
                updated = alt

            l = len(updated)
            seq = previous[-(MAX_PEP_LEN-l):] + updated + after[:(MAX_PEP_LEN-l)]


        elif 'inframe_insertion' in impact:

            if '-' in pos: # 722-723 EE/E
                first, second = pos.split('-')
                first, second = int(first), int(second)
                ref,alt = change.split('/')

            else:  # 229|E/-
                first = int(pos)
                second = int(pos)
                ref,alt = change.split('/')
            
            previous = protein[:first-1]
            after = protein[second:]
            impacted = protein[first-1:second]

            if alt == '-':
                updated = ''
            else:
                updated = alt

            l = len(updated)
            seq = previous[-(MAX_PEP_LEN-l):] + updated + after[:(MAX_PEP_LEN-l)]  

        elif 'frameshift_variant' in impact:
            
            # convert hg19 to hg38 point
            coord = lo_mapping[int(coord)]

            # step1 is to get start,end,before,after,in the context of hg38, mode does not matter, think about as following
            '''
            AC/A forward strand coord is the first nt
            I need two anchors start and end (reference coord), and between them I have before and after (in forward strand)
            '''

            before,after = genetic.split('/')
            start = int(coord) - 1
            end = start + len(before) + 1

            print(pos,change,enst,cat,impact,genetic,chrom,strand,coord)
            print(start,end,before,after)


            # step 2 is get n_codon and p_codon from cds file
            '''
            you will get start_cds_pos,end_cds_pos,n_codon,p_codon, for your two anchors, n_codon and p_codon are for start, easily verifiable, in the sense of translation
            '''
            cds = cds_all.loc[cds_all.index==enst,:]
            if cds.shape[0] == 1:
                s_ = cds.iloc[0]
                strand = s_['strand']
                cds_seq = s_['cds_seq']
                try:
                    coord = literal_eval(s_['coord'])
                    dic_info = literal_eval(s_['dic'])
                except:   # means None coord and dic, the incomplete or problematic CDS in ensembl
                    seq = None
                else:
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
                            if pos % 2 == 1:
                                n_exon = pos // 2 + 1
                            else:
                                seq = None
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
                            if pos % 2 == 1:
                                n_exon = pos // 2 + 1
                            else:
                                seq = None
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
            else:
                seq = None

            print(start_cds_pos,end_cds_pos,n_codon,p_codon)

            # step 3 is to get sequence using dual approaches
            if strand == '+':
                new_cds = cds_seq[:start_cds_pos] + after + cds_seq[end_cds_pos-1:]
            else:
                new_cds = cds_seq[:start_cds_pos] + str(Seq(after).reverse_complement()) + cds_seq[end_cds_pos-1:]
            new_cds_pep = str(Seq(new_cds).translate(to_stop=False))

            if '*' in new_cds_pep:
                variant_seq = new_cds_pep.split('*')[0]
                seq = variant_seq
            else:
                aa_seq = protein
                preceding = aa_seq[:n_codon-1]
                n_aa_needed = MAX_PEP_LEN - 1
                first_part = preceding[-n_aa_needed:]
                if strand == '+':
                    stretch = dict_fa[chrom][start-1-(p_codon-1):start] + after + dict_fa[chrom][end-1:end-1+1000]
                elif strand == '-':
                    stretch = dict_fa[chrom][start-1-1000:start] + after + dict_fa[chrom][end-1:end+(p_codon-1)]
                    stretch = str(Seq(stretch).reverse_complement())
                second_part = str(Seq(stretch).translate(to_stop=False)).split('*')[0]
                if len(second_part) > 0:
                    variant_seq = first_part + second_part
                    seq = variant_seq
                else:
                    seq = None
            print(seq)

        return seq

MAX_PEP_LEN = 15
MIN_PEP_LEN = 8
db = '/gpfs/data/yarmarkovichlab/public/ImmunoVerse/database'
canonical = pd.read_csv('/gpfs/data/yarmarkovichlab/HMF_vaccine/immunoverse_result/canonical.txt',sep='\t')
all_canonical_enst = set(canonical['ENST'].tolist())
ensg2strand = pd.Series(index=canonical['ENSG'].values,data=canonical['strand'].values).to_dict()
ensg2enst = pd.Series(index=canonical['ENSG'].values,data=canonical['ENST'].values).to_dict()

enst2protein = {}   # this contains all non-canonical isoforms if valid
with open(os.path.join(db,'Homo_sapiens.GRCh38.pep.all.fa'),'r') as in_handle:
    for title,seq in SimpleFastaParser(in_handle):
        pat1 = r'gene:(ENSG\d+)\.\d+'
        pat2 = r'transcript:(ENST\d+)\.\d+'
        ensg = re.search(pat1,title).group(1)
        enst = re.search(pat2,title).group(1)
        if (enst not in all_canonical_enst) and ('*' not in seq):
            enst2protein[enst] = seq

can_enst2protein = {}  
with open(os.path.join(db,'ensembl_protein.fasta'),'r') as in_handle:
    for title,seq_ in SimpleFastaParser(in_handle):
        ensg_,enst_,gs_ = title.split('|')
        can_enst2protein[enst_] = seq_

dict_fa = {}  
with open(os.path.join(db,'hg38.fa'),'r') as in_handle:
    for title,seq in SimpleFastaParser(in_handle):
        dict_fa[title] = seq

cds_all = pd.read_csv(os.path.join(db,'hg38_cds_final.txt'),sep='\t',index_col=0)
valid_chrom = set(['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM'])

mapping = {
    1:'chr1',
    2:'chr2',
    3:'chr3',
    4:'chr4',
    5:'chr5',
    6:'chr6',
    7:'chr7',
    8:'chr8',
    9:'chr9',
    10:'chr10',
    11:'chr11',
    12:'chr12',
    13:'chr13',
    14:'chr14',
    15:'chr15',
    16:'chr16',
    17:'chr17',
    18:'chr18',
    19:'chr19',
    20:'chr20',
    21:'chr21',
    22:'chr22',
    'X':'chrX',
    'Y':'chrY',
    'MT':'chrM',
    '1':'chr1',
    '2':'chr2',
    '3':'chr3',
    '4':'chr4',
    '5':'chr5',
    '6':'chr6',
    '7':'chr7',
    '8':'chr8',
    '9':'chr9',
    '10':'chr10',
    '11':'chr11',
    '12':'chr12',
    '13':'chr13',
    '14':'chr14',
    '15':'chr15',
    '16':'chr16',
    '17':'chr17',
    '18':'chr18',
    '19':'chr19',
    '20':'chr20',
    '21':'chr21',
    '22':'chr22',
}


# hg19 - hg38 fs
lo_mapping = {}
with open('fs_hg19.txt','r') as f1, open('fs_hg38.txt','r') as f2:
    for line1,line2 in zip(f1,f2):
        line1 = line1.rstrip('\n')
        _,_,pos1 = line1.split(' ')
        line2 = line2.rstrip('\n')
        _,_,pos2 = line2.split('\t')
        lo_mapping[int(pos1)] = int(pos2)

# start to process
vep_df = pd.read_csv('/gpfs/data/yarmarkovichlab/HMF_vaccine/CPCT02030461T_HC3NNBGX9_S7/CPCT02030461T_vep/CPCT02030461T.purple.somatic.vep.vcf',sep='\t',skiprows=97)
vep_df['#CHROM'] = vep_df['#CHROM'].map(mapping).values
vep_df = vep_df.loc[vep_df['#CHROM'].isin(valid_chrom),:]
cond = []
for item in vep_df['INFO']:
    if ('missense_variant') in item or ('inframe_deletion' in item) or ('inframe_insertion' in item) or ('frameshift_variant' in item):
        cond.append(True)
    else:
        cond.append(False)
vep_df = vep_df.loc[cond,:]
with open('wgs_variant.fasta','w') as f:
    for chrom,coord,info,ref,alt in zip(vep_df['#CHROM'],vep_df['POS'],vep_df['INFO'],vep_df['REF'],vep_df['ALT']):
        # modify, one key dict, some only has key if using dict(pair)
        dic = {}
        for pair in info.split(';'):
            tmp = pair.split('=')
            if len(tmp) == 2 and tmp[0] == 'CSQ':
                dic[tmp[0]] = tmp[1] 
        lis_csq = [item for item in dic['CSQ'].split(',') if ('missense_variant') in item or ('inframe_deletion' in item) or ('inframe_insertion' in item) or ('frameshift_variant' in item)]
        lis_enst = [item.split('|')[6] for item in lis_csq]
        lis_cat = [True if item in all_canonical_enst else False for item in lis_enst]
        csq = None
        for i,cat in enumerate(lis_cat):
            if cat:
                csq = lis_csq[i]
                break
        if csq is not None:
            cat = 'canonical'
            ensg = csq.split('|')[4]
            if ensg in ensg2strand.keys():
                gs = csq.split('|')[3]
                enst = csq.split('|')[6]
                pos = csq.split('|')[14]
                change = csq.split('|')[15]
                impact = csq.split('|')[1]
                genetic = csq.split('|')[16]
                strand = ensg2strand[ensg]
                if 'frameshift' in impact:
                    # print(chrom,int(coord)-1,int(coord))
                    genetic = ref + '/' + alt
                if 'missense' in impact and '-' in pos:
                    pos = int(pos.split('-')[1])
                    ref,alt = change.split('/')
                    ref,alt = ref[-1],alt[-1]
                    change = ref + '/' + alt
                seq = variant_effect_analyzer(pos,change,enst,cat,impact,genetic,chrom,strand,coord,enst2protein,can_enst2protein,dict_fa,cds_all)
                if seq is not None:
                    pos_change = change.split('/')[0] + str(pos) + change.split('/')[1]
                    f.write('>{}|{}|{}|{}|{}|{}|{}|{}|{}|wgs_called\n{}\n'.format(ensg,enst,gs,chrom,coord,strand,pos_change,impact,genetic,seq))