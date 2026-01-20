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

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

from Bio.Blast import NCBIWWW, NCBIXML


def run_blastn(seq):
    hits = []
    seq = 'CTCTGCTTTAGCTGATTTCCCTTCTGCTTTGGCTGATTCTCCCTCCAGTTTGGCTGATTTCACCTCCGCTTTGGC'
    result_handle = NCBIWWW.qblast('blastn','nt',Seq(seq))
    with open('results.xml', 'w') as f: 
        blast_results = result_handle.read() 
        f.write(blast_results)
    for record in NCBIXML.parse(open('results.xml')):
        if record.alignments:
            for align in record.alignments:
                for hsp in align.hsps:
                    if hsp.expect < 1e-20:
                        # print(align.title)
                        # print(hsp.query)
                        # print(hsp.match)
                        # print(hsp.sbjct)
                        hits.append(align.title)
    return hits

dic = {}
with open('aligned_seq_to_ref_genome_and_cds.fasta','r') as in_handle:
    for title,seq in SimpleFastaParser(in_handle):
        dic[title] = seq

with open('check.txt','w') as f:
    for k,v in tqdm(dic.items(),total=len(dic)):
        hits = run_blastn(v)
        f.write('{}\t{}\t{}\n'.format(k,v,'@'.join(hits)))



