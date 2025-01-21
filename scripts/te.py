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
te_gtf_path = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/splicing/GRCh38_Ensembl_rmsk_TE.gtf'
hg38_normal_dir = '/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/NeoVerse/GTEx/selected/hg38_telocal_intron'


def process_te_gtf():
    dic = {}
    te_df = pd.read_csv(te_gtf_path,sep='\t',header=None)
    te_df.columns = ['chrom','source','feature','start','end','score','strand','phase','attribute']
    for item in te_df['attribute']:
        gid,tid,fid,cid,_ = item.split(';')
        gid = gid.split(' "')[1].rstrip('"')
        tid = tid.split(' "')[1].rstrip('"')
        fid = fid.split(' "')[1].rstrip('"')
        cid = cid.split(' "')[1].rstrip('"')
        dic[tid] = cid
    return dic

def process_te_gtex():
    normal_erv = pd.read_csv(os.path.join(hg38_normal_dir,'normal_erv.txt'),sep='\t',index_col=0)
    normal_erv_aux = pd.read_csv(os.path.join(hg38_normal_dir,'normal_erv_aux_df.txt'),sep='\t',index_col=0)
    n_data = normal_erv.values / normal_erv_aux['total_count'].values.reshape(1,-1) * 1e6
    normal_erv_cpm = pd.DataFrame(data=n_data,index=normal_erv.index,columns=normal_erv.columns)
    normal_erv_cpm = normal_erv_cpm.loc[[True if ':' in item else False for item in normal_erv_cpm.index],:].T
    normal_erv_cpm['tissue'] = [item.split(',')[1] for item in normal_erv_cpm.index]
    final_data = np.empty(shape=(normal_erv_cpm.shape[1]-1,31),dtype=np.float64)
    tid2medians = {}
    all_tissues = []
    for i,(t,sub_df) in enumerate(normal_erv_cpm.groupby(by='tissue')):
        sub_df = sub_df.iloc[:,:-1]
        all_tissues.append(t)
        t_meds = np.median(sub_df.values,axis=0).reshape(-1)
        final_data[:,i] = t_meds
    normal_erv_cpm.drop(columns='tissue',inplace=True)
    for i,tid_full in enumerate(normal_erv_cpm.columns):
        tid2medians[tid_full.split(':')[0]] = final_data[i,:].reshape(-1).tolist()
    return tid2medians,all_tissues
    

def process_tumor_te():
    dic = {}
    for c in tqdm(cancers):
        erv_path = os.path.join(root_atlas_dir,c,'ERV.txt')
        erv = pd.read_csv(erv_path,sep='\t',index_col=0)
        erv = erv.loc[[True if ':' in item else False for item in erv.index],:]
        erv.index = [item.split(':')[0] for item in erv.index]
        dic[c] = erv['median_tumor'].to_dict()
    return dic



data = []
for c in cancers:
    final_path = os.path.join(root_atlas_dir,c,'antigen','fdr','final_enhanced.txt')
    final = pd.read_csv(final_path,sep='\t')
    cond = [False if '[]' in item else True for item in final['presented_by_each_sample_hla']]
    final = final.loc[cond,:]
    final = final.loc[final['typ'].isin(['TE_chimeric_transcript','ERV']),:]
    data.append(final)
final = pd.concat(data,axis=0,keys=cancers).reset_index(level=-2).rename(columns={'level_0':'cancer'})


# overreprensentaion in aml
lis = []
with open('all_annotations.txt','r') as f:
    for line in f:
        lis.append(line.rstrip('\n').split('_')[0])
        print(line.rstrip('\n').split('_')[0])
vc = pd.Series(lis).value_counts()
vc.to_csv('aml_te_or.txt',sep='\t')
sys.exit('stop')

# ORF2
orf2 = 'MTGSNSHITILTLNVNGLNSPIKRHRLASWIKSQDPSVCCIQETHLTCRDTHRLKIKGWRKIYQANGKQKKAGVAILVSDKTDFKPTKIKRDKEGHYIMVKGSIQQEELTILNIYAPNTGAPRFIKQVLSDLQRDLDSHTLIMGDFNTPLSILDRSTRQKVNKDTQELNSALHQTDLIDIYRTLHPKSTEYTFFSAPHHTYSKIDHIVGSKALLSKCKRTEIITNYLSDHSAIKLELRIKNLTQSRSTTWKLNNLLLNDYWVHNEMKAEIKMFFETNENKDTTYQNLWDAFKAVCRGKFIALNAYKRKQERSKIDTLTSQLKELEKQEQTHSKASRRQEITKIRAELKEIETQKTLQKINESRSWFFERINKIDRPLARLIKKKREKNQIDTIKNDKGDITTDPTEIQTTIREYYKHLYANKLENLEEMDTFLDTYTLPRLNQEEVESLNRPITGSEIVAIINSLPTKKSPGPDGFTAEFYQRYKEELVPFLLKLFQSIEKEGILPNSFYEASIILIPKPGRDTTKKENFRPISLMNIDAKILNKILANRIQQHIKKLIHHDQVGFIPGMQGWFNIRKSINVIQHINRAKDKNHVIISIDAEKAFDKIQQPFMLKTLNKLGIDGMYLKIIRAIYDKPTANIILNGQKLEAFPLKTGTRQGCPLSPLLFNIVLEVLARAIRQEKEIKGIQLGKEEVKLSLFADDMIVYLENPIVSAQNLLKLISNFSKVSGYKINVQKSQAFLYNNNRQTESQIMGELPFTIASKRIKYLGIQLTRDVKDLFKENYKPLLKEIKEDTNKWKNIPCSWVGRINIVKMAILPKVIYRFNAIPIKLPMTFFTELEKTTLKFIWNQKRARIAKSILSQKNKAGGITLPDFKLYYKATVTKTAWYWYQNRDIDQWNRTEPSEIMPHIYNYLIFDKPEKNKQWGKDSLLNKWCWENWLAICRKLKLDPFLTPYTKINSRWIKDLNVKPKTIKTLEENLGITIQDIGVGKDFMSKTPKAMATKDKIDKWDLIKLKSFCTAKETTIRVNRQPTTWEKIFATYSSDKGLISRIYNELKQIYKKKTNNPIKKWAKDMNRHFSKEDIYAAKKHMKKCSSSLAIREMQIKTTMRYHLTPVRMAIIKKSGNNRCWRGCGEIGTLVHCWWDCKLVQPLWKSVWRFLRDLELEIPFDPAIPLLGIYPKDYKSCCYKDTCTRMFIAALFTIAKTWNQPNCPTMIDWIKKMWHIYTMEYYAAIKNDEFISFVGTWMKLETIILSKLSQEQKTKHRIFSLIGGN'
cond = []
for item in final['source']:
    if 'L1_ORF2' in item:
        cond.append(True)
    else:
        cond.append(False)
final_orf2 = final.loc[cond,:]
final_orf2.to_csv('final_orf2.txt',sep='\t',index=None)
sys.exit('stop')


tid2cid = process_te_gtf()


# find tumor specific one
# tid2medians,all_tissues = process_te_gtex()
# dic = process_tumor_te()

# with open('tid2medians.p','wb') as f:
#     pickle.dump(tid2medians,f)
# with open('all_tissues.p','wb') as f:
#     pickle.dump(all_tissues,f)
# with open('dic.p','wb') as f:
#     pickle.dump(dic,f)

with open('tid2medians.p','rb') as f:
    tid2medians = pickle.load(f)
with open('all_tissues.p','rb') as f:
    all_tissues = pickle.load(f)
with open('dic.p','rb') as f:
    dic = pickle.load(f)

final_ts = final.loc[(final['typ']=='ERV') & (final['unique']) & (~final['source'].str.contains(';')),:]
all_tids = []
tid2full = {}
for item in final_ts['source']:
    tid = item.split('|')[0]
    all_tids.append(tid)
    tid2full.setdefault(tid,[]).append(item)
all_tids = list(set(all_tids))

tmp = {}
for k,v in tid2full.items():
    if len(v) == 1:
        tmp[k] = v[0]
    else:
        lfcs = []
        for item in v:
            lfcs.append(float(item.split('|')[4]))
        tmp[k] = v[np.argmax(lfcs)]
tid2full = tmp


tid2tumors = {}
for tid in all_tids:
    data = []
    for k,v in dic.items():
        data.append(v[tid])
    tid2tumors[tid] = data

all_data = []
for tid in all_tids:
    data = []
    data.extend(tid2tumors[tid])
    data.extend(tid2medians[tid])
    all_data.append(data)

df = pd.DataFrame.from_records(all_data,columns=cancers+all_tissues,index=all_tids)

ori_array = [tuple(['cancer']*21+['normal']*31),tuple(df.columns.tolist())]
mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
df.columns = mi

ori_array = [tuple(df.index.tolist()),
             tuple([tid2full[item].split('|')[4] for item in df.index]),
             tuple([tid2cid[item] for item in df.index]),
             tuple([tid2full[item].split('|')[-1] for item in df.index])]
mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
df.index = mi
df.to_csv('te_holy.txt',sep='\t')

tmp1 = df.index.levels[1]
tmp2 = df.index.levels[0]
ts_te = []
for item1,item2 in zip(tmp1,tmp2):
    if float(item1) > 5:
        ts_te.append(item2)
final_ts = final_ts.loc[[True if item.split('|')[0] in ts_te else False for item in final_ts['source']],:]
final_ts.to_csv('ts_te_antigen.txt',sep='\t',index=None)
sys.exit('stop')


# category
col = []
for item1,item2 in zip(final['typ'],final['source']):

    if 'L1_ORF1' in item2:
        col.append('canonical_ORF1')
    elif 'L1_ORF2' in item2:
        col.append('canonical_ORF2')
    elif 'ENSG' in item2 and 'ERVK' in item2:
        col.append('canonical_ERVK')
    else:

        if item1 == 'TE_chimeric_transcript':
            if 'LINE' in item2:
                col.append('chimera_LINE')
            elif 'LTR' in item2:
                col.append('chimera_LTR')
            elif 'SINE' in item2:
                col.append('chimera_SINE')
            elif 'SVA' in item2:
                col.append('chimera_Retroposon_SVA')
            elif 'DNA' in item2:
                col.append('chimera_DNA_transposon')
            elif 'Satellite' in item2:
                col.append('chimera_Satellite')
            else:
                col.append('chimera_unknown')
        elif item1 == 'ERV':
            uid = item2.split('|')[0]
            cid = tid2cid.get(uid,'unknown')

            if cid != 'unknown':
                if cid == 'Retroposon':
                    cid = 'Retroposon_SVA'
                col.append('self_translate_{}'.format(cid))
            else:
                if 'LINE' in item2:
                    col.append('both_LINE')
                elif 'LTR' in item2:
                    col.append('both_LTR')
                elif 'SINE' in item2:
                    col.append('both_SINE')
                elif 'SVA' in item2:
                    col.append('both_Retroposon_SVA')
                elif 'DNA' in item2:
                    col.append('both_DNA_transposon')
                elif 'Satellite' in item2:
                    col.append('both_Satellite')
                else:
                    col.append('both_unknown')
final['class'] = col
final = final.loc[final['class']!='both_unknown',:]
# final.to_csv('te_all_antigens.txt',sep='\t',index=None)

data = []
all_types = []
for s,sub_df in final.groupby(by='class'):
    all_types.append(s)
    vc = sub_df['cancer'].value_counts().to_dict()
    row = []
    for c in cancers:
        row.append(vc.get(c,0))
    data.append(row)
df = pd.DataFrame.from_records(data=data,index=all_types,columns=cancers)
df.to_csv('te_class_df.txt',sep='\t')




