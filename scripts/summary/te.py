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



# analyze ORF2 and ORF1
df = pd.read_csv('final_all_ts_antigens.txt',sep='\t')
df_orf = df.loc[(df['typ'].isin(['self_translate_te'])) & (df['source'].str.contains('L1_ORF1')),:]
prioritized_peps = list(set(df_orf['pep'].values.tolist()))

orf1 = 'MGKKQNRKTGNSKTQSASPPPKERSSSPATEQSWMENDFDELREEGFRRSNYSELREDIQTKGKEVENFEKNLEECITRITNTEKCLKELMELKTKARELREECRSLRSRCDQLEERVSAMEDEMNEMKREGKFREKRIKRNEQSLQEIWDYVKRPNLRLIGVPESDVENGTKLENTLQDIIQENFPNLARQANVQIQEIQRTPQRYSSRRATPRHIIVRFTKVEMKEKMLRAAREKGRVTLKGKPIRLTADLSAETLQARREWGPIFNILKEKNFQPRISYPAKLSFISEGEIKYFIDKQMLRDFVTTRPALKELLKEALNMERNNRYQPLQNHAKM'
orf2 = 'MTGSNSHITILTLNVNGLNSPIKRHRLASWIKSQDPSVCCIQETHLTCRDTHRLKIKGWRKIYQANGKQKKAGVAILVSDKTDFKPTKIKRDKEGHYIMVKGSIQQEELTILNIYAPNTGAPRFIKQVLSDLQRDLDSHTLIMGDFNTPLSILDRSTRQKVNKDTQELNSALHQTDLIDIYRTLHPKSTEYTFFSAPHHTYSKIDHIVGSKALLSKCKRTEIITNYLSDHSAIKLELRIKNLTQSRSTTWKLNNLLLNDYWVHNEMKAEIKMFFETNENKDTTYQNLWDAFKAVCRGKFIALNAYKRKQERSKIDTLTSQLKELEKQEQTHSKASRRQEITKIRAELKEIETQKTLQKINESRSWFFERINKIDRPLARLIKKKREKNQIDTIKNDKGDITTDPTEIQTTIREYYKHLYANKLENLEEMDTFLDTYTLPRLNQEEVESLNRPITGSEIVAIINSLPTKKSPGPDGFTAEFYQRYKEELVPFLLKLFQSIEKEGILPNSFYEASIILIPKPGRDTTKKENFRPISLMNIDAKILNKILANRIQQHIKKLIHHDQVGFIPGMQGWFNIRKSINVIQHINRAKDKNHVIISIDAEKAFDKIQQPFMLKTLNKLGIDGMYLKIIRAIYDKPTANIILNGQKLEAFPLKTGTRQGCPLSPLLFNIVLEVLARAIRQEKEIKGIQLGKEEVKLSLFADDMIVYLENPIVSAQNLLKLISNFSKVSGYKINVQKSQAFLYNNNRQTESQIMGELPFTIASKRIKYLGIQLTRDVKDLFKENYKPLLKEIKEDTNKWKNIPCSWVGRINIVKMAILPKVIYRFNAIPIKLPMTFFTELEKTTLKFIWNQKRARIAKSILSQKNKAGGITLPDFKLYYKATVTKTAWYWYQNRDIDQWNRTEPSEIMPHIYNYLIFDKPEKNKQWGKDSLLNKWCWENWLAICRKLKLDPFLTPYTKINSRWIKDLNVKPKTIKTLEENLGITIQDIGVGKDFMSKTPKAMATKDKIDKWDLIKLKSFCTAKETTIRVNRQPTTWEKIFATYSSDKGLISRIYNELKQIYKKKTNNPIKKWAKDMNRHFSKEDIYAAKKHMKKCSSSLAIREMQIKTTMRYHLTPVRMAIIKKSGNNRCWRGCGEIGTLVHCWWDCKLVQPLWKSVWRFLRDLELEIPFDPAIPLLGIYPKDYKSCCYKDTCTRMFIAALFTIAKTWNQPNCPTMIDWIKKMWHIYTMEYYAAIKNDEFISFVGTWMKLETIILSKLSQEQKTKHRIFSLIGGN'
data = []
for pep,sub_df in df_orf.groupby(by='pep'):
    cancers_ = sub_df['cancer'].values.tolist()
    pos = orf1.index(pep)
    hla = sub_df['additional_query'].iloc[0]
    data.append((pep,pos,len(cancers_),','.join(cancers_),hla))
final_orf_df = pd.DataFrame.from_records(data,columns=['pep','pos','rec','cancers','hla']).sort_values(by='pos')
final_orf_df.to_csv('orf1_stat.txt',sep='\t',index=None)


safety_screen_df = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/code/hla_ligand_atlas_now_0.05_tesorai.txt',sep='\t')
all_tissues = ['Adrenal gland', 'Aorta', 'Bladder', 'Bone marrow', 'Brain', 'Cerebellum', 'Colon', 'Esophagus', 'Gallbladder', 'Heart', 'Kidney', 'Liver', 
                'Lung', 'Lymph node', 'Mamma', 'Muscle', 'Myelon', 'Ovary', 'Pancreas', 'Prostate', 'Skin', 'Small intestine', 'Spleen', 'Stomach', 'Testis', 
                'Thymus', 'Thyroid', 'Tongue', 'Trachea', 'Uterus']
final = df_orf
store_data = []
for k in prioritized_peps:
    final_p = final.loc[final['pep']==k,:]
    all_occur = final_p['cancer'].values.tolist()
    tmp = []
    for c in cancers:
        if c in all_occur:
            sub = final_p.loc[final_p['cancer']==c,:]
            all_intensity = []
            for item in sub['detailed_intensity']:
                all_intensity.extend(literal_eval(item))
            med_intensity = np.median(all_intensity)
            tmp.append(med_intensity)
        else:
            tmp.append(0)

    all_tissues = np.array(all_tissues)
    tmp_normal = np.full(shape=len(all_tissues),fill_value=0.0)
    tmp_normal_df = safety_screen_df.loc[safety_screen_df['peptide']==k,:]
    for t,sub_df in tmp_normal_df.groupby(by='tissue'):
        med_intensity = np.median(sub_df['percentile'].values)
        indices = np.where(all_tissues == t)[0]
        tmp_normal[indices[0]] = med_intensity
    
    tmp = tmp + tmp_normal.tolist()
    store_data.append(tmp)

df = pd.DataFrame(data=store_data,index=prioritized_peps,columns=cancers+list(all_tissues))
ori_array = [tuple(['cancer']*21+['normal']*30),tuple(df.columns.tolist())]
mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
df.columns = mi
ori_array = [tuple(df.index.tolist())]
mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
df.index = mi
df.to_csv('peptide_view_orf1.txt',sep='\t')



# original main
data = []
for c in cancers:
    final_path = os.path.join(root_atlas_dir,c,'antigen','fdr','final_enhanced.txt')
    final = pd.read_csv(final_path,sep='\t')
    cond = [False if ('[]' in item) and ('(\'HLA-' not in item) else True for item in final['presented_by_each_sample_hla']]
    final = final.loc[cond,:]
    final = final.loc[final['typ'].isin(['TE_chimeric_transcript','ERV']),:]
    data.append(final)
final = pd.concat(data,axis=0,keys=cancers).reset_index(level=-2).rename(columns={'level_0':'cancer'})
final.to_csv('te_all_antigens.txt',sep='\t',index=None)


df = pd.read_csv('./stats/final_all_ts_antigens.txt',sep='\t')
final2 = df.loc[df['typ'].isin(['self_translate_te','TE_chimeric_transcript']),:]
final2.to_csv('all_te_final.txt',sep='\t',index=None)

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
    data.extend(tid2medians.get(tid,[0]*31))
    all_data.append(data)

df = pd.DataFrame.from_records(all_data,columns=cancers+all_tissues,index=all_tids)

ori_array = [tuple(['cancer']*21+['normal']*31),tuple(df.columns.tolist())]
mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
df.columns = mi

tmp1 = tuple(df.index.tolist()) # for later
tmp2 = tuple([tid2full[item].split('|')[4] for item in df.index]) # for later, because multiindex will mess up things
tmp3 = tuple([tid2cid[item] for item in df.index])
tmp4 = tuple([tid2full[item].split('|')[-1] for item in df.index])

ori_array = [tuple(df.index.tolist()),
             tuple([tid2full[item].split('|')[4] for item in df.index]),
             tuple([tid2cid[item] for item in df.index]),
             tuple([tid2full[item].split('|')[-1] for item in df.index])]
mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
df.index = mi
df.to_csv('te_holy.txt',sep='\t')

ts_te = []
ts_tid2cid = {}
ts_tid2lfc = {}
for item1,item2,item3,item4 in zip(tmp1,tmp2,tmp3,tmp4):
    if float(item2) > 5:
        ts_te.append(item1)
        ts_tid2cid[item1] = item3
        ts_tid2lfc[item1] = item2
final_ts = final_ts.loc[[True if item.split('|')[0] in ts_te else False for item in final_ts['source']],:]
final_ts['cid'] = [ts_tid2cid[item.split('|')[0]] for item in final_ts['source']]
final_ts['lfc'] = [ts_tid2lfc[item.split('|')[0]] for item in final_ts['source']]
final_ts.to_csv('ts_te_antigen.txt',sep='\t',index=None)







