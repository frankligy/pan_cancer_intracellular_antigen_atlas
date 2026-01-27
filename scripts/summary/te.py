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



# # analyze ORF2 and ORF1
# df = pd.read_csv('./stats/final_all_ts_antigens.txt',sep='\t')
# df_orf = df.loc[df['source'].str.contains('L1_ORF2'),:]
# prioritized_peps = list(set(df_orf['pep'].values.tolist()))

# orf1 = 'MGKKQNRKTGNSKTQSASPPPKERSSSPATEQSWMENDFDELREEGFRRSNYSELREDIQTKGKEVENFEKNLEECITRITNTEKCLKELMELKTKARELREECRSLRSRCDQLEERVSAMEDEMNEMKREGKFREKRIKRNEQSLQEIWDYVKRPNLRLIGVPESDVENGTKLENTLQDIIQENFPNLARQANVQIQEIQRTPQRYSSRRATPRHIIVRFTKVEMKEKMLRAAREKGRVTLKGKPIRLTADLSAETLQARREWGPIFNILKEKNFQPRISYPAKLSFISEGEIKYFIDKQMLRDFVTTRPALKELLKEALNMERNNRYQPLQNHAKM'
# orf2 = 'MTGSNSHITILTLNVNGLNSPIKRHRLASWIKSQDPSVCCIQETHLTCRDTHRLKIKGWRKIYQANGKQKKAGVAILVSDKTDFKPTKIKRDKEGHYIMVKGSIQQEELTILNIYAPNTGAPRFIKQVLSDLQRDLDSHTLIMGDFNTPLSILDRSTRQKVNKDTQELNSALHQTDLIDIYRTLHPKSTEYTFFSAPHHTYSKIDHIVGSKALLSKCKRTEIITNYLSDHSAIKLELRIKNLTQSRSTTWKLNNLLLNDYWVHNEMKAEIKMFFETNENKDTTYQNLWDAFKAVCRGKFIALNAYKRKQERSKIDTLTSQLKELEKQEQTHSKASRRQEITKIRAELKEIETQKTLQKINESRSWFFERINKIDRPLARLIKKKREKNQIDTIKNDKGDITTDPTEIQTTIREYYKHLYANKLENLEEMDTFLDTYTLPRLNQEEVESLNRPITGSEIVAIINSLPTKKSPGPDGFTAEFYQRYKEELVPFLLKLFQSIEKEGILPNSFYEASIILIPKPGRDTTKKENFRPISLMNIDAKILNKILANRIQQHIKKLIHHDQVGFIPGMQGWFNIRKSINVIQHINRAKDKNHVIISIDAEKAFDKIQQPFMLKTLNKLGIDGMYLKIIRAIYDKPTANIILNGQKLEAFPLKTGTRQGCPLSPLLFNIVLEVLARAIRQEKEIKGIQLGKEEVKLSLFADDMIVYLENPIVSAQNLLKLISNFSKVSGYKINVQKSQAFLYNNNRQTESQIMGELPFTIASKRIKYLGIQLTRDVKDLFKENYKPLLKEIKEDTNKWKNIPCSWVGRINIVKMAILPKVIYRFNAIPIKLPMTFFTELEKTTLKFIWNQKRARIAKSILSQKNKAGGITLPDFKLYYKATVTKTAWYWYQNRDIDQWNRTEPSEIMPHIYNYLIFDKPEKNKQWGKDSLLNKWCWENWLAICRKLKLDPFLTPYTKINSRWIKDLNVKPKTIKTLEENLGITIQDIGVGKDFMSKTPKAMATKDKIDKWDLIKLKSFCTAKETTIRVNRQPTTWEKIFATYSSDKGLISRIYNELKQIYKKKTNNPIKKWAKDMNRHFSKEDIYAAKKHMKKCSSSLAIREMQIKTTMRYHLTPVRMAIIKKSGNNRCWRGCGEIGTLVHCWWDCKLVQPLWKSVWRFLRDLELEIPFDPAIPLLGIYPKDYKSCCYKDTCTRMFIAALFTIAKTWNQPNCPTMIDWIKKMWHIYTMEYYAAIKNDEFISFVGTWMKLETIILSKLSQEQKTKHRIFSLIGGN'
# data = []
# for pep,sub_df in df_orf.groupby(by='pep'):
#     cancers_ = sub_df['cancer'].values.tolist()
#     pos = orf2.index(pep)
#     hla = sub_df['additional_query'].iloc[0]
#     data.append((pep,pos,len(cancers_),','.join(cancers_),hla))
# final_orf_df = pd.DataFrame.from_records(data,columns=['pep','pos','rec','cancers','hla']).sort_values(by='pos')
# final_orf_df.to_csv('orf2_stat.txt',sep='\t',index=None)


# safety_screen_df = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/code/hla_ligand_atlas_now_0.05_tesorai.txt',sep='\t')
# all_tissues = ['Adrenal gland', 'Aorta', 'Bladder', 'Bone marrow', 'Brain', 'Cerebellum', 'Colon', 'Esophagus', 'Gallbladder', 'Heart', 'Kidney', 'Liver', 
#                 'Lung', 'Lymph node', 'Mamma', 'Muscle', 'Myelon', 'Ovary', 'Pancreas', 'Prostate', 'Skin', 'Small intestine', 'Spleen', 'Stomach', 'Testis', 
#                 'Thymus', 'Thyroid', 'Tongue', 'Trachea', 'Uterus']
# final = df_orf
# store_data = []
# for k in prioritized_peps:
#     final_p = final.loc[final['pep']==k,:]
#     all_occur = final_p['cancer'].values.tolist()
#     tmp = []
#     for c in cancers:
#         if c in all_occur:
#             sub = final_p.loc[final_p['cancer']==c,:]
#             all_intensity = []
#             for item in sub['detailed_intensity']:
#                 all_intensity.extend(literal_eval(item))
#             med_intensity = np.median(all_intensity)
#             tmp.append(med_intensity)
#         else:
#             tmp.append(0)

#     all_tissues = np.array(all_tissues)
#     tmp_normal = np.full(shape=len(all_tissues),fill_value=0.0)
#     tmp_normal_df = safety_screen_df.loc[safety_screen_df['peptide']==k,:]
#     for t,sub_df in tmp_normal_df.groupby(by='tissue'):
#         med_intensity = np.median(sub_df['percentile'].values)
#         indices = np.where(all_tissues == t)[0]
#         tmp_normal[indices[0]] = med_intensity
    
#     tmp = tmp + tmp_normal.tolist()
#     store_data.append(tmp)

# df = pd.DataFrame(data=store_data,index=prioritized_peps,columns=cancers+list(all_tissues))
# ori_array = [tuple(['cancer']*21+['normal']*30),tuple(df.columns.tolist())]
# mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
# df.columns = mi
# ori_array = [tuple(df.index.tolist())]
# mi = pd.MultiIndex.from_arrays(arrays=ori_array,sortorder=0)
# df.index = mi
# df.to_csv('peptide_view_orf2.txt',sep='\t')


# main
df = pd.read_csv('./stats/final_all_ts_antigens.txt',sep='\t')
final = df.loc[df['typ'].isin(['self_translate_te','TE_chimeric_transcript']),:]

tid2cid = process_te_gtf()
col1 = []
col2 = []
col3 = []
for items in final['source']:
    tmp1 = []
    tmp2 = []
    for item in items.split(';'):
        if 'TE_info' in item:
            tid = item.split('TE_info:')[1].split(',')[1]
            tmp1.append(tid)
            tmp2.append(tid2cid[tid])
        elif '_dup' in item or item.endswith('sense'):
            tid = item.split('|')[0]
            tmp1.append(tid)
            tmp2.append(tid2cid[tid])
        else:
            continue
    col1.append(','.join(tmp1))
    col2.append(','.join(tmp2))
    if len(set(tmp2)) == 1:
        col3.append(tmp2[0])
    else:
        col3.append('ambiguous')
final['tid'] = col1
final['cid'] = col2
final['verdict'] = col3
final.to_csv('all_te_final.txt',sep='\t',index=None)



col1 = []
col2 = []
col3 = []
for cid, sub_df in final.groupby(by='verdict'):
    n_auto = sub_df.loc[sub_df['typ']=='self_translate_te',:].shape[0]
    n_chi = sub_df.loc[sub_df['typ']=='TE_chimeric_transcript',:].shape[0]
    col1.append(n_auto)
    col2.append(n_chi)
    col3.append(cid)
df = pd.DataFrame(index=col3,data={'autonomous':col1,'chimera':col2})


df.T.plot.pie(subplots=True,figsize=(20,6))
plt.savefig('te_auto_chi_pie.pdf',bbox_inches='tight')
plt.close()

df.to_csv('te_auto_chi.txt',sep='\t')







