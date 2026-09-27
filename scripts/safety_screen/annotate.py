#!/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/bin/python3.7

import pandas as pd
import numpy as np
import sys,os
import subprocess
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import re

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

actual_final = []

# 30 tissues
# website miss stomach
# compared to the peptide.py which is the folder name aman defined, here needs to change a few to match the readme

tissues = [
    'AdrenalGland',
    'Aorta',
    'Bladder',
    'BoneMarrow',
    'Brain',
    'Cerebellum',
    'Colon',
    'Esophagus',
    'Gallbladder',
    'Heart',
    'Kidney',
    'Liver',
    'Lung',
    'LymphNode',
    'Mamma',
    'Muscle',
    'Myelon',
    'Ovary',
    'Pancreas',
    'Prostate',
    'Skin',
    'SmallIntestine',
    'Spleen',
    'Stomach',
    'Testis',
    'Thymus',
    'Thyroid',
    'Tongue',
    'Trachea',
    'Uterus',
]

df = pd.read_csv('all_raw_files.txt',sep='\t',header=None)
df.columns = ['id','file','ftp','type','other']

# do have whiltesplace remember
anno = {
    'DN02':'A*11:01; A*68:01; B*15:01; B*35:03; C*03:03; C*04:01',
    'DN03':'A*01:01; A*11:01; B*15:01; B*35:01; C*03:03; C*04:01',
    'DN04':'A*02:01; A*23:01; B*27:05; B*50:01; C*02:02; C*06:02',
    'DN05':'A*01:01; A*11:01; B*07:02; B*49:01; C*07:01; C*07:02',
    'DN06':'A*03:01; A*68:02; B*07:02; B*14:02; C*07:02; C*08:02',
    'DN08':'A*32:01; A*68:01; B*15:01; B*44:02; C*03:03; C*07:04',
    'DN09':'A*24:02; A*30:01; B*13:02; B*35:08; C*04:01; C*06:02',
    'DN11':'A*01:01; A*69:01; B*37:01; B*49:01; C*06:02; C*07:01',
    'DN12':'A*02:01; A*11:01; B*15:01; B*35:01; C*03:04; C*04:01',
    'DN13':'A*02:05; A*11:01; B*40:02; B*58:01; C*02:02; C*07:01',
    'DN14':'A*02:01; A*68:02; B*14:02; B*27:05; C*02:02; C*08:02',
    'DN15':'A*01:01; A*02:01; B*08:01; B*44:02; C*07:01; C*07:04',
    'DN16':'A*01:01; A*24:02; B*08:01; B*41:01; C*07:01; C*17:01',
    'DN17':'A*03:01; A*24:02; B*35:03; B*45:01; C*04:01; C*16:01',
    'DN278':'A*02:01; B*44:02; C*05:01',
    'DN281':'A*11:01; A*26:01; B*08:01; B*35:01; C*07:02; C*04:01',
    'DN1':'A*03:01; A*29:02; B*07:02; B*44:03; C*07:02; C*16:01',
    'DN3':'A*24:02; A*25:01; B*18:01; B*41:01; C*12:03; C*17:01',
    'DN4':'A*02:01; A*26:08; B*15:01; B*44:02; C*03:04; C*05:01',
    'DN5':'A*01:01; A*03:01; B*07:06; B*07:02; C*07:02; C*15:05',
    'DN6':'A*01:01; A*25:01; B*13:02; B*39:01; C*06:02; C*12:03'
}

pat = re.compile(r'-(DN\d+)_')

data = []
for t in tissues:
    sub = df.loc[df['file'].str.contains(t),:]
    for f in sub['file']:
        uid = re.search(pat,f).group(1)
        hla = anno[uid]
        data.append((t,f,uid,hla))
final = pd.DataFrame.from_records(data,columns=['tissue','file','uid','hla'])
actual_final.append(final)


'''hepatocytes'''
cell_types = [
    "hepatocytes",
    "hepatocytes",
    "hepatocytes",
    "hepatocytes",
    "hepatocytes",
    "hepatocytes",
    "hepatocytes",
    "hepatocytes",
    "hepatocytes",
    "hepatocytes",
    "hepatocytes",
    "hepatocytes",
]

raw_files = [
    "2064_L_MdB072_7mei19_donor_MdB026_OT.raw",
    "2064_L_MdB072_7mei19_donor_MdB026_OT_2.raw",
    "2047_L_MdB_6jan2020_MdB024.raw",
    "2047_L_MdB_6jan2020_MdB024_20200106173725.raw",
    "2062_L_MdB070_3mei19_donor_MdB028_OT.raw",
    "2062_L_MdB070_3mei19_donor_MdB028_2_OT.raw",
    "2047_L_MdB_6jan2020_MdB031.raw",
    "2047_L_MdB_6jan2020_MdB031_20200107003149.raw",
    "2064_L_MdB072_7mei19_donor_MdB051_OT.raw",
    "2064_L_MdB072_7mei19_donor_MdB051_OT_2.raw",
    "2062_L_MdB070_3mei19_donor_MdB052_OT.raw",
    "2062_L_MdB070_3mei19_donor_MdB052_2_OT.raw",
]

sample_ids = [
    "donor01",
    "donor01",
    "donor02",
    "donor02",
    "donor03",
    "donor03",
    "donor04",
    "donor04",
    "donor05",
    "donor05",
    "donor06",
    "donor06",
]

hla_types = [
    "A*01:01; A*03:01; B*07:02; B*15:01; C*04:01; C*07:02",
    "A*01:01; A*03:01; B*07:02; B*15:01; C*04:01; C*07:02",
    "A*02:02; A*33:01; B*14:02; B*53:01; C*04:01; C*08:02",
    "A*02:02; A*33:01; B*14:02; B*53:01; C*04:01; C*08:02",
    "A*02:01; A*03:01; B*07:02; B*44:02; C*07:02; C*07:04",
    "A*02:01; A*03:01; B*07:02; B*44:02; C*07:02; C*07:04",
    "A*02:01; A*30:01; B*35:01; B*42:01; C*16:01; C*17:01",
    "A*02:01; A*30:01; B*35:01; B*42:01; C*16:01; C*17:01",
    "A*01:01; A*24:02; B*13:02; B*27:05; C*01:02; C*06:02",
    "A*01:01; A*24:02; B*13:02; B*27:05; C*01:02; C*06:02",
    "A*24:02; A*26:01; B*15:01; B*38:01; C*03:03; C*12:03",
    "A*24:02; A*26:01; B*15:01; B*38:01; C*03:03; C*12:03",
]

hepatocytes_df = pd.DataFrame({
    "tissue": cell_types,
    "file": raw_files,
    "uid": sample_ids,
    "hla": hla_types,
})

actual_final.append(hepatocytes_df)

'''beta_cell'''
cell_types = [
    "beta_cell",
    "beta_cell",
    "beta_cell",
    "beta_cell",
    "beta_cell",
    "beta_cell",
    "beta_cell",
    "beta_cell",
    "beta_cell",
    "beta_cell",
    "beta_cell",
    "beta_cell",
    "beta_cell",
    "beta_cell",
    "beta_cell",
    "beta_cell",
]

raw_files = [
    "FL1100_IPP132_R1_Basal_even.raw",
    "FL1100_IPP132_R1_Basal_odd.raw",
    "FL1100_IPP132_R1_INFa_even.raw",
    "FL1100_IPP132_R1_INFa_odd.raw",
    "FL1100_IPP132_R2_Basal_even.raw",
    "FL1100_IPP132_R2_Basal_odd.raw",
    "FL1100_IPP132_R2_INFa_even.raw",
    "FL1100_IPP132_R2_INFa_odd.raw",
    "FL1100_IPP132_R3_Basal_even.raw",
    "FL1100_IPP132_R3_Basal_odd.raw",
    "FL1100_IPP132_R3_INFa_even.raw",
    "FL1100_IPP132_R3_INFa_odd.raw",
    "FL1100_IPP132_R4_Basal_even.raw",
    "FL1100_IPP132_R4_Basal_odd.raw",
    "FL1100_IPP132_R4_INFa_even.raw",
    "FL1100_IPP132_R4_INFa_odd.raw",
]

sample_ids = [
    "ECN",
    "ECN",
    "ECN_IFNa",
    "ECN_IFNa",
    "ECN",
    "ECN",
    "ECN_IFNa",
    "ECN_IFNa",
    "ECN",
    "ECN",
    "ECN_IFNa",
    "ECN_IFNa",
    "ECN",
    "ECN",
    "ECN_IFNa",
    "ECN_IFNa",
]

hla_types = [
    "A*03:01; A*02:01; B*40:01; B*49:01; C*03:04; C*07:01",
    "A*03:01; A*02:01; B*40:01; B*49:01; C*03:04; C*07:01",
    "A*03:01; A*02:01; B*40:01; B*49:01; C*03:04; C*07:01",
    "A*03:01; A*02:01; B*40:01; B*49:01; C*03:04; C*07:01",
    "A*03:01; A*02:01; B*40:01; B*49:01; C*03:04; C*07:01",
    "A*03:01; A*02:01; B*40:01; B*49:01; C*03:04; C*07:01",
    "A*03:01; A*02:01; B*40:01; B*49:01; C*03:04; C*07:01",
    "A*03:01; A*02:01; B*40:01; B*49:01; C*03:04; C*07:01",
    "A*03:01; A*02:01; B*40:01; B*49:01; C*03:04; C*07:01",
    "A*03:01; A*02:01; B*40:01; B*49:01; C*03:04; C*07:01",
    "A*03:01; A*02:01; B*40:01; B*49:01; C*03:04; C*07:01",
    "A*03:01; A*02:01; B*40:01; B*49:01; C*03:04; C*07:01",
    "A*03:01; A*02:01; B*40:01; B*49:01; C*03:04; C*07:01",
    "A*03:01; A*02:01; B*40:01; B*49:01; C*03:04; C*07:01",
    "A*03:01; A*02:01; B*40:01; B*49:01; C*03:04; C*07:01",
    "A*03:01; A*02:01; B*40:01; B*49:01; C*03:04; C*07:01",
]

beta_cell_df = pd.DataFrame({
    "tissue": cell_types,
    "file": raw_files,
    "uid": sample_ids,
    "hla": hla_types,
})
actual_final.append(beta_cell_df)

'''iPSC'''
cell_types = [
    "iPSC",
    "iPSC",
    "iPSC",
    "iPSC",
    "iPSC",
    "iPSC",
    "iPSC",
    "iPSC",
    "iPSC",
    "iPSC",
    "iPSC",
    "iPSC",
]

raw_files = [
    "CB_iPSC_Ctl_251119_1.raw",
    "CB_iPSC_Ctl_251119_2.raw",
    "CB_iPSC_IFN_251119_1.raw",
    "CB_iPSC_IFN_251119_2.raw",
    "Fibro_iPSC_Ctl_251119_1.raw",
    "Fibro_iPSC_Ctl_251119_2.raw",
    "Fibro_iPSC_IFN_251119_1.raw",
    "Fibro_iPSC_IFN_251119_2.raw",
    "iPSC22_IFN_280619_1.raw",
    "iPSC22_IFN_280619_2.raw",
    "iPSC_375M_290119_1.raw",
    "iPSC_375M_290119_2.raw",
]

sample_ids = [
    "Fibro-iPSC.2",
    "Fibro-iPSC.2",
    "Fibro-iPSC.2_IFN",
    "Fibro-iPSC.2_IFN",
    "Fibro-iPSC.1",
    "Fibro-iPSC.1",
    "Fibro-iPSC.1_IFN",
    "Fibro-iPSC.1_IFN",
    "hiPSC22_IFN",
    "hiPSC22_IFN",
    "hiPSC22",
    "hiPSC22",
]

hla_types = [
    "A*01:01; A*02:01; B*15:03; B*53:01; C*02:10; C*04:01",
    "A*01:01; A*02:01; B*15:03; B*53:01; C*02:10; C*04:01",
    "A*01:01; A*02:01; B*15:03; B*53:01; C*02:10; C*04:01",
    "A*01:01; A*02:01; B*15:03; B*53:01; C*02:10; C*04:01",
    "A*01:01; A*02:01; B*15:03; B*53:01; C*02:10; C*04:01",
    "A*01:01; A*02:01; B*15:03; B*53:01; C*02:10; C*04:01",
    "A*01:01; A*02:01; B*15:03; B*53:01; C*02:10; C*04:01",
    "A*01:01; A*02:01; B*15:03; B*53:01; C*02:10; C*04:01",
    "A*02:01; B*07:02; B*40:01; C*03:04; C*07:02",
    "A*02:01; B*07:02; B*40:01; C*03:04; C*07:02",
    "A*02:01; B*07:02; B*40:01; C*03:04; C*07:02",
    "A*02:01; B*07:02; B*40:01; C*03:04; C*07:02",
]

iPSC_df = pd.DataFrame({
    "tissue": cell_types,
    "file": raw_files,
    "uid": sample_ids,
    "hla": hla_types,
})
actual_final.append(iPSC_df)

'''Treg'''
cell_types = [
    "Treg",
    "Treg",
    "Treg",
    "Treg",
    "Treg",
    "Treg",
    "Treg",
    "Treg",
]

raw_files = [
    "FL0012779.raw",
    "FL0012783.raw",
    "FL0012803.raw",
    "FL0012807.raw",
    "FL0013327.raw",
    "FL0013331.raw",
    "FL0013353.raw",
    "FL0013357.raw",
]

sample_ids = [
    "donor02",
    "donor02",
    "donor03",
    "donor03",
    "donor04",
    "donor04",
    "donor05",
    "donor05",
]

hla_types = [
    "A*02:01",
    "A*02:01",
    "A*02:01",
    "A*02:01",
    "A*02:01",
    "A*02:01",
    "A*02:01",
    "A*02:01",
]

Treg_df = pd.DataFrame({
    "tissue": cell_types,
    "file": raw_files,
    "uid": sample_ids,
    "hla": hla_types,
})
actual_final.append(Treg_df)

'''immune cells'''
raw_files_by_tissue = {
    "CD14": [
        "20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_CD14_1_R1.raw",
        "20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_CD14_1_R2.raw",
        "20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_CD14_2_R1.raw",
        "20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_CD14_2_R2.raw",
        "20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_CD14_3_R1.raw",
        "20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_CD14_3_R2.raw",
        "20180424_QEh1_LC1_FaMa_HLAIp_D2_CD14_01_R1.raw",
        "20180424_QEh1_LC1_FaMa_HLAIp_D2_CD14_01_R2.raw",
        "20180428_QEh1_LC1_FaMa_HLAIp_D3_CD14_01_R1.raw",
        "20180428_QEh1_LC1_FaMa_HLAIp_D3_CD14_01_R2.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD14_1_R1.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD14_1_R2.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD14_2_R1.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD14_2_R2.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD14_3_R1.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD14_3_R2.raw",
    ],
    "ImmDC": [
        "20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_ImDC_1_R1.raw",
        "20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_ImDC_1_R2.raw",
        "20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_ImDC_2_R1.raw",
        "20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_ImDC_2_R2.raw",
        "20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_ImDC_3_R1.raw",
        "20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_ImDC_3_R2.raw",
        "20180424_QEh1_LC1_FaMa_HLAIp_D2_ImmDc_01_R1.raw",
        "20180424_QEh1_LC1_FaMa_HLAIp_D2_ImmDc_01_R2.raw",
        "20180428_QEh1_LC1_FaMa_HLAIp_D3_ImmDc_01_R1.raw",
        "20180428_QEh1_LC1_FaMa_HLAIp_D3_ImmDc_01_R2.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_ImmDC_1_R1.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_ImmDC_1_R2.raw",
    ],
    "MatureDC": [
        "20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_MDC_1_R1.raw",
        "20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_MDC_1_R2.raw",
        "20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_MDC_2_R1.raw",
        "20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_MDC_2_R2.raw",
        "20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_MDC_3_R1.raw",
        "20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_MDC_3_R2.raw",
        "20180424_QEh1_LC1_FaMa_HLAIp_D2_MaDc_01_R1.raw",
        "20180424_QEh1_LC1_FaMa_HLAIp_D2_MaDc_01_R2.raw",
        "20180428_QEh1_LC1_FaMa_HLAIp_D3_MaDc_01_R1.raw",
        "20180428_QEh1_LC1_FaMa_HLAIp_D3_MaDc_01_R2.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_MaDC_1_R1.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_MaDC_1_R2.raw",
    ],
    "CD4": [
        "20180424_QEh1_LC1_FaMa_HLAIp_D2_CD4_01_R1.raw",
        "20180424_QEh1_LC1_FaMa_HLAIp_D2_CD4_02_R1.raw",
        "20180424_QEh1_LC1_FaMa_HLAIp_D2_CD4_02_R2.raw",
        "20180424_QEh1_LC1_FaMa_HLAIp_D2_CD4_03_R1.raw",
        "20180424_QEh1_LC1_FaMa_HLAIp_D2_CD4_03_R2.raw",
        "20180428_QEh1_LC1_FaMa_HLAIp_D3_CD4_01_R1.raw",
        "20180428_QEh1_LC1_FaMa_HLAIp_D3_CD4_01_R2.raw",
        "20180428_QEh1_LC1_FaMa_HLAIp_D3_CD4_02_R1.raw",
        "20180428_QEh1_LC1_FaMa_HLAIp_D3_CD4_02_R2.raw",
        "20180428_QEh1_LC1_FaMa_HLAIp_D3_CD4_03_R1.raw",
        "20180428_QEh1_LC1_FaMa_HLAIp_D3_CD4_03_R2.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD4_1_R1.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD4_1_R2.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD4_2_R1.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD4_2_R2.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD4_3_R1.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD4_3_R2.raw",
    ],
    "CD8": [
        "20180424_QEh1_LC1_FaMa_HLAIp_D2_CD8_01_R1.raw",
        "20180424_QEh1_LC1_FaMa_HLAIp_D2_CD8_01_R2.raw",
        "20180424_QEh1_LC1_FaMa_HLAIp_D2_CD8_02_R1.raw",
        "20180424_QEh1_LC1_FaMa_HLAIp_D2_CD8_02_R2.raw",
        "20180424_QEh1_LC1_FaMa_HLAIp_D2_CD8_03_R1.raw",
        "20180428_QEh1_LC1_FaMa_HLAIp_D3_CD8_01_R1.raw",
        "20180428_QEh1_LC1_FaMa_HLAIp_D3_CD8_01_R2.raw",
        "20180428_QEh1_LC1_FaMa_HLAIp_D3_CD8_02_R1.raw",
        "20180428_QEh1_LC1_FaMa_HLAIp_D3_CD8_02_R2.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD8_1_R1.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD8_1_R2.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD8_2_R1.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD8_2_R2.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD8_3_R1.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD8_3_R2.raw",
    ],
    "CD19": [
        "20180424_QEh1_LC1_FaMa_HLAIp_D2_CD19_01_R1.raw",
        "20180424_QEh1_LC1_FaMa_HLAIp_D2_CD19_01_R2.raw",
        "20180428_QEh1_LC1_FaMa_HLAIp_D3_CD19_01_R1.raw",
        "20180428_QEh1_LC1_FaMa_HLAIp_D3_CD19_01_R2.raw",
        "20180428_QEh1_LC1_FaMa_HLAIp_D3_CD19_02_R1.raw",
        "20180428_QEh1_LC1_FaMa_HLAIp_D3_CD19_02_R2.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD19_1_R1.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD19_1_R2.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD19_2_R1.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD19_2_R2.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD19_3_R1.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD19_3_R2.raw",
    ],
    "CD4_Act": [
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD4_Act_1_R1.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD4_Act_1_R2.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD4_Act_2_R1.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD4_Act_2_R2.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD4_Act_3_R1.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD4_Act_3_R2.raw",
    ],
    "CD8_Act": [
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD8_Act_1_R1.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD8_Act_1_R2.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD8_Act_2_R1.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD8_Act_2_R2.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD8_Act_3_R1.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD8_Act_3_R2.raw",
    ],
    "CD19_Act": [
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD19_Act_1_R1.raw",
        "20180814_QEh1_LC1_SA_FaMa_HLAIp_CD19_Act_1_R2.raw",
    ],
}

date_to_donor = {
    "20171216": "D1",
    "20180424": "D2",
    "20180428": "D3",
    "20180814": "D4",
}

donor_hla = {
    "D1": "A*02:01; A*11:01; B*15:01; B*51:01; C*03:04; C*14:02",
    "D2": "A*02:01; B*40:01; B*44:03; C*02:02; C*03:04",
    "D3": "A*01:01; A*11:01; B*08:01; B*44:03; C*04:01; C*07:21",
    "D4": "A*01:01; A*32:01; B*08:01; B*40:01; C*03:04; C*07:01",
}

cell_types = []
raw_files = []
sample_ids = []
hla_types = []

for tissue, files in raw_files_by_tissue.items():
    for file_name in files:
        donor = date_to_donor[file_name[:8]]

        cell_types.append(tissue)
        raw_files.append(file_name)
        sample_ids.append(donor)
        hla_types.append(donor_hla[donor])

immune_cell_df = pd.DataFrame({
    "tissue": cell_types,
    "file": raw_files,
    "uid": sample_ids,
    "hla": hla_types,
})
actual_final.append(immune_cell_df)

actual_final = pd.concat(actual_final,axis=0)
actual_final.to_csv('final.txt',sep='\t',index=None)
actual_final.to_csv('/gpfs/data/yarmarkovichlab/public/ImmunoVerse/database/final.txt',sep='\t',index=None)


