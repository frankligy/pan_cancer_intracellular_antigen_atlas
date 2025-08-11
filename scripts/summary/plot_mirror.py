#!/gpfs/data/yarmarkovichlab/Frank/pan_cancer/antigen_portal/spectrum_env/bin/python3.8

import sys
import os
import numpy as np
import pandas as pd
from pyteomics import mzml
import matplotlib.pyplot as plt
import matplotlib as mpl
import multiprocessing as mp
import subprocess
from tqdm import tqdm
import json
import argparse 
from dash import Dash, Input, Output, callback, dash_table, html, dcc, State
import re
from ast import literal_eval
from scipy.stats import pearsonr

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus


def annotate_ion_type(annotation, ion_types="abym"):
    if annotation.ion_type[0] in ion_types:
        if abs(annotation.isotope) == 1:
            iso = "+i" if annotation.isotope > 0 else "-i"
        elif annotation.isotope != 0:
            iso = f"{annotation.isotope:+}i"
        else:
            iso = ""
        nl = {"-NH3": "*", "-H2O": "o"}.get(annotation.neutral_loss, "")
        return f"{annotation.ion_type}{iso}{'+' * annotation.charge}{nl}"
    else:
        return ""


def draw_spectrum_mirror(mzml_path1,scan1,pep1,ma1,mzml_path2,scan2,pep2,ma2,taskname):

    # exp spectrum
    for i,spectrum in enumerate(mzml.read(mzml_path1)):
        pat = r'scan=(\d+)'
        spectrum_id = re.search(pat,spectrum['id']).group(1)
        if int(spectrum_id) == scan1:
            mz = spectrum['m/z array']
            intensity = spectrum['intensity array']
            intensity = [np.interp(item,xp=[intensity.min(),intensity.max()],fp=[0,1]) for item in intensity]
            rt = spectrum['scanList']['scan'][0]['scan start time']
            precursor_mz = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
            precursor_charge = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state']
            sus_spectrum = sus.MsmsSpectrum(spectrum['id'], precursor_mz, precursor_charge, mz, intensity, rt)
            break
    exp_spectrum = sus_spectrum

    if ma1 == 'ITMS':
        tol_mass, tol_mode = [0.5,'Da']
    elif ma1 == 'FTMS':
        tol_mass, tol_mode = [20,'ppm']
    elif ma1 == 'TOF':
        tol_mass, tol_mode = [25,'ppm']
    elif ma1 == 'Unknown':
        tol_mass, tol_mode = [20,'ppm']


    p_fragment_tol_mass, p_fragment_tol_mode = tol_mass, tol_mode
    f_fragment_tol_mass, f_fragment_tol_mode = tol_mass, tol_mode
    peptide = pep1

    exp_spectrum = (
        exp_spectrum.set_mz_range(min_mz=0, max_mz=5000)
        .filter_intensity(min_intensity=0.0000, max_num_peaks=5000)
        .annotate_proforma(
            peptide, f_fragment_tol_mass, f_fragment_tol_mode, ion_types="abymp"
        )
    )

    # syn spectrum
    for i,spectrum in enumerate(mzml.read(mzml_path2)):
        pat = r'scan=(\d+)'
        spectrum_id = re.search(pat,spectrum['id']).group(1)
        if int(spectrum_id) == scan2:
            mz = spectrum['m/z array']
            intensity = spectrum['intensity array']
            intensity = [np.interp(item,xp=[intensity.min(),intensity.max()],fp=[0,1]) for item in intensity]
            rt = spectrum['scanList']['scan'][0]['scan start time']
            precursor_mz = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
            precursor_charge = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state']
            sus_spectrum = sus.MsmsSpectrum(spectrum['id'], precursor_mz, precursor_charge, mz, intensity, rt)
            break
    syn_spectrum = sus_spectrum

    if ma2 == 'ITMS':
        tol_mass, tol_mode = [0.5,'Da']
    elif ma2 == 'FTMS':
        tol_mass, tol_mode = [20,'ppm']
    elif ma2 == 'TOF':
        tol_mass, tol_mode = [25,'ppm']
    elif ma2 == 'Unknown':
        tol_mass, tol_mode = [20,'ppm']


    p_fragment_tol_mass, p_fragment_tol_mode = tol_mass, tol_mode
    f_fragment_tol_mass, f_fragment_tol_mode = tol_mass, tol_mode
    peptide = pep2

    syn_spectrum = (
        syn_spectrum.set_mz_range(min_mz=0, max_mz=5000)
        .filter_intensity(min_intensity=0.0000, max_num_peaks=5000)
        .annotate_proforma(
            peptide, f_fragment_tol_mass, f_fragment_tol_mode, ion_types="abymp"
        )
    )

    fig, ax = plt.subplots(figsize=(12, 6))
    sup.mirror(exp_spectrum, syn_spectrum, ax=ax)
    ax.grid(False)
    plt.savefig('mirror_{}.pdf'.format(taskname), bbox_inches="tight", transparent=True)
    plt.close()


# batch1
# znf749_5utr
mzml_path1 = 'PRL-OGB-1550-b1/SK-N-AS_PDX_W632_1_10086_16612.mzML'
scan1 = 15996
pep1 = 'RYLPSSVFL'
ma1 = 'FTMS'
mzml_path2 = 'PRL-OGB-1550-b1/20240627_E_OdinLC_IC_synpep_2pmol_oldstock.mzML'
scan2 = 15185
pep2 = 'RYLPSSVFL'
ma2 = 'FTMS'
draw_spectrum_mirror(mzml_path1,scan1,pep1,ma1,mzml_path2,scan2,pep2,ma2,'znf749_5utr')

# c19orf48_pseudo
mzml_path1 = 'PRL-OGB-1550-b1/SK-N-AS_PDX_W632_1_10086_16612.mzML'
scan1 = 14179
pep1 = 'AYPASLQTL'
ma1 = 'FTMS'
mzml_path2 = 'PRL-OGB-1550-b1/20240627_E_OdinLC_IC_synpep_2pmol_oldstock.mzML'
scan2 = 13233
pep2 = 'AYPASLQTL'
ma2 = 'FTMS'
draw_spectrum_mirror(mzml_path1,scan1,pep1,ma1,mzml_path2,scan2,pep2,ma2,'c19orf48_pseudo')

# sPMEL
mzml_path1 = 'PRL-OGB-1550-b1/YG20190108_AA_MEL_sel_Pat2_01_27626_28593.mzML'
scan1 = 28052
pep1 = 'KTWDQVPFSV'
ma1 = 'FTMS'
mzml_path2 = 'PRL-OGB-1550-b1/20240627_E_OdinLC_IC_synpep_2pmol_oldstock.mzML'
scan2 = 15183
pep2 = 'KTWDQVPFSV'
ma2 = 'FTMS'
draw_spectrum_mirror(mzml_path1,scan1,pep1,ma1,mzml_path2,scan2,pep2,ma2,'sPMEL')

# polr2k_5utr
mzml_path1 = 'PRL-OGB-1550-b1/SK-N-AS_PDX_W632_1_10086_16612.mzML'
scan1 = 12056
pep1 = 'LYLETRSEF'
ma1 = 'FTMS'
mzml_path2 = 'PRL-OGB-1550-b1/20240627_E_OdinLC_IC_synpep_2pmol_oldstock.mzML'
scan2 = 11466
pep2 = 'LYLETRSEF'
ma2 = 'FTMS'
draw_spectrum_mirror(mzml_path1,scan1,pep1,ma1,mzml_path2,scan2,pep2,ma2,'polr2k_5utr')

# tmem203_3utr
mzml_path1 = 'PRL-OGB-1550-b1/SK-N-AS_PDX_W632_1_10086_16612.mzML'
scan1 = 10478
pep1 = 'STIRVLSGY'
ma1 = 'FTMS'
mzml_path2 = 'PRL-OGB-1550-b1/20240627_E_OdinLC_IC_synpep_2pmol_oldstock.mzML'
scan2 = 9890
pep2 = 'STIRVLSGY'
ma2 = 'FTMS'
draw_spectrum_mirror(mzml_path1,scan1,pep1,ma1,mzml_path2,scan2,pep2,ma2,'tmem203_3utr')


# slc45a2_splice
mzml_path1 = 'PRL-OGB-1550-b1/20180222_QE_HFX_LC2_HLAIp_JMI_SA_Me290_R2_15324_15742.mzML'
scan1 = 15471
pep1 = 'FTDSQGNDIK'
ma1 = 'FTMS'
mzml_path2 = 'PRL-OGB-1550-b1/20240627_E_OdinLC_IC_synpep_2pmol_oldstock.mzML'
scan2 = 4152
pep2 = 'FTDSQGNDIK'
ma2 = 'FTMS'
draw_spectrum_mirror(mzml_path1,scan1,pep1,ma1,mzml_path2,scan2,pep2,ma2,'slc45a2_splice')


# batch2
# fuso_tid
mzml_path1 = 'PRL-OGB-1667-b2/HCT116MG_D190212_S2_14578.mzML'
scan1 = 14578
pep1 = 'TIDELQKI'
ma1 = 'FTMS'
mzml_path2 = 'PRL-OGB-1667-b2/20250213_E_Evo_IC_synpep_14_2pmol_rerun2.mzML'
scan2 = 12124
pep2 = 'TIDELQKI'
ma2 = 'FTMS'
draw_spectrum_mirror(mzml_path1,scan1,pep1,ma1,mzml_path2,scan2,pep2,ma2,'fuso_tid')

# fuso_lsd
mzml_path1 = 'PRL-OGB-1667-b2/20171109_QEh1_LC1_HLApI_SA_JMI_PDX146_2_R1_28094.mzML'
scan1 = 28094
pep1 = 'LSDLGSGIYR'
ma1 = 'FTMS'
mzml_path2 = 'PRL-OGB-1667-b2/20250213_E_Evo_IC_synpep_14_2pmol_rerun2.mzML'
scan2 = 11700
pep2 = 'LSDLGSGIYR'
ma2 = 'FTMS'
draw_spectrum_mirror(mzml_path1,scan1,pep1,ma1,mzml_path2,scan2,pep2,ma2,'fuso_lsd')

# fuso_iek
mzml_path1 = 'PRL-OGB-1667-b2/20120518_EXQ0_MiBa_SA_HCT116_mHLA-2_14582.mzML'
scan1 = 14582
pep1 = 'IEKEVISKY'
ma1 = 'FTMS'
mzml_path2 = 'PRL-OGB-1667-b2/20250213_E_Evo_IC_synpep_14_2pmol_rerun2.mzML'
scan2 = 7579
pep2 = 'IEKEVISKY'
ma2 = 'FTMS'
draw_spectrum_mirror(mzml_path1,scan1,pep1,ma1,mzml_path2,scan2,pep2,ma2,'fuso_iek')

# nia_QPK
mzml_path1 = 'PRL-OGB-1667-b2/OvCa104_classI_Rep4_8594.mzML'
scan1 = 8594
pep1 = 'QPKTKLLLL'
ma1 = 'ITMS'
mzml_path2 = 'PRL-OGB-1667-b2/20250606_E_EvoJ_c04_IC_synpep_rerun_multiactivation_20250607133024.mzML'
scan2 = 12416
pep2 = 'QPKTKLLLL'
ma2 = 'ITMS'
draw_spectrum_mirror(mzml_path1,scan1,pep1,ma1,mzml_path2,scan2,pep2,ma2,'nia_qpk')

# nia_ldi
mzml_path1 = 'PRL-OGB-1667-b2/OvCa45_classI_Rep2_16702.mzML'
scan1 = 16702
pep1 = 'LDIHTFGLYY'
ma1 = 'ITMS'
mzml_path2 = 'PRL-OGB-1667-b2/20250213_E_Evo_IC_synpep_14_2pmol_rerun2.mzML'
scan2 = 21112
pep2 = 'LDIHTFGLYY'
ma2 = 'ITMS'
draw_spectrum_mirror(mzml_path1,scan1,pep1,ma1,mzml_path2,scan2,pep2,ma2,'nia_ldi')

# nia_klk
mzml_path1 = 'PRL-OGB-1667-b2/OvCa53_classI_Rep1_1456.mzML'
scan1 = 1456
pep1 = 'KLKPGILKK'
ma1 = 'ITMS'
mzml_path2 = 'PRL-OGB-1667-b2/20250606_E_EvoJ_c04_IC_synpep_rerun_multiactivation_20250607133024.mzML'
scan2 = 2652
pep2 = 'KLKPGILKK'
ma2 = 'ITMS'
draw_spectrum_mirror(mzml_path1,scan1,pep1,ma1,mzml_path2,scan2,pep2,ma2,'nia_klk')

# nia_lvg
mzml_path1 = 'PRL-OGB-1667-b2/20170608_QEh1_LC1_ChCh_FAMA_SA_HLAIp_UWB1289_ctrl_3_R2_18302.mzML'
scan1 = 18302
pep1 = 'LVGPNGVGK'
ma1 = 'FTMS'
mzml_path2 = 'PRL-OGB-1667-b2/20250213_E_Evo_IC_synpep_14_2pmol_rerun2.mzML'
scan2 = 6470
pep2 = 'LVGPNGVGK'
ma2 = 'FTMS'
draw_spectrum_mirror(mzml_path1,scan1,pep1,ma1,mzml_path2,scan2,pep2,ma2,'nia_lvg')

# orf2_ria
mzml_path1 = 'PRL-OGB-1667-b2/20171107_QEh1_LC1_HLApI_SA_JMI_CRC-04_3_R1_15694.mzML'
scan1 = 15694
pep1 = 'RIAKSILSQK'
ma1 = 'FTMS'
mzml_path2 = 'PRL-OGB-1667-b2/20250213_E_Evo_IC_synpep_14_2pmol_rerun2.mzML'
scan2 = 4768
pep2 = 'RIAKSILSQK'
ma2 = 'FTMS'
draw_spectrum_mirror(mzml_path1,scan1,pep1,ma1,mzml_path2,scan2,pep2,ma2,'orf2_ria')

# orf2_ilp
mzml_path1 = 'PRL-OGB-1667-b2/20180831_QEh1_LC1_SA_JMI_HLAIp_CRC-04_IFN2_R01_37367.mzML'
scan1 = 37367
pep1 = 'ILPKVIYRF'
ma1 = 'FTMS'
mzml_path2 = 'PRL-OGB-1667-b2/20250213_E_Evo_IC_synpep_14_2pmol_rerun2.mzML'
scan2 = 13483
pep2 = 'ILPKVIYRF'
ma2 = 'FTMS'
draw_spectrum_mirror(mzml_path1,scan1,pep1,ma1,mzml_path2,scan2,pep2,ma2,'orf2_ilp')

# orf2_sgy
mzml_path1 = 'PRL-OGB-1667-b2/20170708_QEh1_LC1_ChCh_SA_HLAIp_CRC-08_1_R1_19326.mzML'
scan1 = 19326
pep1 = 'SGYKINVQK'
ma1 = 'FTMS'
mzml_path2 = 'PRL-OGB-1667-b2/20250213_E_Evo_IC_synpep_14_2pmol_rerun2.mzML'
scan2 = 6029
pep2 = 'SGYKINVQK'
ma2 = 'FTMS'
draw_spectrum_mirror(mzml_path1,scan1,pep1,ma1,mzml_path2,scan2,pep2,ma2,'orf2_sgy')

'''batch3'''
# HAVCR1 +2
mzml_path1 = 'PRL-OGB-1740-b3/YJ20180316_SK_HLA_RCC_Pat9_IFNy_4Ips_noAPD_R1_02.mzML'
scan1 = 14091
pep1 = 'DLSRRDVSL'
ma1 = 'FTMS'
mzml_path2 = 'PRL-OGB-1740-b3/20250606_E_EvoJ_c04_IC_synpep_20250607124157.mzML'
scan2 = 10650
pep2 = 'DLSRRDVSL'
ma2 = 'FTMS'
draw_spectrum_mirror(mzml_path1,scan1,pep1,ma1,mzml_path2,scan2,pep2,ma2,'havcr1_2')

# HAVCR1 +3
mzml_path1 = 'PRL-OGB-1740-b3/YJ20180316_SK_HLA_RCC_Pat9_IFNy_4Ips_noAPD_R1_02.mzML'
scan1 = 14108
pep1 = 'DLSRRDVSL'
ma1 = 'FTMS'
mzml_path2 = 'PRL-OGB-1740-b3/20250606_E_EvoJ_c04_IC_synpep_20250607124157.mzML'
scan2 = 10704
pep2 = 'DLSRRDVSL'
ma2 = 'FTMS'
draw_spectrum_mirror(mzml_path1,scan1,pep1,ma1,mzml_path2,scan2,pep2,ma2,'havcr1_3')

# LARP4B_fs
mzml_path1 = 'PRL-OGB-1740-b3/20180911_QEh1_LC1_SA_JMI_HLAIp_CRC-05_TRAM2_R01.mzML'
scan1 = 22057
pep1 = 'SAYLGRTL'
ma1 = 'FTMS'
mzml_path2 = 'PRL-OGB-1740-b3/20250606_E_EvoJ_c04_IC_synpep_20250607124157.mzML'
scan2 = 12526
pep2 = 'SAYLGRTL'
ma2 = 'FTMS'
draw_spectrum_mirror(mzml_path1,scan1,pep1,ma1,mzml_path2,scan2,pep2,ma2,'LARP4B_fs')

# CTNNB1_mis
mzml_path1 = 'PRL-OGB-1740-b3/2060_L_MdB069_1mei19_donor_MdB029_HCC_OT.mzML'
scan1 = 16818
pep1 = 'YLDSGIHSGATA'
ma1 = 'FTMS'
mzml_path2 = 'PRL-OGB-1740-b3/20250606_E_EvoJ_c04_IC_synpep_20250607124157.mzML'
scan2 = 8845
pep2 = 'YLDSGIHSGATA'
ma2 = 'FTMS'
draw_spectrum_mirror(mzml_path1,scan1,pep1,ma1,mzml_path2,scan2,pep2,ma2,'CTNNB1_mis')
