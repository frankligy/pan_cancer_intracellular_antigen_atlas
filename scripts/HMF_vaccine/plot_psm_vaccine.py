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


def cosine_similarity(a,b):
    dot = np.dot(a, b)
    norma = np.linalg.norm(a)
    normb = np.linalg.norm(b)
    cos = dot / (norma * normb)
    return cos




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

def draw_spectrum(mzml_path,scan,pep,ma):

    for i,spectrum in enumerate(mzml.read(mzml_path)):
        pat = r'scan=(\d+)'
        spectrum_id = re.search(pat,spectrum['id']).group(1)
        if int(spectrum_id) == scan:
            mz = spectrum['m/z array']
            intensity = spectrum['intensity array']
            rt = spectrum['scanList']['scan'][0]['scan start time']
            precursor_mz = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
            precursor_charge = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state']
            sus_spectrum = sus.MsmsSpectrum(spectrum['id'], precursor_mz, precursor_charge, mz, intensity, rt)
            break
    spectrum = sus_spectrum

    if ma == 'ITMS':
        tol_mass, tol_mode = [0.5,'Da']
    elif ma == 'FTMS':
        tol_mass, tol_mode = [20,'ppm']
    elif ma == 'TOF':
        tol_mass, tol_mode = [25,'ppm']
    elif ma == 'Unknown':
        tol_mass, tol_mode = [20,'ppm']


    p_fragment_tol_mass, p_fragment_tol_mode = tol_mass, tol_mode
    f_fragment_tol_mass, f_fragment_tol_mode = tol_mass, tol_mode
    peptide = pep


    spectrum = (
        spectrum.set_mz_range(min_mz=0, max_mz=5000)
        .remove_precursor_peak(p_fragment_tol_mass, p_fragment_tol_mode)
        .filter_intensity(min_intensity=0.0001, max_num_peaks=5000)
        .annotate_proforma(
            peptide, f_fragment_tol_mass, f_fragment_tol_mode, ion_types="aby", neutral_losses={"NH3": -17.026549, "H2O": -18.010565}
        )
    )

    fig, ax = plt.subplots(figsize=(12, 6))
    sup.spectrum(spectrum, annot_fmt=annotate_ion_type, grid=False, ax=ax)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_title('{};{};{}'.format(pep,pep2tpm[pep],ma))
    plt.savefig(os.path.join('./spectra','spectrum_{}_{}.png'.format(pep,scan)), bbox_inches='tight')
    plt.close()

def draw_spectrum_mirror(mzml_path,scan,pep,ma):

    # exp spectrum
    for i,spectrum in enumerate(mzml.read(mzml_path)):
        pat = r'scan=(\d+)'
        spectrum_id = re.search(pat,spectrum['id']).group(1)
        if int(spectrum_id) == scan:
            mz = spectrum['m/z array']
            intensity = spectrum['intensity array']
            intensity = [np.interp(item,xp=[intensity.min(),intensity.max()],fp=[0,1]) for item in intensity]
            rt = spectrum['scanList']['scan'][0]['scan start time']
            precursor_mz = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
            precursor_charge = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state']
            sus_spectrum = sus.MsmsSpectrum(spectrum['id'], precursor_mz, precursor_charge, mz, intensity, rt)
            break
    spectrum = sus_spectrum

    if ma == 'ITMS':
        tol_mass, tol_mode = [0.5,'Da']
    elif ma == 'FTMS':
        tol_mass, tol_mode = [20,'ppm']
    elif ma == 'TOF':
        tol_mass, tol_mode = [25,'ppm']
    elif ma == 'Unknown':
        tol_mass, tol_mode = [20,'ppm']


    p_fragment_tol_mass, p_fragment_tol_mode = tol_mass, tol_mode
    f_fragment_tol_mass, f_fragment_tol_mode = tol_mass, tol_mode
    peptide = pep

    spectrum = (
        spectrum.set_mz_range(min_mz=0, max_mz=5000)
        .remove_precursor_peak(p_fragment_tol_mass, p_fragment_tol_mode)
        .filter_intensity(min_intensity=0.0000, max_num_peaks=5000)
        .annotate_proforma(
            peptide, f_fragment_tol_mass, f_fragment_tol_mode, ion_types="by"
        )
    )

    exp_anno = {}
    for i1,i2 in zip(list(spectrum.annotation),list(spectrum.intensity)):
        i1 = str(i1)
        if i1 != '?':
            i1 = i1.split('/')[0]
            if i1 not in exp_anno.keys():
                exp_anno[i1] = i2
            else:
                if i2 > exp_anno[i1]:  
                    exp_anno[i1] = i2

    # synthetic spectrum
    df = pd.read_csv(pep2path[pep],sep='\t')
    mz = df['mz'].values
    intensity = df['predicted'].values
    intensity = [np.interp(item,xp=[intensity.min(),intensity.max()],fp=[0,1]) for item in intensity]
    sus_spectrum = sus.MsmsSpectrum('placehoder', precursor_mz, precursor_charge, mz, intensity, rt)
    syn_spectrum = sus_spectrum

    syn_spectrum = (
        syn_spectrum.set_mz_range(min_mz=0, max_mz=5000)
        .remove_precursor_peak(p_fragment_tol_mass, p_fragment_tol_mode)
        .filter_intensity(min_intensity=0.0000, max_num_peaks=5000)
        .annotate_proforma(
            peptide, f_fragment_tol_mass, f_fragment_tol_mode, ion_types="aby", neutral_losses={"NH3": -17.026549, "H2O": -18.010565}
        )
    )

    syn_anno = {}
    for i1,i2 in zip(list(syn_spectrum.annotation),list(syn_spectrum.intensity)):
        i1 = str(i1)
        if i1 != '?':
            i1 = i1.split('/')[0]
            if i1 not in syn_anno.keys():
                syn_anno[i1] = i2
            else:
                if i2 > syn_anno[i1]:  
                    syn_anno[i1] = i2

    # calcualte metrics
    print(exp_anno)
    print(syn_anno)
    store_intensity_exp = []
    store_intensity_syn = []
    common_by = set(exp_anno.keys()).intersection(set(syn_anno.keys()))
    for common in common_by:
        store_intensity_exp.append(exp_anno[common])
        store_intensity_syn.append(syn_anno[common])
    p, s = pearsonr(x=store_intensity_exp,y=store_intensity_syn)
    cos = cosine_similarity(store_intensity_exp,store_intensity_syn)

    fig, ax = plt.subplots(figsize=(12, 6))
    sup.mirror(spectrum, syn_spectrum, ax=ax)
    fig.suptitle('{}\nPearson r:{}\np-value:{}\ncosine:{}'.format(pep,round(p,2),round(s,2),cos),fontsize=6)
    ax.grid(False)
    plt.savefig('mirror_{}_{}.pdf'.format(pep,scan), bbox_inches="tight", transparent=True)
    plt.close()






final = pd.read_csv('/gpfs/data/yarmarkovichlab/HMF_vaccine/immunoverse_result/antigen/other_alg/final_enhanced.txt',sep='\t')
cond = [False if ('[]' in item) and ('(\'HLA-' not in item) else True for item in final['presented_by_each_sample_hla']]
final = final.loc[cond,:]
final = final.loc[final['typ']=='variant',:]
final = final.loc[final['source'].str.contains('wgs_'),:]
gene = pd.read_csv('/gpfs/data/yarmarkovichlab/HMF_vaccine/immunoverse_result/genes_lfc.txt',sep='\t',index_col=0)
ensg2tpm = gene['median_tumor'].to_dict()
final['tpm'] = [ensg2tpm[item.split('|')[0]] for item in final['source']]
final.to_csv('check.txt',sep='\t',index=None)
valid_peps = final['pep'].values.tolist()
pep2tpm = pd.Series(index=final['pep'].values,data=final['tpm'].values).to_dict()


mzml_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/immunopeptidome/melanoma/mzml_dir'
cmd = 'find {} -maxdepth 2 -type f -name "*.mzML" -exec echo {{}} \;'.format(mzml_dir)
all_mzml = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
mapping = {}
for mzml_ in all_mzml:
    f = mzml_.split('/')[-1]
    d = mzml_.split('/')[-2]
    mapping[f.split('.mzML')[0]] = d

master = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas/SKCM/antigen/fdr/msmsScans_all_add_tesorai.txt'


variant_df = pd.read_csv('/gpfs/data/yarmarkovichlab/HMF_vaccine/immunoverse_result/antigen/other_alg/variant_neoantigen.txt',sep='\t')
variant_df = variant_df.loc[variant_df['Sequence'].isin(valid_peps),:]

query_peps = ['KSAIETVL','KIRKSNVL','FLRTFLGGRF','APTRTAAL']
variant_df = variant_df.loc[variant_df['Sequence'].isin(query_peps),:]

pep2path = {
    'KSAIETVL':'ms2pip_prediction_KSAIETVL_2.tsv',
    'KIRKSNVL':'ms2pip_prediction_KIRKSNVL_3.tsv',
    'FLRTFLGGRF':'ms2pip_prediction_FLRTFLGGRF_3.tsv',
    'APTRTAAL':'ms2pip_prediction_APTRTAAL_2.tsv',
}

for scan,raw,pep in zip(variant_df['Scan number'],variant_df['Raw file'],variant_df['Sequence']):
    mzml_path = os.path.join(mzml_dir,mapping[raw],'{}.mzML'.format(raw))
    cmd = 'grep {} {} | head -n 5 > tmp.txt'.format(raw,master)
    subprocess.run(cmd,shell=True)
    try:
        tmp_df = pd.read_csv('tmp.txt',sep='\t')
    except:
        ma = 'ITMS'
    else:
        ma = tmp_df.iloc[:,23].iloc[0]
    draw_spectrum_mirror(mzml_path,scan,pep,ma)






