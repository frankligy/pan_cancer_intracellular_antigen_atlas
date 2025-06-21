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
import dash

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'




def hla_formatting(pre,pre_type,post_type):
    if pre_type == 'netMHCpan_output' and post_type == 'netMHCpan_input':  # HLA-A*01:01 to HLA-A01:01
        post = [hla.replace('*','') for hla in pre]
    elif pre_type == 'netMHCpan_input' and post_type == 'netMHCpan_output':  # HLA-A01:01 to HLA-A*01:01
        post = [hla[:5] + '*' + hla[5:] for hla in pre]
    elif pre_type == 'netMHCpan_output' and post_type == 'deepimmuno':  # HLA-A*01:01 to HLA-A*0101
        post = [hla.replace(':','') for hla in pre]
    elif pre_type == 'deepimmuno' and post_type == 'netMHCpan_output': # HLA-A*0101 to HLA-A*01:01
        post = []
        for hla in pre:
            first = hla[:6] # HLA-A*
            second = hla[6:8] # 01, it is so far always 2 digit, no case exceed 2 digit
            third = hla[8:] # 01 or 101, it could be 3 digit
            now = first + second + ':' + third
            post.append(now)
    elif pre_type == 'deepimmuno_nostar' and post_type == 'netMHCpan_input':   # HLA-A0101 to HLA-A01:01
        post = []
        for hla in pre:
            first = hla[:5] # HLA-A
            second = hla[5:7] # 01, it is so far always 2 digit, no case exceed 2 digit
            third = hla[7:] # 01 or 101, it could be 3 digit
            now = first + second + ':' + third
            post.append(now)
    return post


# @callback(Output('url','search'),State('cancer_dropdown','value'),State('query_type','value'),State('query_peptide','value'),State('query_source','value'),State('sort_dropdown','value'),Input('submit_button','n_clicks'))
# def update_url(cancer,query_type,query_peptide,query_source,sort_by_column,n_clicks):
#     return '?cancer={}&type={}&peptide={}&source={}&sort={}'.format(cancer,query_type,query_peptide,query_source,sort_by_column)

# @callback(Output('cancer_dropdown', 'value'),Output('query_type', 'value'),Output('query_peptide', 'value'),Output('query_source', 'value'),Output('sort_dropdown', 'value'),Input('url', 'search'))
# def load_from_url(search):

#     parsed = dict(item.split('=') for item in search.lstrip('?').split('&'))

#     return parsed['cancer'],parsed['type'],parsed['peptide'],parsed['source'],parsed['sort']


@callback(Output('candidate','data'),Output('downloader','data'),State('cancer_dropdown','value'),State('query_type','value'),State('query_peptide','value'),State('query_source','value'),State('query_hla','value'),State('sort_dropdown','value'),State('download_table','value'),Input('submit_button','n_clicks'))
def filter_table(cancer,query_type,query_peptide,query_source,query_hla,sort_by_column,to_download,n_clicks):

    if cancer == 'All':
        lis = []
        all_cancers = ['BRCA','KIRC','COAD','STAD','MESO','LIHC','ESCA','CESC','BLCA','RT','AML','DLBC','GBM','NBL','PAAD','HNSC','OV','LUSC','LUAD','CHOL','SKCM']
        for c in all_cancers:
            tmp = pd.read_csv('./static/{}_final_enhanced.txt'.format(c),sep='\t')
            lis.append(tmp)
        final = pd.concat(lis,axis=0,keys=all_cancers).reset_index(level=-2).rename(columns={'level_0':'cancer'})

    else:
        final = pd.read_csv('./static/{}_final_enhanced.txt'.format(cancer),sep='\t')
        final.insert(0,'cancer',np.full(shape=final.shape[0],fill_value=cancer))

    cond = [False if ('[]' in item) and ('(\'HLA-' not in item) else True for item in final['presented_by_each_sample_hla']]
    final = final.loc[cond,:]
    selected_columns = ['cancer','pep','typ','source','highest_score','depmap_median','sc_pert','presented_by_each_sample_hla','additional_query','detailed_intensity']
    final = final.loc[:,selected_columns]
    final.columns = ['Cancer','Peptide','Antigen Class','Source','Spectral Score','Essentiality','Homogeneity','presented_by_each_sample_hla','additional_query','detailed_intensity']


    # round
    for column in ['Essentiality','Homogeneity']:
        col = []
        for item in final[column]:
            if isinstance(item,str):   # not empty
                item = literal_eval(item)
                if isinstance(item,float):  # single
                    item = round(item,2)
                elif isinstance(item,list):  # multiple
                    item = str([i if i is None else round(i,2) for i in item])
            col.append(item)
        final[column] = col

    # simplify self gene source
    col = []
    for t,s in zip(final['Antigen Class'],final['Source']):
        if t == 'self_gene':
            col.append(simplify_self_gene_source(s))
        else:
            col.append(s)
    final.insert(3,'Simplified Source',col)

    # add intensity
    col = [round(np.median(literal_eval(item)),2) for item in final['detailed_intensity']]
    final.insert(5,'Abundance',col)
    final.drop(columns='detailed_intensity',inplace=True)

    # remove underscore for type
    final['Antigen Class'] = [item.replace('_',' ') for item in final['Antigen Class']]


    # now act
    if n_clicks == 0:
        return final.to_dict('records'), dash.no_update

    # filter
    if query_type != 'All':
        final = final.loc[final['Antigen Class']==query_type,:]

    if isinstance(query_source,str):
        query_source = query_source.upper()
        final = final.loc[final['Source'].str.contains(query_source),:]

    if isinstance(query_peptide,str):
        query_peptide = query_peptide.upper()
        final = final.loc[final['Peptide']==query_peptide,:]

    if isinstance(query_hla,str):
        query_hla = query_hla.upper()
        final = final.loc[(final['presented_by_each_sample_hla'].str.contains(query_hla)) | (final['additional_query'].str.contains(query_hla)),:]


    # sort
    is_ascending = True
    if sort_by_column in ['Spectral Score','Abundance']:
        is_ascending = False
    final = final.sort_values(by=sort_by_column,ascending=is_ascending)

    # download
    if to_download == 'No':
        second_return = dash.no_update
    elif to_download == 'Yes':
        second_return = dcc.send_data_frame(final.to_csv,"ImmunoVerse_table.csv")

    return final.to_dict('records'),second_return


def get_hla_info(meta,lists,hlas):
    for k,vs in hlas.items():
        if len(vs) == 0:
            continue
        else:
            for v in vs:
                if v[0] is not None:
                    lists.append(v)
    lists = list(set(lists))
    df = pd.DataFrame.from_records(lists,columns=['hla','rank_pert','nM','id'])
    hla = [item.replace('*','') for item in df['hla'].tolist()]
    df['freq'] = [hla2freq.get(item,float('nan')) for item in hla]

    # add hla-agnostic
    tmp = 1
    for f in df['freq']:
        if np.isnan(f):
            f = 0
        tmp *= (1-f)
    all_freq = 1 - tmp
    all_freq = round(all_freq,2)
    added_df = pd.DataFrame({'hla':['all_hla'],'rank_pert':[None],'nM':[None],'id':[None],'freq':[round(all_freq,5)]})
    df = pd.concat([df,added_df],axis=0)

    # adding recurrency
    hla2info = {}
    considered = {}  # {hla:[s1,s2]}
    for k,vs in hlas.items():
        for v in vs:
            if len(v) == 4:
                hla,pert,nm,id_ = v
                if hla is not None:
                    considered.setdefault(hla,[]).append(k)
    for v in lists:
        hla,pert,nm,id_ = v
        if hla not in considered.keys():
            considered[hla] = []
    considered = {k:list(set(v)) for k,v in considered.items()}
    for k,vs in considered.items():
        hla = k.split('HLA-')[1]
        sub_meta = meta.loc[meta['HLA'].notna(),:]
        samples = sub_meta.loc[sub_meta['HLA'].str.contains(hla),:]['biology'].unique().tolist()
        try:
            frac = len(vs)/len(samples)
        except ZeroDivisionError:
            frac = 0
        hla2info[k] = (len(vs),len(samples),round(frac,2))

    # add hla agnostic
    all_hla = [item.split('HLA-')[1] for item in considered.keys()]
    all_sample = []
    for item in considered.values():
        all_sample.extend(item)
    all_sample = list(set(all_sample))
    sub_meta = meta.loc[meta['HLA'].notna(),:]
    cond = []
    for item in sub_meta['HLA']:
        flag = False
        for hla in all_hla:
            if hla in item:
                flag = True
                break
        cond.append(flag)
    samples = sub_meta.loc[cond,:]['biology'].unique().tolist()
    try:
        frac = len(all_sample)/len(samples)
    except ZeroDivisionError:
        frac = 0
    hla2info['all_hla'] = (len(all_sample),len(samples),round(frac,2))

    col1,col2,col3 = [],[],[]
    for hla in df['hla']:
        info = hla2info.get(hla,(0,0,0))
        col1.append(info[0])
        col2.append(info[1])
        col3.append(info[2])
    df['n_detected'] = col1
    df['n_total'] = col2
    df['frac'] = col3

    df_every = pd.DataFrame({'hla':['all_sample'],'rank_pert':[None],'nM':[None],'id':[None],'freq':[None],'n_detected':[len(hlas)],'n_total':[len(meta['biology'].unique())],'frac':[round(len(hlas)/len(meta['biology'].unique()),2)]})
    df_all = df.loc[df['hla']=='all_hla',:]
    df_rest = df.loc[df['hla']!='all_hla',:].sort_values(by='frac',ascending=False)
    df = pd.concat([df_every,df_all,df_rest],axis=0)

    return df




def simplify_self_gene_source(source):
    pat = re.compile(r'ENSG(\d+)\|ENST(\d+)\|')
    sources = source.split(';')
    tmp = []
    for item in sources:
        match = re.search(pat,item)
        if match:
            tmp.append(item)
    gs_list = []
    for item in tmp:
        gs = item.split('|')[-1]
        gs_list.append(gs)
    return ';'.join(gs_list)

    

@callback(Output('output_header','children'),Output('output_text','children'),Input('candidate','active_cell'),Input('candidate', 'data'),Input('candidate','page_current'),Input('candidate','page_size'))
def click_table(active_cell,data,page_current,page_size):
    if active_cell:
        row = active_cell['row'] + page_current * page_size
        antigen_value = data[row]['Peptide']
        cell_value = data[row][active_cell['column_id']]
        return 'You selected antigen: {}'.format(antigen_value),'and clicked the cell with following content: {}'.format(cell_value)
    else:
        return 'You selected antigen: None', 'and clicked the cell with following content: None'


@callback(Output('psm','src'),Input('cancer_dropdown','value'),Input('candidate','active_cell'),Input('candidate', 'data'),Input('candidate','page_current'),Input('candidate','page_size'))
def draw_psm(cancer,active_cell,data,page_current,page_size):
    if active_cell:
        row = active_cell['row'] + page_current * page_size
        cell_value = data[row]['Peptide']  # must be pep
        pep = cell_value
        cancer = data[row]['Cancer']
        image_name = '{}_spectrum_{}.png'.format(cancer,pep)
        cmd = 'curl -H "Cache-Control: no-cache" -o {} -k {}'.format(os.path.join(assets_dir,image_name),os.path.join(cloud_dir,image_name))
        subprocess.run(cmd,shell=True)
        return app.get_asset_url('{}_spectrum_{}.png'.format(cancer,pep))
   

@callback(Output('differential_1','src'),Output('differential_2','src'),Output('intensity_1','src'),Output('intensity_2','src'),Input('cancer_dropdown','value'),Input('candidate','active_cell'),Input('candidate', 'data'),Input('candidate','page_current'),Input('candidate','page_size'))
def draw_diffential_and_intensity(cancer,active_cell,data,page_current,page_size):
    if active_cell:
        row = active_cell['row'] + page_current * page_size
        cell_value = data[row]['Peptide']  # must be pep
        pep = cell_value
        typ = data[row]['Antigen Class']
        cancer = data[row]['Cancer']
        source = data[row]['Source']
        source = modify_source_string(source)

        print(source)

        exist_percentile_plot_path = os.path.join(assets_dir,'{}_{}_{}.png'.format(cancer,pep,'percentile'))
        exist_rank_abundance_plot_path = os.path.join(assets_dir,'{}_{}_rank_abundance.png'.format(cancer,pep))

        exist_percentile_plot_file_name = '{}_{}_{}.png'.format(cancer,pep,'percentile')
        exist_rank_abundance_plot_file_name = '{}_{}_rank_abundance.png'.format(cancer,pep)

        if source == 'not_unique':
            typ = 'not_unique'

        if typ == 'self gene':
            ensg,enst,symbol = source.split('|')   
            exist_plot_gene_file_name =  '{}_{}_{}_expr_{}.png'.format(cancer,ensg,symbol,'boxplot+boxplot')
            cmd1 = 'curl -H "Cache-Control: no-cache" -o {} -k {}'.format(os.path.join(assets_dir,exist_plot_gene_file_name),os.path.join(cloud_dir,exist_plot_gene_file_name))
            cmd2 = 'curl -H "Cache-Control: no-cache" -o {} -k {}'.format(os.path.join(assets_dir,exist_percentile_plot_file_name),os.path.join(cloud_dir,exist_percentile_plot_file_name))
            cmd3 = 'curl -H "Cache-Control: no-cache" -o {} -k {}'.format(os.path.join(assets_dir,exist_rank_abundance_plot_file_name),os.path.join(cloud_dir,exist_rank_abundance_plot_file_name))
            for cmd in [cmd1,cmd2,cmd3]:
                subprocess.run(cmd,shell=True)
            return app.get_asset_url(exist_plot_gene_file_name),None,app.get_asset_url(exist_percentile_plot_file_name),app.get_asset_url(exist_rank_abundance_plot_file_name)

        elif typ == 'splicing':
            coords = source.split('|')[0]
            exist_plot_splicing_file_name = '{}_{}_splicing.png'.format(cancer,coords)
            cmd1 = 'curl -H "Cache-Control: no-cache" -o {} -k {}'.format(os.path.join(assets_dir,exist_plot_splicing_file_name),os.path.join(cloud_dir,exist_plot_splicing_file_name))
            cmd2 = 'curl -H "Cache-Control: no-cache" -o {} -k {}'.format(os.path.join(assets_dir,exist_percentile_plot_file_name),os.path.join(cloud_dir,exist_percentile_plot_file_name))
            cmd3 = 'curl -H "Cache-Control: no-cache" -o {} -k {}'.format(os.path.join(assets_dir,exist_rank_abundance_plot_file_name),os.path.join(cloud_dir,exist_rank_abundance_plot_file_name))
            for cmd in [cmd1,cmd2,cmd3]:
                subprocess.run(cmd,shell=True)
            return app.get_asset_url(exist_plot_splicing_file_name),None,app.get_asset_url(exist_percentile_plot_file_name),app.get_asset_url(exist_rank_abundance_plot_file_name)
        elif typ == 'TE chimeric transcript':
            coords = source.split('|')[0]
            erv = source.split('|')[7].split(',')[1]
            exist_plot_splicing_file_name = '{}_{}_splicing.png'.format(cancer,coords)
            exist_plot_erv_file_name = '{}_{}_expr.png'.format(cancer,erv)
            cmd1 = 'curl -H "Cache-Control: no-cache" -o {} -k {}'.format(os.path.join(assets_dir,exist_plot_splicing_file_name),os.path.join(cloud_dir,exist_plot_splicing_file_name))
            cmd2 = 'curl -H "Cache-Control: no-cache" -o {} -k {}'.format(os.path.join(assets_dir,exist_plot_erv_file_name),os.path.join(cloud_dir,exist_plot_erv_file_name))
            cmd3 = 'curl -H "Cache-Control: no-cache" -o {} -k {}'.format(os.path.join(assets_dir,exist_percentile_plot_file_name),os.path.join(cloud_dir,exist_percentile_plot_file_name))
            cmd4 = 'curl -H "Cache-Control: no-cache" -o {} -k {}'.format(os.path.join(assets_dir,exist_rank_abundance_plot_file_name),os.path.join(cloud_dir,exist_rank_abundance_plot_file_name))
            for cmd in [cmd1,cmd2,cmd3,cmd4]:
                subprocess.run(cmd,shell=True)
            return app.get_asset_url(exist_plot_splicing_file_name),app.get_asset_url(exist_plot_erv_file_name),app.get_asset_url(exist_percentile_plot_file_name),app.get_asset_url(exist_rank_abundance_plot_file_name)

        elif typ == 'ERV':
            erv = source.split('|')[0]  
            exist_plot_erv_file_name = '{}_{}_expr.png'.format(cancer,erv)
            cmd1 = 'curl -H "Cache-Control: no-cache" -o {} -k {}'.format(os.path.join(assets_dir,exist_plot_erv_file_name),os.path.join(cloud_dir,exist_plot_erv_file_name))
            cmd2 = 'curl -H "Cache-Control: no-cache" -o {} -k {}'.format(os.path.join(assets_dir,exist_percentile_plot_file_name),os.path.join(cloud_dir,exist_percentile_plot_file_name))
            cmd3 = 'curl -H "Cache-Control: no-cache" -o {} -k {}'.format(os.path.join(assets_dir,exist_rank_abundance_plot_file_name),os.path.join(cloud_dir,exist_rank_abundance_plot_file_name))
            for cmd in [cmd1,cmd2,cmd3]:
                subprocess.run(cmd,shell=True)
            return  app.get_asset_url(exist_plot_erv_file_name),None,app.get_asset_url(exist_percentile_plot_file_name),app.get_asset_url(exist_rank_abundance_plot_file_name)
        
        elif typ == 'fusion' or typ == 'variant' or typ == 'pathogen' or typ == 'nuORF' or typ == 'intron retention' or typ == 'not_unique':  # no diff plot
            cmd1 = 'curl -H "Cache-Control: no-cache" -o {} -k {}'.format(os.path.join(assets_dir,exist_percentile_plot_file_name),os.path.join(cloud_dir,exist_percentile_plot_file_name))
            cmd2 = 'curl -H "Cache-Control: no-cache" -o {} -k {}'.format(os.path.join(assets_dir,exist_rank_abundance_plot_file_name),os.path.join(cloud_dir,exist_rank_abundance_plot_file_name))
            for cmd in [cmd1,cmd2]:
                subprocess.run(cmd,shell=True)
            return  None,None,app.get_asset_url(exist_percentile_plot_file_name),app.get_asset_url(exist_rank_abundance_plot_file_name)

def modify_source_string(source):
    pat = re.compile(r'ENSG(\d+)\|ENST(\d+)\|')

    # first check out self_gene
    if ';' in source:
        sources = source.split(';')
        tmp = []
        for item in sources:
            match = re.search(pat,item)
            if match:
                tmp.append(item)
        if len(tmp) == 1:
            source = tmp[0]
        elif len(tmp) > 1:
            source = 'not_unique'
        else:
            source = source

    # if not self_gene, further get rid of nc and nuorf
    if ';' in source:  
        sources = source.split(';')
        tmp = []
        for item in sources:
            if ('nuORF' not in item) and ('nc|' not in item):
                tmp.append(item)
        source = ';'.join(tmp)

        if ';' in source:
            source = 'not_unique'

    return source


@callback(Output('hla_table','data'),Input('cancer_dropdown','value'),Input('candidate','active_cell'),Input('candidate', 'data'),Input('candidate','page_current'),Input('candidate','page_size'))   
def display_hla_table(cancer,active_cell,data,page_current,page_size):

    if active_cell:
        row = active_cell['row'] + page_current * page_size
        if cancer == 'All':
            cancer = data[row]['Cancer']
        meta = pd.read_csv('./static/{}_metadata.txt'.format(cancer),sep='\t')
        lists = literal_eval(data[row]['additional_query'])
        hlas = literal_eval(data[row]['presented_by_each_sample_hla'])
        df = get_hla_info(meta,lists,hlas)

        # add deepimmuno
        try:
            sub_mapping = deepimmuno_dic[data[row]['Peptide']]  # hla to immunogenicity
        except:
            df['immunogenicity'] = [None for hla in df['hla']]   # not predictable by deepimmuno
        else:
            df['immunogenicity'] = [sub_mapping.get(hla.replace(':',''),None) for hla in df['hla']]

        # rename
        df.columns = ['HLA','Rank(%)','Binding(nM)','Identity','Frequency','#Detected','#Total','Recurrence','Immunogenicity']


        return df.to_dict('records')

# @callback(Output('download_image', 'data'),Input('btn_download_img', 'n_clicks'), Input('psm','src'))
# def func(n_clicks,src):
#     return dcc.send_file('.' + src)


if __name__ == '__main__':

    us_hla = pd.read_csv('./static/US_HLA_frequency.csv',sep=',',index_col=0)
    us_hla['Percent US population'] = us_hla['Percent US population'].round(decimals=2)
    us_hla.index = hla_formatting(us_hla.index.to_list(),'deepimmuno_nostar','netMHCpan_input')
    hla2freq = us_hla['Percent US population'].to_dict()
    assets_dir = './assets'
    cloud_dir = 'https://genome.med.nyu.edu/public/yarmarkovichlab/ImmunoVerse/assets'

    # preprocess deepimmuno
    deepimmuno_df = pd.read_csv('./static/all_deepimmuno_immunogenicity.txt',sep='\t')
    deepimmuno_df['immunogenicity'] = deepimmuno_df['immunogenicity'].round(decimals=2)
    deepimmuno_dic = {}
    for pep,sub_df in deepimmuno_df.groupby(by='peptide'):
        mapping = sub_df.loc[:,['HLA','immunogenicity']].set_index(keys='HLA')['immunogenicity'].to_dict()
        deepimmuno_dic[pep] = mapping

    # start to build app
    app = Dash(__name__,assets_folder='./assets')
    app.title = 'ImmunoVerse'
    app.layout = html.Div([

        # # url
        # dcc.Location(id='url', refresh=False),

        # header
        html.Div([
        html.Div(html.Img(src='./static/logo.png', className='logo'),className='logo_div'),
        html.Div(children='Pan-Cancer Atlas of T Cell Therapeutic Targets',className='title_div')
        ],className='header_div'),
        
        # body
        html.Div([

            # query row
            html.Div([
                html.Div([html.Div('Cancer:',className='body_query_note'),dcc.Dropdown(id='cancer_dropdown',options=['All','BRCA','KIRC','COAD','STAD','MESO','LIHC','ESCA','CESC','BLCA','RT','AML','DLBC','GBM','NBL','PAAD','HNSC','OV','LUSC','LUAD','CHOL','SKCM'],value='NBL',className='body_query_func')],className='body_query_row'),
                html.Div([html.Div('Antigen Class:',className='body_query_note'),dcc.Dropdown(id='query_type',options=['All','self_gene','splicing','TE_chimeric','TE_self_translate','nuORF','intron_retention','fusion','variant','pathogen'],value='All',className='body_query_func')],className='body_query_row'),
                html.Div([html.Div('Peptide:',className='body_query_note'),dcc.Input(id='query_peptide',placeholder='QYNPIRTTF',className='body_query_func')],className='body_query_row'),
                html.Div([html.Div('Gene/Source:',className='body_query_note'),dcc.Input(id='query_source',placeholder='PHOX2B',className='body_query_func')],className='body_query_row'),
                html.Div([html.Div('HLA:',className='body_query_note'),dcc.Input(id='query_hla',placeholder='A*24:02',className='body_query_func')],className='body_query_row'),
                html.Div([html.Div('Sort by:',className='body_query_note'),dcc.Dropdown(id='sort_dropdown',options=['Peptide','Spectral Score','Abundance'],value='Spectral Score',className='body_query_func')],className='body_query_row'),
                html.Div([html.Div('Download Table:',className='body_query_note'),dcc.Dropdown(id='download_table',options=['Yes','No'],value='No',className='body_query_func'),dcc.Download(id='downloader')],className='body_query_row'),
                html.Button('Submit', id='submit_button', n_clicks=0, className='body_query_row'),
            ],className='body_query_div'),

            # display row
            html.Div([
                html.Div(html.H5('Please leave Peptide, Gene/Source, HLA blank if not specific (not case-sensitive), then click submit button'),style={'text-align':'left','color':'#3B5998'}),
                html.H5("Please click any row/cell in the table below to expand visuals for each antigen",style={'text-align': 'left','color':'#3B5998'}),
                html.H5("When selecting All cancers, please allow up to 10 seconds to render the table",style={'text-align': 'left','color':'#3B5998'}),
                html.H5(['Due to periodic HPC maintainence, please check back again if encountering rendering issues or ', html.A('Contact us',href='mailto:guangyuan.li@nyulangone.org',title='Click to send an email')],style={'text-align': 'left','color':'#3B5998'}),

                html.Div(dash_table.DataTable(id='candidate',
                                              page_size=10, page_current=0, page_action='native',hidden_columns=['presented_by_each_sample_hla','additional_query','Source'],
                                              tooltip_header = {
                                                 'Cancer':'Abbreviation for each cancer',
                                                 'Peptide':'Peptide Sequence',
                                                 'Antigen Class':'Class of Molecular Event',
                                                 'Simplified Source':'Source Aberration',
                                                 'Spectral Score':'Confidence of MS identification, the higher the better',
                                                 'Abundance':'Percentile in descending order, the higher the more abundant',
                                                 'Essentiality':'Depmap Chrono Score, the lower the more essential',
                                                 'Homogeneity':'Single Cell coverage in percentile, the higher the better, ranging from 0-1'
                                              },
                                              style_table={'overflowX': 'auto'},
                                              style_cell={
                                                'minWidth': '165px', 'width': '165px', 'maxWidth': '165px',
                                                'textAlign': 'left',
                                                'textOverflow': 'ellipsis',
                                                'fontFamily': 'Inter, sans-serif',
                                                'cursor': 'pointer'
                                                },

                                              style_header={
                                                'fontWeight': 'bold',
                                                'cursor':'help'
                                            
                                              },
                                              style_data_conditional=[
                                                {
                                                    'if': {'row_index': 'odd'},
                                                    'backgroundColor': 'rgb(242, 242, 242)',
                                                },
                                              ]
                                              ),className='main_table'),


                html.Div([html.H2(id='output_header'),html.P(id='output_text')]),

                html.Div([
                        html.H2('Differential plots'),
                        html.Img(id='differential_1',width='45%',style={'float':'left'}),
                        html.Img(id='differential_2',width='45%',style={'float':'right'})
                    ],style={'overflow':'hidden'}),

                html.Div([html.H2('PSM plot'),html.Img(id='psm',width='45%')]),

                html.Div([
                        html.H2('Intensity plots'),
                        html.Img(id='intensity_1',width='45%',style={'float':'left'}),
                        html.Img(id='intensity_2',width='45%',style={'float':'right'}),
                    ],style={'overflow':'hidden'}),

                html.Div(html.H2('HLA table')),
                html.Div(dash_table.DataTable(id='hla_table',
                                            page_size=20,page_current=0,page_action='native',   
                                            tooltip_header = {
                                                'HLA':'HLA Allele',
                                                'Rank(%)':'NetMHCpan4.1 ranking percentile',
                                                'Binding(nM)':'NetMHCpan4.1 predicted binding (nM)',
                                                'Identity':'NetMHCpan4.1 predicted Weak Binder (WB) or Strong Binder (SB)',
                                                'Frequency':'U.S. Population coverage',
                                                '#Detected':'Number of sample where the peptide gets identified',
                                                '#Total':'Total number of sample possessing this HLA',
                                                'Recurrence':'Recurrence of the peptide',
                                                "Immunogenicity":'DeepImmuno predicted immunogenicity'
                                              },
                                              style_table={'overflowX': 'auto'},
                                              style_cell={
                                                'minWidth': '165px', 'width': '165px', 'maxWidth': '165px',
                                                'textAlign': 'left',
                                                'textOverflow': 'ellipsis',
                                                'fontFamily': 'Inter, sans-serif',
                                                },
                                              style_header={
                                                'fontWeight': 'bold',
                                                'cursor':'help'
                                              },
                                              style_data_conditional=[
                                                {
                                                    'if': {'row_index': 'odd'},
                                                    'backgroundColor': 'rgb(242, 242, 242)',
                                                }]
                                            )),


                
                ],className='body_display_div')



            ],className='body_div'),

        # stats
        html.Div([
            html.Div([
                html.Div([
                    html.H3("27,541", className="stat-number"),
                    html.P("Tumor-specific Antigens", className="stat-label")
                ], className="stat-box"),

                html.Div([
                    html.H3("21", className="stat-number"),
                    html.P("Cancer Types", className="stat-label")
                ], className="stat-box"),

                html.Div([
                    html.H3("11", className="stat-number"),
                    html.P("Molecular Sources", className="stat-label")
                ], className="stat-box"),
            ], className="stats-row"),

            html.Div([
                html.A("Read the Research Paper", href="https://www.biorxiv.org/content/10.1101/2025.01.22.634237v1", className="cta-link", target="_blank")
            ], className="stats-cta")
        ], className="stats-section"),

        
        # footer
        html.Div([
            # Row 1: Full-width branding
            html.Div([
                html.H2("ImmunoVerse", style={"color": "white"}),
                html.P(
                    "A comprehensive atlas of tumor-specific T cell targets across human cancers, providing resources for immunotherapy development.",
                    style={"color": "lightgray", "maxWidth": "600px"}
                )
            ], className="footer-top"),

            # Row 2: Side-by-side columns
            html.Div([
                html.Div([
                    html.H4("RESOURCES"),
                    html.A("Documentation", href="https://github.com/frankligy/pan_cancer_intracellular_antigen_atlas", className="footer-link", target="_blank"),
                    html.Br(),
                    html.A("API Reference", href="https://github.com/frankligy/pan_cancer_intracellular_antigen_atlas", className="footer-link", target="_blank")
                ], className="footer-column"),

                html.Div([
                    html.H4("CONNECT"),
                    html.A("GitHub", href="https://github.com/frankligy/pan_cancer_intracellular_antigen_atlas", className="footer-link", target="_blank"),
                    html.Br(),
                    html.A("Contact Us", href="https://github.com/frankligy/pan_cancer_intracellular_antigen_atlas", className="footer-link", target="_blank")
                ], className="footer-column")
            ], className="footer-middle"),

            # Row 3: Bottom copyright
            html.Div("© 2025 ImmunoVerse. All rights reserved.", className="footer-bottom")
        ], className="footer")
        



    ],className='whole_div')

    host = subprocess.run(['hostname'],stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[0]
    port = 8050
    app.run(host=host,port=port)



