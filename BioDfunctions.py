import base64
import os
import json
import pickle
import uuid
import re
import requests
import time
import warnings
import io
import math

from glob import glob
from fnmatch import fnmatch
from scipy.optimize import curve_fit
from scipy.stats import linregress
from scipy.integrate import quad

import plotly.graph_objs as go
import streamlit as st
import numpy as np
import pandas as pd
from PIL import Image

vault_id=5938
tumor_spheres_file = pd.read_csv("spheres-310g-tumorDoses.csv")
irdc_version = "IRDC v1.1"

#############################################################
# variable definitions used for dosimetry:
#############################################################
time_unit = 'hours'

isotopes_halflives = {"Lu-177": 159.53,
                    "Ac-225": 238.08,
                    "I-125": 1427.76,
                    "I-131": 192.47,
                    "Ga-68": 1.13,
                    "Cu-64": 12.7}
                    
isotopes_BioD_halflives = {"Lu-177": 159.53,
                        "Ac-225": 238.08,
                        "Pb-203": 51.87,                                  
                        "Pb-212": 10.64,                                  
                        "I-125": 1427.76,
                        "I-131": 192.47,
                        "Ga-68": 1.13,
                        "Cu-64": 12.7}

isotopes_human_halflives = {"Lu-177": 159.6,          #halflives in hours   
                        "Ac-225": 238.08,   
                        "In-111": 67.2,                               
                        "Pb-203": 51.87,                                  
                        "Pb-212": 10.64,                                  
                        "Ra-223": 11.43*24,                                  
                        "Y-90": 64.1,                                  
                        "I-124": 100.32,                                    
                        "I-131": 192.48,                                    
                        "Bi-213": 45.59/60,                                    
                        "Ga-68": 1.13,                                     
                        "At-211": 7.21,                                      
                        "F-18": 1.83,                                      
                        "Zr-89": 78.41,                                      
                        "Cu-67": 61.83,                                      
                        "Cu-64": 12.7}  

doselimits_file = pd.read_csv("organ-doselimits-calculator.csv")
tissue_masses_mouse_file = pd.read_csv("tissue-masses-mouse.csv")
tissue_masses_human_file = pd.read_csv("tissue-masses-human.csv")
wholebody_mouse = 25. # gram 
mouse_keyword_tissues = 'Available tissues for mouse'
mouse_keyword_masses = 'Mass 30g mouse model [g]'  
human_keyword_tissues = 'Available tissues for human'
human_keyword_masses = 'Mass Male human model [g]'  
wholebody_human = 70. # kilogram   
wholebody_human = float(list(tissue_masses_human_file[human_keyword_masses][tissue_masses_human_file[human_keyword_tissues] == 'total body'])[0])
key_human_organ_weight = 'human organ weight [g]'
key_mouse_organ_weight = 'mouse organ weight [g]'
key_human_wb_weight = 'wb_human [g]'
key_mouse_wb_weight = 'wb_mouse [g]'
htiac_key = 'hTIAC_g [MBq * h/MBq]'
htiac_org_key = 'hTIAC_org [MBq * h/MBq]'
keyword_mtiac_g = 'mTIAC_g [MBq * h/MBq]'
keyword_mtiac_org = 'mTIAC_org [MBq * h/MBq]'
keyword_Radioisotope1 = 'radioisotope 1 (used in BioD)'
time_keyword = 'Time'
injdose_keyword = 'id/g '
keyword_projected_bm = 'bone marrow (from blood)'


alpha_scal= 'metabolic scaling (2) (alpha weighted scaling - Stephen Graves)' 
time_mass_fda = 'combined time and relative mass scaling - FDA guideline'
rel_mass_scal = 'relative mass scaling (FDA)'
rel_mass_scal_fda = 'relative mass scaling - FDA guideline'
no_scaling = 'no scaling'
time_scal = 'time scaling'
time_mass_scal = 'combined time and relative mass scaling'
allom_scal = 'allometric scaling'
metabol_scal = 'metabolic scaling (1) with exponent 0.25'
scaling_options = alpha_scal, metabol_scal, rel_mass_scal, time_scal, time_mass_scal, no_scaling

monofit, bifit, biexpelim, biexpelim2, trapfit, linexpfit, linphysdecay = ['monoexponential - linear regression on log(y)','biexponential','biexponential with absorption and elimination (clearance)','biexponential with absorption and elimination (A1)','trapezoidal','linear extrapolation','trapezoidal with physical decay']
monoexp_dec_fit = 'monoexponential direct method'
biexp_abs_elim_two_components = 'biexponential with absorption and elimination (A1 and A2)'
fitmodel_options = [monoexp_dec_fit, monofit, bifit, biexpelim2, biexp_abs_elim_two_components, trapfit, linexpfit,linphysdecay] # excluded biexpelim from David Huang

bm_sgourus = 'Sgourus (1993)'
bm_sgourus_greg = 'Sgourus (1993) - Greg'
bm_traino = 'Traino (2007)'
bm_grudzinski = 'Joe Grudzinski'
bm_graves = 'Stephen Graves'
bm_options = [bm_traino, bm_graves, bm_sgourus, bm_sgourus_greg, bm_grudzinski]


mta_keyword = 'MTA [GBq]'
mta_organ_keyword = 'MTA per organ [GBq]'
doselimit_keyword = 'Dose limit [Gy]'
keyword_dose = 'absorbed dose [mSv/MBq = mGy/MBq]'
key_olindainput = 'Olinda Input Organ Name'
keyword_tumor_dose_mta = 'tumor exposure at MTA [Gy]'



#################################################################################
#####  Functions needed for dosimetry projections:
#################################################################################

def update_molecule_request():
    st.session_state.searchBioD = False
    st.session_state.tissues_updated = False      
    st.session_state.BioD_bigDataSet = False        
    st.session_state.calc_scaling = False
    st.session_state.continueDosimetry = False
    st.session_state.continue_with_dosimetry = False
    st.session_state.first_molecule_requested = True
    st.session_state.moleculefound = False

def biod_input_file(df_input):
    dragdropinput = False
    data_input = []
    keywords_list = []
    tissues = []
    injdose_keyword = ''
    st.header('Data input from file')
    tissues_all_options = (tissue_masses_mouse_file['Available tissues for mouse'].to_numpy())
    keyword_donotimport = 'do not import'
    column_options = [keyword_donotimport,time_keyword]
    for element in tissues_all_options:
        column_options.append(element)
    data_input = {}
    tissues = []
    input_expander = st.expander('See dataset rawdata and assign column header names')
    input_expander.dataframe(df_input)
    original_headers = input_expander.checkbox("Use original header names for all tissues?", value=True)
    for i_inputelement,element in enumerate(df_input.head()):
        if i_inputelement ==0 :
            assigned_keyword = [time_keyword]
            column_options_extended = column_options
        else:
            element_short = element.split(' ')[0]
            for writing in [element,element_short]:
                tissue_all_spellings = get_all_possible_tissue_spellings(writing)
                assigned_keyword = list(set(tissue_all_spellings).intersection(set(column_options)))
                if len(assigned_keyword) != 0:
                    break
            column_options_extended = column_options
            column_options_extended.append(element)                    

        if original_headers:
            try:
                preselect_index = column_options.index(element)
            except:
                preselect_index = column_options.index(assigned_keyword[0])
            selected_keyword = input_expander.selectbox(label = f'{element} is used as header',options=column_options,index = preselect_index, key=str(i_inputelement)+'_column_dragdrop')                
        elif len(assigned_keyword) == 0:
            preselect_index = column_options.index(keyword_donotimport)
            selected_keyword = input_expander.selectbox(label = f'Assign column names: {element} was not found in tissue list',options=column_options,index = preselect_index, key=str(i_inputelement)+'_column_dragdrop')
        else:
            preselect_index = column_options.index(assigned_keyword[0])
            selected_keyword = input_expander.selectbox(label = f'Assign column names: {element} is {assigned_keyword}, please change if not',options=column_options,index = preselect_index, key=str(i_inputelement)+'_column_dragdrop')

        if selected_keyword == keyword_donotimport:
            pass
        elif not selected_keyword == time_keyword:
            injd_keyword = injdose_keyword+ selected_keyword
            tissues.append(selected_keyword)
            try:
                data_input[injd_keyword] = df_input[element].map(lambda element: element.split('±')[0])
            except:
                data_input[injd_keyword] = df_input[element]
        else:
            data_input[selected_keyword] = df_input[element]

    keywords_list = []
    for ie,tissue_ie in enumerate(tissues):
        keyword = injdose_keyword+tissue_ie
        keywords_list.append(keyword)
    if len(keywords_list) > 0:
        dragdropinput = True
    else:
        dragdropinput = False


    # get only numeric time value if column Time is in format '2h, 24h, etc'
    try:
        data_input[time_keyword] = data_input[time_keyword].str.split('h').str[0]
    except:
        pass

    return [data_input,tissues,keywords_list,dragdropinput]

def biod_input_manually(submit_biod_input_manual):
    data_input = []
    keywords_list = []
    tissues = []
    with st.expander("Manual input of %ID/g"):
        nrtimepoints = st.slider('Number of timepoints', min_value=2, max_value=12, step=1,value=4)
        tissues = get_selected_tissues_list(tissue_masses_mouse_file)

        timepoints_initial  = [0.5, 1,  2,  4, 6, 10] # values from Stabin, The Practice of Internal Dosimetry, page 88
        injdose_initial     = [72, 35, 24, 20, 15, 12]
        with st.form(key='columns_in_form'):
            data_input = {}

            keywords_list = []
            data_input[time_keyword] = list()
            for ie,tissue_ie in enumerate(tissues):
                keyword = tissue_ie
                data_input[keyword] = list()
                keywords_list.append(keyword)

            columns = st.columns(len(data_input))
            for i in range(int(nrtimepoints)):
                for j, (k, v) in enumerate(data_input.items()):
                    with columns[j]:
                        if k == time_keyword:
                            timepoint = timepoints_initial[i] if i < len(timepoints_initial) else timepoints_initial[-1]
                            v.append(st.text_input(f'timepoint {i+1} in {time_unit} for each tissue [{time_unit}]',timepoint))
                        else:
                            for ie,tissue_ie in enumerate(tissues):
                                k = k.replace(injdose_keyword, '')
                                if fnmatch(k, tissue_ie):
                                    injd = injdose_initial[i] if i < len(injdose_initial) else injdose_initial[-1]
                                    injd += (ie)*(0.1*injd)
                                    v.append(st.text_input(f'{k} [%ID/g] at timepoint {i+1}',injd))     
            submitted = st.form_submit_button("Submit")
            if submitted:
                submit_biod_input_manual = True

    return [data_input,tissues,keywords_list,submit_biod_input_manual]

def search_moleculebatchID_CDD(continueCDDsearch,moleculefound,mol_batch_id="RB-0003401-006"):
    with st.form("Search CDD Vault with Molecule Batch-ID"):
        molecule_batch_ID = str(st.text_input('Search CDD Vault with Molecule Batch-ID', mol_batch_id)).lstrip()
        molecule_submitted = st.form_submit_button(f'Search in CDD Vault')
        if molecule_submitted:
            continueCDDsearch = True
            if 'RB' in molecule_batch_ID:
                molecule_batch_ID=molecule_batch_ID.split('RB')
                molecule_batch_ID='RB'+molecule_batch_ID[-1][:12]
                if len(molecule_batch_ID) != 14:
                    st.error('error: exactly RB-00000-000 format is needed')
                    continueCDDsearch = False
                    moleculefound = False
                else:
                    continueCDDsearch = True
                    update_molecule_request()

            else:
                st.error('error: no RB number')
                continueCDDsearch = False
                moleculefound = False
    return [molecule_batch_ID,continueCDDsearch,moleculefound]

def reset_session_state_dragndropInput():
    st.session_state.dragdropinput = False    
    st.session_state.calc_scaling = False
    st.session_state.decay_correction = False

def change_tissue_input(all_tissues,selected_tissues):
    with st.form('Select tissues'):
        tissue_checkbox = {}
        tissue_columns = st.columns(3)
        for i, tissue in enumerate(all_tissues):
            with tissue_columns[i%3]:
                if tissue in selected_tissues:
                    tissue_checkbox[tissue] = st.checkbox(tissue, value=True)
                else:
                    tissue_checkbox[tissue] = st.checkbox(tissue, value=False)
        tissues = [t for t in tissue_checkbox.keys() if tissue_checkbox[t]]

        submit_tissues = st.form_submit_button('Select tissues and submit here')
        if submit_tissues:
            st.session_state.tissues_updated = True
    return [st.session_state.tissues_updated,tissues]

def decay_correction_changed():
    if st.session_state.decay_correction == False:
        st.session_state.decay_correction = True

def extract_values_dict(dictionary):
    try:
        value = dictionary.get('value',None)
        note = dictionary.get('note',None)
        outlier = dictionary.get('outlier',None)
        outlier_type = dictionary.get('outlier_type',None)
        return value, note, outlier, outlier_type
    except:
        pass

def get_selected_tissues_list(tissue_masses_mouse_file):
    tissues_list = ['kidney','tumor', 'blood','muscle']
    selectable_tissues = np.asarray(tissue_masses_mouse_file['Available tissues for mouse'])
    selected_tissues = st.multiselect('add more tissues',selectable_tissues)
    tissues_list.extend(xtissue for xtissue in selected_tissues if xtissue not in tissues_list)
    tissue_checkbox = {}
    tissue_columns = st.columns(3)
    for i, tissue in enumerate(tissues_list):
        with tissue_columns[i%3]:
            if tissue in ["tumor", "kidney"]:
                tissue_checkbox[tissue] = st.checkbox(tissue, value=True)
            else:
                tissue_checkbox[tissue] = st.checkbox(tissue, value=False)
    tissues = [t for t in tissue_checkbox.keys() if tissue_checkbox[t]]
    # to make blood the second keyword
    if 'blood' in tissues:
        if 'kidney' in tissues:
            order = [tissues.index('kidney')]
        else:
            order = list()
        order.append(tissues.index('blood'))
        for element_keywords in tissues:
            if element_keywords == 'kidney':
                pass
            elif element_keywords != 'blood':
                order.append(tissues.index(element_keywords))
        tissues = [tissues[i_order] for i_order in order]
    return tissues

def select_conditions_BioD(parameters_experiment_BioD,readout_definitions,bioD_data):
    with st.form('Select conditions'):
        selected_conditions_all = dict()
        parameter_columns = st.columns(len(parameters_experiment_BioD))
        for i, parameter_name in enumerate(parameters_experiment_BioD):
            selected_conditions = dict()
            ignore_parameter = False
            with parameter_columns[i]:
                parameter_name_id = get_readout_name_id(parameter_name,readout_definitions)
                all_parameters_experiment = get_parameters_fromBioD(parameter_name_id,bioD_data)
                strings_all_parameter_readouts = [str(x) for x in all_parameters_experiment] 
                try:    
                    strings_all_parameter_readouts = sorted(strings_all_parameter_readouts, key=float)
                    all_parameters_experiment = sorted(all_parameters_experiment, key=float)
                except:
                    pass
                if fnmatch(parameter_name,'Time point'):
                    st.metric(parameter_name+' [hr]', ' / '.join(strings_all_parameter_readouts))
                else:
                    st.metric(parameter_name, ' / '.join(strings_all_parameter_readouts))
                for c_nr,condition_name in enumerate(all_parameters_experiment):
                    if condition_name == 'not defined':
                        st.write('not defined')
                        ignore_parameter = True
                    else:
                        try:
                            condition_name = float(condition_name)
                            selected_conditions[condition_name] = st.checkbox(f'{condition_name:.2f}', value=True, key=f'{c_nr}-{condition_name}')
                        except:
                            selected_conditions[condition_name] = st.checkbox(f'{condition_name}', value=True)
                if not ignore_parameter:
                    selected_conditions_all[parameter_name] = selected_conditions
        st.form_submit_button('Select conditions and submit here')
    return selected_conditions_all


def request_molecule_batch_ID():
    col1, col2 = st.columns([1,1])
    with col1:
        molecule_ID_entered = st.checkbox('Search CDD Vault with Molecule Batch-ID')
        if molecule_ID_entered: 
            molecule_batch_ID = st.text_input('Molecule Batch ID', "RB-0001419-001")

    with col2:
        rayz_ID_entered = st.checkbox('Search CDD Vault with RAYZ-ID and Batch-ID')
        if rayz_ID_entered:
            rayz_id = st.text_input('RAYZ ID', "RAYZ-06116-Lu")
            try:
                batches_rayz = batches_from_rayzID(rayz_id)

                st.table(data = pd.DataFrame(pd.DataFrame(batches_rayz)))
                batch_name = st.text_input('Batch name', "001")
                name_index = batches_rayz['Batch names'].index(batch_name)
                batch_id = batches_rayz['Batch IDs'][name_index]  
                molecule_batch_ID = batches_rayz['Molecule Batch-ID'][name_index]

            except:
                st.write('Please enter a valid ID for the search')

    if molecule_ID_entered:
        st.write('You chose: ', molecule_batch_ID)
        return [molecule_batch_ID,molecule_ID_entered,rayz_ID_entered]
    elif rayz_ID_entered:
        try:
            st.write(f'You chose {molecule_batch_ID}: {rayz_id} with Batch ID {batch_id}')
        except:
            st.write(f'Please enter both RAYZ ID and Batch ID')
        return [molecule_batch_ID,molecule_ID_entered,rayz_ID_entered]

    else:
        return ['n/a',molecule_ID_entered,rayz_ID_entered]    

def download_button(object_to_download, download_filename, button_text, pickle_it=False):
    """
    Generates a link to download the given object_to_download.
    Params:
    ------
    object_to_download:  The object to be downloaded.
    download_filename (str): filename and extension of file. e.g. mydata.csv,
    some_txt_output.txt download_link_text (str): Text to display for download
    link.
    button_text (str): Text to display on download button (e.g. 'click here to download file')
    pickle_it (bool): If True, pickle file.
    Returns:
    -------
    (str): the anchor tag to download object_to_download
    Examples:
    --------
    download_link(your_df, 'YOUR_DF.csv', 'Click to download data!')
    download_link(your_str, 'YOUR_STRING.txt', 'Click to download text!')
    """
    if pickle_it:
        try:
            object_to_download = pickle.dumps(object_to_download)
        except pickle.PicklingError as e:
            st.write(e)
            return None

    else:
        if isinstance(object_to_download, bytes):
            pass

        elif isinstance(object_to_download, pd.DataFrame):
            object_to_download = object_to_download.to_csv(index=False)

        # Try JSON encode for everything else
        else:
            object_to_download = json.dumps(object_to_download)

    try:
        # some strings <-> bytes conversions necessary here
        b64 = base64.b64encode(object_to_download.encode()).decode()

    except AttributeError as e:
        b64 = base64.b64encode(object_to_download).decode()

    button_uuid = str(uuid.uuid4()).replace('-', '')
    button_id = re.sub('\d+', '', button_uuid)

    custom_css = f""" 
        <style>
            #{button_id} {{
                background-color: rgb(17, 119, 187);
                color: rgb(255, 255, 255);
                padding: 0.25em 0.38em;
                position: relative;
                text-decoration: none;
                border-radius: 4px;
                border-width: 1px;
                border-style: solid;
                border-color: rgb(230, 234, 241);
                border-image: initial;
            }} 
            #{button_id}:hover {{
                border-color: rgb(246, 51, 102);
                color: rgb(30, 30, 30);
            }}
            #{button_id}:active {{
                box-shadow: none;
                background-color: rgb(246, 51, 102);
                color: white;
                }}
        </style> """

    dl_link = custom_css + f'<a download="{download_filename}" id="{button_id}" href="data:file/txt;base64,{b64}">{button_text}</a><br></br>'
    return dl_link

def find_molecule_batch_id(molecule_id,token,bigDataSet=True):
    if bigDataSet:
        export_id, identified_molecule = get_molecule_batch_id_bigDataset(molecule_id, token)
        pid = os.getpid()
        return [export_id, identified_molecule,pid]
    else:
        base_url = f"https://app.collaborativedrug.com/api/v1/vaults/{vault_id}/"
        headers = {'X-CDD-token':f'{token}'}
        url = base_url + f"batches"
        parameters_search = {'molecule_batch_identifier': molecule_id}
        response = requests.request("GET", url, headers=headers, params=parameters_search).json()
        return response

def find_molecule_batch_id_bigDataset(molecule_id,token):
    base_url = f"https://app.collaborativedrug.com/api/v1/vaults/{vault_id}/"
    headers = {'X-CDD-token':f'{token}'}
    url = base_url + f"batches"
    parameters_search = {'molecule_batch_identifier': molecule_id, 'async': True}
    response = requests.request("GET", url, headers=headers, params=parameters_search).json()
    try:
        return response["id"]
    except:
        st.write(response)

def get_molecule_batch_id_bigDataset(molecule_id, token):
    export_id = find_molecule_batch_id_bigDataset(molecule_id,token)
    with st.spinner('Wait for database query'):
        i = 0
        status = "new"
        while True:
            print(f"Export status is {status}, checking in {2**i} seconds...")
            if i < 5:
                time.sleep(2**i)
            else:
                time.sleep(2*i)
            status = check_export_status(export_id, token)
            if status == "finished":
                print("Export ready!")
                break
            i += 1
    identified_molecule = get_export(export_id, token)
    return [export_id, identified_molecule]

def get_single_molecule(molecule_id,token):
    base_url = f"https://app.collaborativedrug.com/api/v1/vaults/{vault_id}/"
    headers = {'X-CDD-token':f'{token}'}
    url = base_url + f"molecules/{molecule_id}"
    response = requests.request("GET", url, headers=headers).json()
    return response

def get_molecule_batch_name(molecule_id,batch_id,token):
    molecule_request = get_single_molecule(molecule_id,token)
    for element in molecule_request['batches']:
        if element['id'] == batch_id:
            molecule_batch_name = element['molecule_batch_identifier']
            break
        else:
            molecule_batch_name = 'n/a'
    return molecule_batch_name

def check_export_status(export_id, token):
    base_url = f"https://app.collaborativedrug.com/api/v1/vaults/{vault_id}/"
    headers = {'X-CDD-token':f'{token}'}
    url = base_url + f"export_progress/{export_id}"

    response = requests.request("GET", url, headers=headers).json()
    return response["status"]

def get_BioD_data_bigDataset(molecule_id,batch_id,token,BioD_id = 57565):
    base_url = f"https://app.collaborativedrug.com/api/v1/vaults/{vault_id}/"
    headers = {'X-CDD-token':f'{token}'}
    url = base_url + f"protocols/{BioD_id}/data"
    parameters_bioD_request = {'molecules':molecule_id, 'async': True, 'page_size': 1000}
    bioD_request = requests.request("GET", url, headers=headers, params=parameters_bioD_request).json()
    return bioD_request["id"]

def get_molecules_from_BioD_bigDataset(molecule_id,batch_id, token, BioD_id = 57565):
    export_id = get_BioD_data_bigDataset(molecule_id,batch_id, token, BioD_id)
    with st.spinner('Wait for database query'):
        i = 0
        status = "new"
        while True:
            print(f"Export status is {status}, checking in {2**i} seconds...")
            time.sleep(2**i)
            status = check_export_status(export_id, token)
            if status == "finished":
                print("Export ready!")
                break
            i += 1
    identified_bioD = get_export(export_id, token)
    bioD_data = list(filter(lambda item: item['molecule'] == molecule_id and item['batch'] == batch_id, identified_bioD['objects']))
    return [export_id, bioD_data]

def get_export(export_id, token):
    base_url = f"https://app.collaborativedrug.com/api/v1/vaults/{vault_id}/"
    headers = {'X-CDD-token':f'{token}'}
    url = base_url + f"exports/{export_id}"

    response = requests.request("GET", url, headers=headers)
    return response.json()

def protocol_molecules_query(protocol_id, token):
    base_url = f"https://app.collaborativedrug.com/api/v1/vaults/{vault_id}/"
    headers = {'X-CDD-token':f'{token}'}
    url = base_url + f"protocols/{protocol_id}/data"
    parameters_bioD_request = {'async': True, 'page_size': 1000}
    protocol_request = requests.request("GET", url, headers=headers, params=parameters_bioD_request).json()
    return protocol_request["id"]

def get_molecules_from_protocol(protocol_id, token):
    export_id = protocol_molecules_query(protocol_id, token)
    with st.spinner('Wait for database query'):
        i = 0
        status = "new"
        while True:
            print(f"Export status is {status}, checking in {2**i} seconds...")
            time.sleep(2**i)
            status = check_export_status(export_id, token)
            if status == "finished":
                print("Export ready!")
                break
            i += 1
    return get_export(export_id, token)

def get_all_uploaded_BioD(token):
    BioD_id = 57565
    BioD_request = get_molecules_from_protocol(BioD_id, token)
    BioD_columns_list = BioD_request.columns.values

    BioD_results_molecule_batch_ids = list()
    new_biod_id = 0
    new_batch_id = 0
    BioD_molecule = 0
    BioD_batch_id = 0

    all_molecule_ids=[]
    for item_nr,item in enumerate(BioD_columns_list):
        item = item.split(':')
        id = item[-1].split('.')[0]
        if fnmatch(item[0], '*molecule*'):
            new_biod_id = int(id)
            all_molecule_ids.append(new_biod_id)
        elif fnmatch(item[0], '*batch*'):
            new_batch_id = int(id)

        if int(new_biod_id) != int(BioD_molecule):
            BioD_molecule = new_biod_id
        if int(new_batch_id) != int(BioD_batch_id):
            BioD_batch_id = new_batch_id   
            mol_batch_name = get_molecule_batch_name(new_biod_id,new_batch_id,token)

            # get RAYZ ID
            # molecule_request = find_molecule_batch_id(mol_batch_name,token)
            # batch_id = molecule_request['objects'][0]['id']  
            # batch_information = get_batch_information(batch_id,token)
            # RAYZ_internal_ID = batch_information['molecule_fields']['RAYZ Internal ID']

            # BioD_results_molecule_batch_ids[mol_batch_name] = RAYZ_internal_ID
            BioD_results_molecule_batch_ids.append(mol_batch_name)

    # results = {
    # 'Molecule-Batch-ID': BioD_results_molecule_batch_ids
    #     }
    # results_df = pd.DataFrame(results)
    without_duplicates = [*set(BioD_results_molecule_batch_ids)]
    return without_duplicates 

def get_uploaded_BioD(token):
    base_url = f"https://app.collaborativedrug.com/api/v1/vaults/{vault_id}/"
    headers = {'X-CDD-token':f'{token}'}
    BioD_id = 57565
    url = base_url + f"protocols/{BioD_id}/data"
    parameters_bioD_request = {'async': False, 'page_size': 1000}
    bioD_request = requests.request("GET", url, headers=headers, params=parameters_bioD_request).json()

    BioD_results_molecule_batch_ids = list()
    batch_id_list = list()
    BioD_molecule = 0
    BioD_batch_id = 0
    for item_nr,item in enumerate(bioD_request['objects']):
        new_biod_id = item['molecule']
        new_batch_id = item['batch']
        if int(new_biod_id) != int(BioD_molecule):
            BioD_molecule = new_biod_id
        if int(new_batch_id) != int(BioD_batch_id):
            BioD_batch_id = new_batch_id   
            mol_batch_name = get_molecule_batch_name(new_biod_id,new_batch_id,token)

            # get RAYZ ID
            molecule_request = find_molecule_batch_id(mol_batch_name,token)
            # batch_id = molecule_request['objects'][0]['id']  
            batch_id = item_nr 
            # batch_information = get_batch_information(batch_id,token)
            # RAYZ_internal_ID = batch_information['molecule_fields']['RAYZ Internal ID']

            # BioD_results_molecule_batch_ids[mol_batch_name] = RAYZ_internal_ID
            BioD_results_molecule_batch_ids.append(mol_batch_name)
            batch_id_list.append(batch_id)
            break # to only return 1 molecule

    results = {
    'Molecule-Batch-ID': BioD_results_molecule_batch_ids
        }
    results_df = pd.DataFrame(results)
    return results_df 

def get_batch_information(batch_id,token):
    base_url = f"https://app.collaborativedrug.com/api/v1/vaults/{vault_id}/"
    headers = {'X-CDD-token':f'{token}'}
    url = base_url + f"batches/{batch_id}"
    parameters_request = {'async': False, 'page_size': 1000}
    request = requests.request("GET", url, headers=headers).json()
    return request['molecule']

def batches_from_rayzID(rayz_id,token):
    base_url = f"https://app.collaborativedrug.com/api/v1/vaults/{vault_id}/"
    headers = {'X-CDD-token':f'{token}'}
    url = base_url + f"molecules"
    parameters_request = {'names': rayz_id}
    request = requests.request("GET", url, headers=headers, params=parameters_request).json()
    batches_from_rayzID = dict()
    batchnames_list = []
    batchIDs_list = []
    moleculebatchIDs_list = []
    rayz_id_objects = request['objects']
    for batch in rayz_id_objects[0]['batches']:
        batchnames_list.append(batch['name'])
        batchIDs_list.append(batch['batch_fields']['Batch ID'])
        moleculebatchIDs_list.append(batch['molecule_batch_identifier'])

    batches_from_rayzID['Batch names'] = batchnames_list
    batches_from_rayzID['Batch IDs'] = batchIDs_list
    batches_from_rayzID['Molecule Batch-ID'] = moleculebatchIDs_list

    return batches_from_rayzID

def check_for_BioD_data(molecule_id,token):
    base_url = f"https://app.collaborativedrug.com/api/v1/vaults/{vault_id}/"
    headers = {'X-CDD-token':f'{token}'}
    BioD_id = 57565
    url = base_url + f"protocols/{BioD_id}/data"
    parameters_bioD_request = {'async': False, 'page_size': 1000}
    bioD_request = requests.request("GET", url, headers=headers, params=parameters_bioD_request).json()

    BioD_molecules_list = []
    for object in bioD_request['objects']:
        molecule_ids_inBioD = object['molecule']
        BioD_molecules_list.append(molecule_ids_inBioD)

    if molecule_id in BioD_molecules_list:
        return True
    else:
        return False

def get_BioD_data(molecule_id,batch_id,token):
    base_url = f"https://app.collaborativedrug.com/api/v1/vaults/{vault_id}/"
    headers = {'X-CDD-token':f'{token}'}
    BioD_id = 57565
    url = base_url + f"protocols/{BioD_id}/data"
    parameters_bioD_request = {'molecules':molecule_id, 'async': False, 'page_size': 1000}
    bioD_request = requests.request("GET", url, headers=headers, params=parameters_bioD_request).json()

    bioD_request_objects = bioD_request['objects']
    identified_bioD = list(filter(lambda item: item['molecule'] == molecule_id and item['batch'] == batch_id, bioD_request_objects))
    return identified_bioD

def get_readout_name_id(readout_name,readout_definitions):
    try:
        readout_name_id = str(list(filter(lambda item: item['name'] == readout_name, readout_definitions))[0]['id'])
        return readout_name_id
    except:
        st.write(readout_name)

def get_tissue_and_time(tissue, timepoint, bioD_data,parameter_tissue_id, parameter_timepoint_id): #standard values for BioD from cut&count
    alldata = bioD_data[bioD_data[f'{parameter_tissue_id}-value'].str.contains(tissue)]
    alldata = alldata[alldata[f'{parameter_timepoint_id}-value'] == timepoint]
    return alldata

def filter_BioD_with_keyword(filter_keyword, bioD_dataset, parameter_keyword, readout_definitions):
    parameter_keyword_id = get_readout_name_id(parameter_keyword,readout_definitions)
    alldata = list(filter(lambda item: item['readouts'][parameter_keyword_id]['value'] == filter_keyword, bioD_dataset))
    return alldata

def get_parameters_fromBioD(readout_name_id,bioD_data):
    all_parameters = list()
    for element in bioD_data:
        try:
            all_parameters.append(element['readouts'][readout_name_id]['value'])
        except:
            all_parameters.append('not defined')

    all_parameters = list(dict.fromkeys(all_parameters))
    return all_parameters

def imaging_protocol_request(token):
    base_url = f"https://app.collaborativedrug.com/api/v1/vaults/{vault_id}/"
    headers = {'X-CDD-token':f'{token}'}
    BioD_id = 84205
    url = base_url + f"protocols/{BioD_id}"
    bioD_request = requests.request("GET", url, headers=headers).json()
    return bioD_request

@st.cache_data
def BioD_protocol_request(token):
    base_url = f"https://app.collaborativedrug.com/api/v1/vaults/{vault_id}/"
    headers = {'X-CDD-token':f'{token}'}
    BioD_id = 57565
    url = base_url + f"protocols/{BioD_id}"
    bioD_request = requests.request("GET", url, headers=headers).json()
    return bioD_request

def get_BioD_output(readout_name,readout_definitions,data):
    readout_id = get_readout_name_id(readout_name,readout_definitions)
    st.write(f'{readout_name}: {readout_id}')
    data_filtered = list()
    for element in data:
        try:
            data_filtered.append(element['readouts'][readout_id]['value'])
        except:
            pass
    data_filtered = list(data[f'{readout_id}-value'])
    # st.write(data)
    # st.write(data_filtered)
    return data_filtered

@st.cache_data
def filter_BioD_conditions(bioD_data,selected_conditions_all,readout_definitions):
    bioD_data_filtered = bioD_data
    for element in selected_conditions_all.keys():
        condition_list = []
        for condition in selected_conditions_all[element]:
            if selected_conditions_all[element][condition]:
                condition_list.append(condition)
        condition_id = get_readout_name_id(element,readout_definitions)
        bioD_data_filtered = bioD_data_filtered[bioD_data_filtered[f'{condition_id}-value'].isin(condition_list)]  #Filter only elements that are listed in condition_list for each condition
    return bioD_data_filtered

def get_condition_raw_aver(tissues,timepoints,bioD_data_filtered,readout_definitions): 
    keywordlist = []

    readouts_expander = st.expander('Select readouts for reporting:')
    readouts_expander.write('%ID/g or %ID/cc or kBq/cc are required for dosimetry (and Sample mass/body weight for BioD from cut and count, VOI volume/body weight/Inj. act. for BioD from Imaging)')
    for u_nr,uploaded_condition in enumerate(list(bioD_data_filtered.columns)):
        if '-' in uploaded_condition:
            pass
        elif 'Avg' in uploaded_condition:
            pass
        else:
            condition_uploaded_id = uploaded_condition
            readout_definitions_df = pd.DataFrame.from_dict(readout_definitions)
            readout_name = list(readout_definitions_df[readout_definitions_df['id'] == int(condition_uploaded_id)]['name'])[0]
            if readout_name in ['Inj dose per gram',"%ID/g",'tumor volume','body weight','Sample mass',"%ID/cc","kBq/cc",'VOI volume','Tissue','Time point','Subject','Injected activity']:
                preselect_condition = True
            else: 
                preselect_condition = False
            condition_selected = readouts_expander.checkbox(f'{u_nr}: {readout_name}',key=f'{u_nr}-conditionselect', value=preselect_condition)
            if condition_selected:
                keywordlist.append(readout_name)

    if "%ID/g" in keywordlist: keyword_condition_to_get = "%ID/g"  
    elif "%ID/cc" in keywordlist: keyword_condition_to_get = "%ID/cc"  
    elif "kBq/cc" in keywordlist: keyword_condition_to_get = "kBq/cc"  

    readouts_expander.write('****Readouts uploaded:****')
    results_dict = {}
    readouts_cols = readouts_expander.columns(len(keywordlist))
    readouts_cols_df = readouts_expander.columns(len(keywordlist))
    readouts_cols_df_n = readouts_expander.columns(len(keywordlist))
    for k_nr, keyword_selected in enumerate(keywordlist):
        keyword_id = str(get_readout_name_id(keyword_selected, readout_definitions=readout_definitions))
        selected_column = bioD_data_filtered[f'{keyword_id}-value']
        results_dict[keyword_selected] = selected_column
        readouts_cols[k_nr].write(f'{keyword_selected}: ')
        readouts_cols_df[k_nr].write(selected_column.unique())

        readouts_cols_df_n[k_nr].write(f'unique entries:')
        readouts_cols_df_n[k_nr].write(f'N = {len(selected_column.unique())}')

        readouts_cols_df_n[k_nr].write(f'all entries:')
        readouts_cols_df_n[k_nr].write(f'N = {len(selected_column)}')

    results_filtered = pd.DataFrame.from_dict(results_dict, orient='index').T.reset_index()

    results_for_dosimetry_Calc = dict()
    results_for_dosimetry_Calc['Time']=timepoints

    parameter_tissue_id = str(get_readout_name_id('Tissue', readout_definitions=readout_definitions))
    parameter_timepoint_id = str(get_readout_name_id('Time point', readout_definitions=readout_definitions))

    for tissue in tissues:
        inj_dose_averages=[]
        for i,timepoint in enumerate(timepoints):
            # data_per_tissue = bioD_data_filtered[bioD_data_filtered[f'{parameter_tissue_id}-value'].str.contains(str(tissue))]
            data_per_tissue = bioD_data_filtered[bioD_data_filtered[f'{parameter_tissue_id}-value'] == (str(tissue))]
            data_per_tissue = data_per_tissue[data_per_tissue[f'{parameter_timepoint_id}-value'] == timepoint]
            readout_id = get_readout_name_id(keyword_condition_to_get,readout_definitions)
            if f'{readout_id}-value' in list(data_per_tissue.columns):
                results_keyword_pertimepoint = list(data_per_tissue[f'{readout_id}-value'])
            else:
                results_keyword_pertimepoint = list(data_per_tissue[f'{readout_id}'])
            inj_dose_averages.append(np.average(results_keyword_pertimepoint))

        # dosimetry dataframe contains only injected doses (averaged)
        results_for_dosimetry_Calc[tissue]=inj_dose_averages

    results_for_dosimetry_Calc_df = pd.DataFrame.from_dict(results_for_dosimetry_Calc, orient='index')
    results_for_dosimetry_Calc_df.sort_values(by=[time_keyword],axis=1,inplace=True)

    results_for_dosimetry_Calc_df = results_for_dosimetry_Calc_df.T
    # results_for_dosimetry_Calc_df = results_for_dosimetry_Calc_df.replace(to_replace='None', value=np.nan).dropna(how='all')
    results_for_dosimetry_Calc_df = results_for_dosimetry_Calc_df.dropna(axis=1, how='all')
    # results_for_dosimetry_Calc_df = results_for_dosimetry_Calc_df.T

    # remove all timepoints that don't have data associated with: 
    # 1. drop Time axis and transform matrix
    results_for_dosimetry_Calc_cleaned_df = results_for_dosimetry_Calc_df.drop(['Time'],axis=1).T
    # 2. drop all columns where all entries are 'NA' and transform matrix back to previous shape
    results_for_dosimetry_Calc_cleaned_df = results_for_dosimetry_Calc_cleaned_df.dropna(axis=1, how='all').T
    # 3. keep only indices from original frame where tissue data were entered
    results_for_dosimetry_Calc_df = results_for_dosimetry_Calc_df.loc[results_for_dosimetry_Calc_cleaned_df.index].reset_index()
    results_for_dosimetry_Calc_df = results_for_dosimetry_Calc_df.drop(['index'],axis=1)

    # format tissues list: remove Time column and make kidney and blood first organs listed
    tissues = list(results_for_dosimetry_Calc_df.keys())
    tissues.remove('Time')
    # to make kidney the first keyword, blood the second keyword
    if 'kidney' in tissues:
        order = [tissues.index('kidney')]
    else:
        order =[]
    if 'blood' in tissues:
        order.append(tissues.index('blood'))
    for element_keywords in tissues:
        if element_keywords == 'kidney':
            pass
        elif element_keywords != 'blood':
            order.append(tissues.index(element_keywords))
    tissues = [tissues[i_order] for i_order in order]

    return [tissues,results_for_dosimetry_Calc_df,results_filtered,keyword_condition_to_get]

def get_masses_CDD(results_filtered,cdd_weight_key = 'Sample mass',cdd_wb_m_key = 'body weight', cdd_tissue_key= 'Tissue', cdd_tumor_key = 'tumor volume'):
    tissues = set(results_filtered['Tissue'])
    averaged_weight_all_timepoints = list()
    for tissue in tissues:
        results_Averaged_df=results_filtered[results_filtered[cdd_tissue_key] == tissue]
        if 'Subject' in list(results_Averaged_df.columns):
            unique_list = (list(results_Averaged_df['Subject'].drop_duplicates().index))
            results_Averaged_df = results_Averaged_df[['Subject',cdd_weight_key]].loc[results_Averaged_df.index.intersection(unique_list)]

        averaged_weight_all_timepoints.append(   {'tissue': tissue,
                                                  'organ weight': results_Averaged_df[cdd_weight_key].mean(),
                                                  'organ weight (n)': len(results_Averaged_df[cdd_weight_key])})

    # Get unique subjects only (to avoid counting body weights in multiple times)
    if 'Subject' in list(results_Averaged_df.columns):
        unique_list = (list(results_filtered['Subject'].drop_duplicates().index))
        results_Averaged_df = results_filtered.loc[results_filtered.index.intersection(unique_list)]

    if cdd_wb_m_key in list(results_Averaged_df.columns):
        averaged_weight_all_timepoints.append(   {'tissue': 'total body',
                                                    'organ weight': results_Averaged_df[cdd_wb_m_key].mean(),
                                                    'organ weight (n)': len(results_Averaged_df[cdd_wb_m_key])})
    if cdd_tumor_key in list(results_Averaged_df.columns):
        averaged_weight_all_timepoints.append(   {'tissue': 'tumor (volume)',
                                                    'volume': results_Averaged_df[cdd_tumor_key].mean(),
                                                    'volume (n)': len(results_Averaged_df[cdd_tumor_key])})
    
    averaged_weight_all_timepoints_df = pd.DataFrame.from_records(averaged_weight_all_timepoints)

    with st.expander(f'organ weights'):
        st.write(averaged_weight_all_timepoints_df)
    return averaged_weight_all_timepoints_df

def select_decay_correction_input():
    decay_correction = st.checkbox('Perform decay correction of input data? - Check if yes',value=True,on_change=decay_correction_changed)
    with st.expander('Instructions for decay correction'):
        st.write('''If your data input is biological decay only (decay corrected to time of injection) you should add the isotope specific decay to get the real number of disintegrations. (Check box)
        If your data input is not decay corrected (decay corrected only back to the time the animal was taken down) your data will represent both physical and biological decay. Do not perform decay correction again (Uncheck box).
        ''')
    return decay_correction

def select_radioisotope1(rayz_id, text = 'Radioisotope 1: used in BioD study',preselect=''):
    try: #try to get isotope from RAYZ ID
        isotope_BioD = rayz_id.split('-')[-1][-2:]+'-'+rayz_id.split('-')[-1][:-2]
        isotope_BioD_nr = list(isotopes_BioD_halflives.keys()).index(isotope_BioD)
    except:
        isotope_BioD_nr = 0

    if len(preselect) > 0:
        isotope_BioD_nr = list(isotopes_BioD_halflives.keys()).index(preselect)

    radioisotope1 = st.radio(text+f' for {rayz_id}', options=(isotopes_BioD_halflives.keys()), index=isotope_BioD_nr, key=text)                                
    st.write(f'You selected {radioisotope1}')
    return radioisotope1

def admAct_corr_kBqcc(results_for_dosimetry_Calc_df,administered_activity):
    rawdata = results_for_dosimetry_Calc_df.copy(deep=False)
    for datakey in results_for_dosimetry_Calc_df.keys():
        if datakey != time_keyword:
            results_for_dosimetry_Calc_df[datakey] = 100* results_for_dosimetry_Calc_df[datakey] / (administered_activity)
    else:
        results_for_dosimetry_Calc_df = results_for_dosimetry_Calc_df.copy(deep=False)
    results_for_dosimetry_Calc_df = results_for_dosimetry_Calc_df.sort_index()
    return [rawdata,results_for_dosimetry_Calc_df]

def decay_corr_kBqcc(results_for_dosimetry_Calc_df,radioisotope1):
    rawdata = results_for_dosimetry_Calc_df.copy(deep=False)
    for datakey in results_for_dosimetry_Calc_df.keys():
        if datakey != time_keyword:
            # results_for_dosimetry_Calc_df[datakey+' [kBqcc]'] = results_for_dosimetry_Calc_df[datakey]
            results_for_dosimetry_Calc_df[datakey] = results_for_dosimetry_Calc_df[datakey] / np.power(0.5, results_for_dosimetry_Calc_df[time_keyword]/isotopes_BioD_halflives[radioisotope1])
    else:
        results_for_dosimetry_Calc_df = results_for_dosimetry_Calc_df.copy(deep=False)
    results_for_dosimetry_Calc_df = results_for_dosimetry_Calc_df.sort_index()
    return [rawdata,results_for_dosimetry_Calc_df]

def decay_corr_input(results_for_dosimetry_Calc_df,radioisotope1):
    rawdata = results_for_dosimetry_Calc_df.copy(deep=False)
    if st.session_state.decay_correction:
        for datakey in results_for_dosimetry_Calc_df.keys():
            if datakey != time_keyword:
                # results_for_dosimetry_Calc_df[datakey+' not corrected'] = results_for_dosimetry_Calc_df[datakey]
                results_for_dosimetry_Calc_df[datakey] = results_for_dosimetry_Calc_df[datakey] * np.power(0.5, results_for_dosimetry_Calc_df[time_keyword]/isotopes_BioD_halflives[radioisotope1])
    else:
        results_for_dosimetry_Calc_df = results_for_dosimetry_Calc_df.copy(deep=False)
    results_for_dosimetry_Calc_df = results_for_dosimetry_Calc_df.sort_index()
    return [rawdata,results_for_dosimetry_Calc_df]

def download_cdd_rawdata(rayz_id,radioisotope1,results_for_dosimetry_Calc_df,rawdata,averaged_weight_all_timepoints_df,tissues,results_filtered,rawdata_kBqcc=[],rawdata_IDcc=[]):
    output_data = io.BytesIO()
    with pd.ExcelWriter(output_data, engine='xlsxwriter') as writer:  
        if st.session_state.decay_correction:
            sheet_title = f"{radioisotope1}-DecayCorrected"
            results_for_dosimetry_Calc_df.to_excel(writer, sheet_name = sheet_title, index = True)
        if len(rawdata_kBqcc) > 0 and len(rawdata_IDcc):
            rawdata.to_excel(writer, sheet_name = 'Input IDcc Biol Decay Only', index = True)
            rawdata_IDcc.to_excel(writer, sheet_name = 'Input kBqcc Biol Decay Only', index = True)
            rawdata_kBqcc.to_excel(writer, sheet_name = 'CDD input kBqcc', index = True)
        else:
            rawdata.to_excel(writer, sheet_name = 'CDD input', index = True)

        averaged_weight_all_timepoints_df.to_excel(writer, sheet_name = 'aver weights from CDD', index = True)
        results_filtered.to_excel(writer, sheet_name = 'All data', index = True)
        if 'Tissue' in results_filtered.columns:
            for tissue in tissues:
                results_df = results_filtered[results_filtered['Tissue']==tissue]
                results_df.to_excel(writer, sheet_name = tissue+'-input', index = True)


    st.download_button(label='Download CDD rawdata',
                            data=output_data.getvalue(),
                            file_name= f'{rayz_id}-CDD_export.xlsx')

def plot_input_data_biod(results_for_dosimetry_Calc_df,condition,rayz_id):
    timepoints = results_for_dosimetry_Calc_df['Time']
    results_for_dosimetry_Calc_df = results_for_dosimetry_Calc_df.T

    results_for_plotting = results_for_dosimetry_Calc_df
    results_for_plotting.columns = results_for_plotting.iloc[0]
    results_for_plotting.drop(results_for_plotting.index[0])

    traces=[]
    colorrange = int(255./(len(timepoints)+2)) if len(timepoints) < 10 else 10
    
    for it,timepoint in enumerate(timepoints):
        timepoint = float(timepoint)
        tissue_names = np.asarray(results_for_plotting.index[1:])
        inj_dose_averages = np.asarray(results_for_plotting[timepoint][1:])

        # create RGB color range to plot timepoints in different colors
        color_r = 0+it*colorrange if it*colorrange < 255 else 255
        color_r = color_r if color_r >= 0 and color_r <= 255 else 10
        color_g = 50+it*colorrange if 50+it*colorrange < 255 else 255
        color_g = color_g if color_g >= 0 and color_g <= 255 else 10
        color_b = 50+it*colorrange if 50+it*colorrange < 255 else 50
        color_b = color_b if color_b >= 0 and color_b <= 255 else 10


        trace_barplot = go.Bar(x = tissue_names, y = inj_dose_averages,
                        marker = dict(color='rgb('+str(color_r)+','+str(color_g)+','+str(color_b)+')'),
                        name = f'{timepoint}h',
                        offset=(it - 2) * 0.1,
                        width=0.1,                                               
                        )
        traces.append(trace_barplot)   

    fig = go.Figure()
    fig.add_traces(traces)
    fig.update_layout(
            title=f'{condition} for {rayz_id} from BioD',
            yaxis_title=f'{condition}',
        )
    st.plotly_chart(fig, use_container_width=True)


@st.cache_resource
def convert_df_to_csv(df):
   return df.to_csv().encode('utf-8')

def monoexpdecay(x, A1, lambda1):
    warnings.filterwarnings('ignore')
    return A1 * np.exp(-lambda1 * x)

def biexpdecay(x, A1, lambda1, A2, lambda2):
    return A1 * np.exp(-lambda1 * x) + A2 * np.exp(-lambda2 * x)

def biexp_abs_elim(x, ka, ke, cl):
    return ((ke*ka)/(cl*(ka-ke)))*(np.exp(-1*ke*x)-np.exp(-1*ka*x))

def monofit_func(x,y,x_fit):
    monoexp_fit = np.polyfit(x,np.log(y),1)
    A1      = np.exp(monoexp_fit[1])
    lambda_fit  = - monoexp_fit[0]
    thalf = np.log(2)/lambda_fit
    y_fit = A1 * np.exp(- lambda_fit * x_fit)
    AUC_calculated = (A1/lambda_fit)    # from book M.Stabin: Basic Principles of Internal Dosimetry Calculations, page 7
    # AUC_numerical, AUC_err = quad(monoexpdecay,0,np.inf,args=(A1,lambda_fit))
    last_timepoint = float(list(x)[-1])
    AUC_numerical_extrapolated, AUC_err = quad(monoexpdecay,last_timepoint,np.inf,args=(A1,lambda_fit))
    return [A1, lambda_fit, thalf, y_fit, AUC_calculated,AUC_numerical_extrapolated]

def monofit_dec_func(x,y,x_fit) :
    popt, pcov = curve_fit(monoexpdecay, x, y, p0 = (1.,1.))
    A1,lambda_fit      = popt
    thalf = np.log(2)/lambda_fit
    AUC_calculated = (A1/lambda_fit)    # from book M.Stabin: Basic Principles of Internal Dosimetry Calculations, page 7
    AUC_numerical, AUC_err = quad(monoexpdecay,0,np.inf,args=(A1,lambda_fit))
    last_timepoint = float(list(x)[-1])
    AUC_numerical_extrapolated, AUC_err = quad(monoexpdecay,last_timepoint,np.inf,args=(A1,lambda_fit))
    y_fit = A1 * np.exp(- lambda_fit * x_fit) 
    return [A1, lambda_fit, thalf, y_fit, AUC_calculated,AUC_numerical_extrapolated]

def biexp_func(x,y,x_fit):
    popt, pcov = curve_fit(biexpdecay, x, y, p0 = (1.,1.,1.,1.), bounds=(0,np.inf))
    A1,lambda1,A2,lambda2 = popt
    last_timepoint = float(list(x)[-1])
    if lambda1 < 1e-4:
        AUC_calculated = (A2/lambda2)
        AUC_numerical, AUC_err = quad(monoexpdecay,0,np.inf,args=(A2,lambda2))
        AUC_numerical_extrapolated, AUC_err = quad(monoexpdecay,last_timepoint,np.inf,args=(A2,lambda2))
    elif lambda2 < 1e-4:
        AUC_calculated = (A1/lambda1)
        AUC_numerical, AUC_err = quad(monoexpdecay,0,np.inf,args=(A1,lambda1))
        AUC_numerical_extrapolated, AUC_err = quad(monoexpdecay,last_timepoint,np.inf,args=(A1,lambda1))
    else:
        AUC_calculated = (A1/lambda1) + (A2/lambda2)
        AUC_numerical, AUC_err = quad(biexpdecay,0,np.inf,args=(A1,lambda1, A2, lambda2))
        AUC_numerical_extrapolated, AUC_err = quad(biexpdecay,last_timepoint,np.inf,args=(A1,lambda1, A2, lambda2))
    # st.write(AUC_calculated)
    # st.write(AUC_numerical)

    y_fit = biexpdecay(x_fit,A1, lambda1, A2, lambda2)
    return [A1,lambda1,A2,lambda2,y_fit,AUC_calculated,AUC_numerical_extrapolated]

def trapezoidal_func(x,y):
    xtrap = pd.DataFrame(np.insert(x.values, 0, values=[0.], axis = 0))[0]

    try:
        ytrap = pd.DataFrame(np.insert(y.values, 0, values=[list(y)[0]], axis = 0))[0]
    except:
        ytrap = pd.DataFrame(np.insert(y.values, 0, values=[0.], axis = 0))[0]
        # st.error('no rawdata uploaded for this organ, please uncheck from tissue list')

    AUC_calculated = 0
    for ie,element_ie in enumerate(xtrap):
        if ie == 0:
            pass
        else:
            y2 = ytrap[ie]
            y1 = ytrap[ie-1]

            x2 = xtrap[ie]
            x1 = xtrap[ie-1]

            AUC_calculated += ((y2+y1)/2)*(x2-x1)  
    return [xtrap, ytrap, AUC_calculated]

def trap_linextrapol_func(x_raw,y_raw, halflive, write_error_warning = True):
    try:
        b, a, r, p, std = linregress(x_raw[-2:],y_raw[-2:]) # y=a+b*x
    except:
        try:
            b, a = [-2,list(y_raw)[-1]]
        except:
            st.error('no rawdata uploaded for this organ, please uncheck from tissue list')

    if b > 0:
        if write_error_warning:
            st.error('slope is not negative, physical decay to 2% of last y value was used for interpolation')
        decaytime = np.log(0.02)/(np.log(0.5)/halflive)  #decay to 2% of last y value
        x_end = list(x_raw)[-1] + decaytime
    else:
        x_end = -a/b

    xtrap = pd.DataFrame(np.insert(x_raw.values, 0, values=[0.], axis = 0))
    xtrap = pd.DataFrame(np.insert(xtrap.values, len(xtrap.index), values=[x_end], axis = 0))[0]
    ytrap = pd.DataFrame(np.insert(y_raw.values, 0, values=[list(y_raw)[0]], axis = 0))
    ytrap = pd.DataFrame(np.insert(ytrap.values, len(ytrap.index), values=[0.], axis = 0))[0]

    AUC_calculated = 0
    for ie,element_ie in enumerate(xtrap):
        if ie == 0:
            pass
        else:
            y2 = ytrap[ie]
            y1 = ytrap[ie-1]
            x2 = xtrap[ie]
            x1 = xtrap[ie-1]
            AUC_calculated += ((y2+y1)/2)*(x2-x1) 
            # st.write(((y2+y1)/2)*(x2-x1))
    AUC_extrapolated = ((list(ytrap)[-1]+list(ytrap)[-2])/2)*(list(xtrap)[-1]-list(xtrap)[-2]) 
    return [xtrap, ytrap, AUC_calculated,AUC_extrapolated]


def trap_physdecayextrapol_func(x_raw,y_raw, halflive):
    try:
        b, a, r, p, std = linregress(x_raw[-2:],y_raw[-2:]) # y=a+b*x
    except:
        try:
            b, a = [-2,list(y_raw)[-1]]
        except:
            st.error('no rawdata uploaded for this organ, please uncheck from tissue list')

    x_end = list(x_raw)[-1]

    
    xtrap = pd.DataFrame(np.insert(x_raw.values, 0, values=[0.], axis = 0))
    xtrap = pd.DataFrame(np.insert(xtrap.values, len(xtrap.index), values=[x_end], axis = 0))[0]
    ytrap = pd.DataFrame(np.insert(y_raw.values, 0, values=[list(y_raw)[0]], axis = 0))
    ytrap = pd.DataFrame(np.insert(ytrap.values, len(ytrap.index), values=[0.], axis = 0))[0]

    AUC_calculated = 0
    for ie,element_ie in enumerate(xtrap):
        if ie == 0:
            pass
        else:
            y2 = ytrap[ie]
            y1 = ytrap[ie-1]
            x2 = xtrap[ie]
            x1 = xtrap[ie-1]
            AUC_calculated += ((y2+y1)/2)*(x2-x1) 
    AUC_meas = AUC_calculated
    AUC_extrapolated = list(y_raw)[-1] / (np.log(2)/halflive)
    AUC_calculated += AUC_extrapolated
    # AUC_extrapolated_ratio = 100* AUC_extrapolated/AUC_calculated
    # st.write(f'measured AUC: {AUC_meas:.2f}; extrapolated AUC: {AUC_extrapolated:.2f}')
    # st.write(f'={AUC_extrapolated_ratio:.2f}% of AUC is extrapolated)')
    return [xtrap, ytrap, AUC_calculated,AUC_extrapolated]

def biexpelim_func(x_raw,y_raw,x_fit): # from David Huang, reference needed
    popt, pcov = curve_fit(biexp_abs_elim, x_raw, y_raw, p0=(0.01, 0.1, 0.001))
    ka,ke,cl = popt
    y_fit = biexp_abs_elim(x_fit,ka,ke,cl)
    AUC_calculated = (1/ke-1/ka)*(ke*ka)/(cl*(ka-ke))
    last_timepoint = float(list(x_raw)[-1])
    AUC_numerical_extrapolated, AUC_err = quad(biexp_abs_elim,last_timepoint,np.inf,args=(ka,ke,cl))
    return [ka, ke, cl, y_fit, AUC_calculated,AUC_numerical_extrapolated]

def biexp_abs_elim2(x, lambda1, lambda2, A1):
    return (A1*(np.exp(-lambda1*x)-np.exp(-lambda2*x)))

def biexp_abs_elim3(x, lambda1, lambda2, A1, A2):
    return (A1*(np.exp(-lambda1*x))-(A2*np.exp(-lambda2*x)))

def biexpelim_func2(x_raw,y_raw,x_fit): # from Voximetry
    popt, pcov = curve_fit(biexp_abs_elim2, x_raw, y_raw, p0=(0.01, 0.1, 0.001))
    lambda1, lambda2, A1 = popt
    y_fit = biexp_abs_elim2(x_fit,lambda1, lambda2, A1)
    AUC_calculated = A1*((1/lambda1) - (1/lambda2))
    last_timepoint = float(list(x_raw)[-1])
    AUC_numerical_extrapolated, AUC_err = quad(biexp_abs_elim2,last_timepoint,np.inf,args=(lambda1, lambda2, A1))
    return [lambda1, lambda2, A1, y_fit, AUC_calculated,AUC_numerical_extrapolated]

def biexpelim_func3(x_raw,y_raw,x_fit): # from Voximetry
    popt, pcov = curve_fit(biexp_abs_elim3, x_raw, y_raw, p0=(0.01, 0.1, 1, 0.001))
    lambda1, lambda2, A1, A2 = popt
    y_fit = biexp_abs_elim3(x_fit,lambda1, lambda2, A1, A2)
    AUC_calculated = (A1*(1/lambda1)) - (A2* (1/lambda2))
    last_timepoint = float(list(x_raw)[-1])
    AUC_numerical_extrapolated, AUC_err = quad(biexp_abs_elim3,last_timepoint,np.inf,args=(lambda1, lambda2, A1, A2))
    return [lambda1, lambda2, A1,A2, y_fit, AUC_calculated,AUC_numerical_extrapolated]


def definition_scaling_method(scaling_method, fitmodel):
    if scaling_method == no_scaling:
        st.write(f'Method #1 in https://doi.org/10.1155/2019/6438196')
        st.latex(r'''  
            TIAC_{human,organ} = TIAC_{mouse,organ}
                ''')
        st.latex(r'''  TIAC = {time-integrated} \space {activity} \space {coefficient}''')
    elif scaling_method == rel_mass_scal:
        st.write(f'Method #2 in https://doi.org/10.1155/2019/6438196')
        st.latex(r'''  
            TIAC_{human,organ} = TIAC_{mouse,organ} * \left(\frac{\left(\frac{m_{organ}}{m_{WB}}\right)_{human}}{\left(\frac{m_{organ}}{m_{WB}}\right)_{mouse}}\right)
                ''')            
        st.latex(r'''  WB = {whole} \space {body}''')
        st.latex(r'''  TIAC = {time-integrated} \space {activity} \space {coefficient}''')

    # elif scaling_method == rel_mass_scal_fda:
        st.write(f'Also defined by FDA: https://www.fda.gov/media/129547/download (Oncology Therapeutic Radiopharmaceuticals - Guideline)')
        # st.latex(r'''  
        #     \tau_{human,organ} = \tau_{mouse,organ} * \left(\frac{\left(\frac{m_{organ}}{m_{WB}}\right)_{human}}{\left(\frac{m_{organ}}{m_{WB}}\right)_{mouse}}\right)
        #         ''')   
        st.latex(r'''  
            \%ID_{human,organ} = \%ID_{mouse,organ} * \left(\frac{\left(\frac{m_{organ}}{m_{WB}}\right)_{human}}{\left(\frac{m_{organ}}{m_{WB}}\right)_{mouse}}\right)
                ''')   

        # st.latex(r'''  \tau = {residence} \space {time} = TIAC''')
        st.latex(r'''  \%ID = {injected} \space {dose} ''')
        if st.checkbox('Show FDA guideline'):
            st.image(Image.open('./pics/FDA-guideline-radiopharmaceuticals.png'), caption='Image 1: Oncology Therapeutic Radiopharmaceuticals - Guideline')  
    elif scaling_method == time_scal:
        st.write(f'You selected: {scaling_method} with fit model {fitmodel}')
        st.write(f'Method #3 in https://doi.org/10.1155/2019/6438196')
        st.latex(r'''  
            t_{human,organ} = t_{mouse,organ} *  {\left(\frac{m_{WB,human}}{m_{WB,mouse}}\right)}^{0.25}
                ''')            
        st.latex(r'''  
            TIAC_{human,organ} = TIAC_{mouse,organ}
                ''')
        st.latex(r'''  WB = {whole} \space {body}''')
        st.latex(r'''  TIAC = {time-integrated} \space {activity} \space {coefficient}''')
    elif scaling_method == 'combined time and relative mass scaling - FDA guideline':
        st.write(f'You selected: {scaling_method} with fit model {fitmodel}')

        st.write(f'Also used by FDA: https://www.fda.gov/media/129547/download (Oncology Therapeutic Radiopharmaceuticals - Guideline)')
        st.latex(r'''  
            \tau_{human,organ} (= \frac{1}{\lambda}) = \tau_{mouse,organ} * \left(\frac{\left(\frac{m_{organ}}{m_{WB}}\right)_{human}}{\left(\frac{m_{organ}}{m_{WB}}\right)_{mouse}}\right)
                ''')            
        st.latex(r'''  
            TIAC_{human,organ} = TIAC_{mouse,organ} * \left(\frac{\left(\frac{m_{organ}}{m_{WB}}\right)_{human}}{\left(\frac{m_{organ}}{m_{WB}}\right)_{mouse}}\right)
                ''')  
        st.latex(r'''  WB = {whole} \space {body}''')
        st.latex(r'''  TIAC = {time-integrated} \space {activity} \space {coefficient}''')
        if st.checkbox('Show FDA guideline'):
            st.image(Image.open('./pics/FDA-guideline-radiopharmaceuticals.png'), caption='Image 1: Oncology Therapeutic Radiopharmaceuticals - Guideline')  
    elif scaling_method == 'combined time and relative mass scaling':
        st.write(f'You selected: {scaling_method} with fit model {fitmodel}')

        st.write(f'Method #4 in https://doi.org/10.1155/2019/6438196')
        st.latex(r'''  
            t_{human,organ} = t_{mouse,organ} *  {\left(\frac{m_{WB,human}}{m_{WB,mouse}}\right)}^{0.25}
                ''')            
        st.latex(r'''  
            TIAC_{human,organ} = TIAC_{mouse,organ} * \left(\frac{\left(\frac{m_{organ}}{m_{WB}}\right)_{human}}{\left(\frac{m_{organ}}{m_{WB}}\right)_{mouse}}\right)
                ''')  
        st.latex(r'''  WB = {whole} \space {body}''')
        st.latex(r'''  TIAC = {time-integrated} \space {activity} \space {coefficient}''')
    elif scaling_method == 'allometric scaling':
        st.write(f'Method #5 in https://doi.org/10.1155/2019/6438196 (BUT misreferenced in this publication - has to be revisited)')
    elif scaling_method == alpha_scal:
        st.write(f'Method by Stephen Graves (consultant 01/2022) ')
        st.write('1. Decay correction of %ID/g (corrected for decaying nucleotide during experimental time course)') # from report 'RayzeBio-SG_PSMA617_BioD_Report_1.0.pdf' page 2
        st.latex(r'''  
            A(t) = \%ID/g_{mouse, decay \space corrected}(t) = \%ID/g_{mouse}(t) * {e^{-\lambda_{radioisotope1} * t}} 
                ''')  
        st.latex(r'''  radioisotope  \space 1 : \space used \space in \space BioD \space study ''')
        st.write('2. Calculate AUC, TIAC(mouse)') # from slides '20220113_RayzeBio-Dosimetry_Graves_v1.0.pptx' slide 4
        st.latex(r'''  
                TIAC_{mouse,organ} = \int_{0}^{\infty} A_{organ} *  {e^{-\lambda_{organ} * t}} dt = \left(\frac{A_{organ}}{\lambda_{organ}}\right)
                ''')            
        st.write('3. Extrapolate to TIAC(human)') # from slides '20220113_RayzeBio-Dosimetry_Graves_v1.0.pptx' slide 4
        st.latex(r'''  
            TIAC_{human,organ} = \alpha * TIAC_{mouse,organ} * \left(\frac{\left(\frac{m_{organ}}{m_{WB}}\right)_{human}}{\left(\frac{m_{organ}}{m_{WB}}\right)_{mouse}}\right)
                ''')

        st.write('For clearance organs (kidney):') # from slides '20220113_RayzeBio-Dosimetry_Graves_v1.0.pptx' slide 7
        st.latex(r'''  
            TIAC_{human,organ} = \alpha * TIAC_{mouse,organ} 
                ''')
            # TIAC_{human,organ} = \frac{\alpha}{\lambda_{mouse, organ}} * TIAC_{mouse,organ} 
        
    elif scaling_method == alpha_scal:
        st.write(f'Method by Stephen Graves (consultant 01/2022) ')          
        st.write('Extrapolation to TIAC(human)') # from slides '20220113_RayzeBio-Dosimetry_Graves_v1.0.pptx' slide 4
        st.latex(r'''  
            TIAC_{human,organ} = \alpha * TIAC_{mouse,organ} * \left(\frac{\left(\frac{m_{organ}}{m_{WB}}\right)_{human}}{\left(\frac{m_{organ}}{m_{WB}}\right)_{mouse}}\right)
                ''')
        st.latex(r'''  
            \alpha = {\left(\frac{m_{WB,human}}{m_{WB,mouse}}\right)}^{0.25}
                ''')
        st.write('For clearance organs (kidney, liver):') # from slides '20220113_RayzeBio-Dosimetry_Graves_v1.0.pptx' slide 7
        st.latex(r'''  
            TIAC_{human,organ} = \alpha * TIAC_{mouse,organ} 
                ''')

        st.latex(r'''  WB = {whole} \space {body}''')
        st.latex(r'''  m = mass''')
        st.latex(r'''  TIAC = {time-integrated} \space {activity} \space {coefficient}''')

    else:
        st.write('Not defined yet')

def get_mass_from_file(tissue, tissue_masses_file_df, keyword_tissues, keyword_masses, if_nan_option):
    tissue_spellings = get_all_possible_tissue_spellings(tissue)
    try:
        mass = np.nan
        for option in tissue_spellings:
            if option in list(tissue_masses_file_df[keyword_tissues]):
                mass = tissue_masses_file_df[tissue_masses_file_df[keyword_tissues] == option].iloc[0][keyword_masses]
                break
    except:
        st.error(f'{tissue} mass cannot be assigned automatically for {keyword_tissues}, please choose from options below')
        tissue_new = st.selectbox(f'Which tissue mass do you want to use for {tissue} ({keyword_tissues})',tissue_masses_file_df[keyword_tissues],key=f'{tissue}-selectmass-{keyword_tissues}')
        mass = tissue_masses_file_df[tissue_masses_file_df[keyword_tissues] == tissue_new].iloc[0][keyword_masses]
    if pd.isna(mass):
        mass = if_nan_option
    return mass

def get_all_possible_tissue_spellings(tissue_input):
    tissue_input = tissue_input.lower()
    tissue_input = tissue_input.strip()
    tissue_option = tissue_input.split(' ')
    results=[]
    if tissue_input not in tissue_option:
        tissue_option.insert(0,tissue_input)
    for tissue in tissue_option:
        poss_singular = []
        poss_singular.append(tissue)
        if tissue == 'adrenal glands':
            poss_singular.append('Adrenals')
        elif tissue == 'salivary glands':
            poss_singular.append('Salivary')
        elif tissue == 'bladder wall':
            poss_singular.append('UB Cont')
            poss_singular.append('Urinary Bladder Wall')
        elif tissue == 'bone marrow' or tissue == keyword_projected_bm or fnmatch(tissue,'*bone marrow*') or fnmatch(tissue,'*Red Mar*') or fnmatch(tissue,'*Red Mar.*'):
            poss_singular.append('Red Marrow') 
            poss_singular.append('red marrow') 
            poss_singular.append('red mar.') 
            poss_singular.append('Bone Marrow') 
        elif tissue == 'bone' and tissue_input != keyword_projected_bm and not fnmatch(tissue_input,'*bone marrow*'):
            poss_singular.append('CortBone')
            poss_singular.append('TrabBone')        
            poss_singular.append('Osteogenic Cells') 
        elif tissue == 'gallbladder':
            poss_singular.append('GB Cont')
            poss_singular.append('Gallbladder Wall')        
        elif tissue == 'heart':
            poss_singular.append('HeartCon')
            poss_singular.append('Heart Wall')  
        elif tissue == 'stomach' or tissue == 'stomach contents' or tissue == 'stomcont':
            poss_singular.append('Stomach Contents') 
            poss_singular.append('StomCont')
            poss_singular.append('Stomach Wall') 
        elif fnmatch(tissue,'*small intestine*') or fnmatch(tissue,'*small intestine contents*') or fnmatch(tissue,'SI'):
            poss_singular.append('Small Intestine') 
            poss_singular.append('SI Cont')
        elif tissue == 'large intestine' or tissue == 'large intestine contents':
            poss_singular.append('Left Colon') 
            poss_singular.append('LLI Wall')
            poss_singular.append('LLI Cont ')
            poss_singular.append('ULI Cont') 
            poss_singular.append('ULI Wall') 
            # poss_singular.append('Rectum') 
        elif tissue == 'cecum':
            poss_singular.append('Right Colon') 
            poss_singular.append('ULI Cont') 
            poss_singular.append('ULI Wall') 
        elif tissue == 'vena cava':
            poss_singular.append('blood')   
        elif tissue == 'total body':
            poss_singular.append('Tot Body')   
        for element in poss_singular:
            elements = element+'s'
            results.extend([element,elements,element.title(),elements.title(),element.title()+' ',elements.title()+' '])
            results.extend([element[:-2].lower()])
            results.extend([element.lower(),elements.lower()])
    return results

def get_matching_organnames(tissue,csv_organs_list):
    tissue_all_spellings = get_all_possible_tissue_spellings(str(tissue).strip())
    match_res_firstRow = set(csv_organs_list) & set(tissue_all_spellings)
    if len(match_res_firstRow) > 0:
        tissue_found_firstRow_list = str(match_res_firstRow)[1:-1].split(',')
        sourcelist = []
        for fi, foundresult in enumerate(tissue_found_firstRow_list):
            if fi == 0:
                sourcelist.append(foundresult[1:-1])
            else:
                sourcelist.append(foundresult[2:-1])
        source_organ_defined = True
    else:
        tissue_all_spellings = get_all_possible_tissue_spellings('total body')
        match_res_TotBody = set(csv_organs_list) & set(tissue_all_spellings)
        
        sourcelist = []
        for fi, foundresult in enumerate(match_res_TotBody):
                sourcelist.append(foundresult)
        source_organ_defined = False
    return [source_organ_defined,sourcelist]

def select_bonemarrow_projection(fitresults):
    if 'blood' in list(fitresults.columns):
        blood_key = 'blood'
    elif len(set(['L1','L2','L3','L4','L5'])&set(list(fitresults.columns))) > 0:
        blood_key = st.multiselect('Which tissue is a readout for blood?', options=list(fitresults.columns)[2:])
    else:
        blood_key = st.selectbox('Which tissue is a readout for blood?', options=list(fitresults.columns)[2:])
    with st.form('Bone Marrow Projection Options'):
        bm_method = st.radio(
                f'Bone Marrow Projection Options',
                (bm_options))

        st.form_submit_button(f'Confirm Bone Marrow Scaling')

    definition_boneMarrow_methods(bm_method)
    return [bm_method,blood_key]

def project_bonemarrow(fitresults,blood_key,bm_method):
    if type(blood_key) == list:
        tiac_g_blood = 0
        for blood_key_el in blood_key:
            tiac_g_blood += float(fitresults[blood_key_el].T[keyword_mtiac_g])
        tiac_g_blood/=len(blood_key) #Average over all submitted keys
    else:
        tiac_g_blood = float(fitresults[blood_key].T[keyword_mtiac_g])
    bm_params = dict()
    if bm_method == bm_sgourus:
        RMECFF = float(st.text_input('RMECFF',value='0.19'))
        HCT = float(st.text_input('HCT',value='0.47'))
        tiac_g_rm = tiac_g_blood * (RMECFF/(1 - HCT))
        bm_params['RMECFF'] = RMECFF
        bm_params['HCT'] = HCT
        
    elif bm_method == bm_sgourus_greg:
        RMECFF = float(st.text_input('RMECFF',value='0.19'))
        HCT = float(st.text_input('HCT',value='0.47'))
        m_RM_patient = float(st.text_input('Red Marrow mass patient (g)',value='1170'))
        m_BL_patient = float(st.text_input('Blood mass patient (g)',value='6363'))
        tiac_g_rm = tiac_g_blood * (RMECFF/(1 - HCT)) * m_RM_patient/m_BL_patient
        bm_params['RMECFF'] = RMECFF
        bm_params['HCT'] = HCT
        bm_params['m_RM_patient'] = m_RM_patient
        bm_params['m_BL_patient'] = m_BL_patient

    elif bm_method == bm_traino:
        RMBLR = float(st.text_input('RMBLR',value='1'))
        m_RM_patient = float(st.text_input('Red Marrow mass patient (g)',value='1170'))
        tiac_g_rm = tiac_g_blood * RMBLR * m_RM_patient/1000
        bm_params['RMBLR'] = RMBLR
         
    elif bm_method == bm_grudzinski:
        m_RM_patient = float(st.text_input('Red Marrow mass patient (g)',value='1170'))
        m_wb_animal = float(st.text_input('Total mass animal (g)',value='25'))
        m_wb_human = float(st.text_input('Total mass human (g)',value='73000'))
        tiac_g_rm = tiac_g_blood * m_wb_animal/m_wb_human * m_RM_patient
        bm_params['m_RM_patient'] = m_RM_patient
        bm_params['m_wb_animal'] = m_wb_animal
        bm_params['m_wb_human'] = m_wb_human

    elif bm_method == bm_graves:
        RMBLR = float(st.text_input('RMBLR',value='1'))
        tiac_g_rm = tiac_g_blood * RMBLR
        bm_params['RMBLR'] = RMBLR

    else:
        st.error('This method is not defined yet') 

    # # add parameters for bone marrow projection to fitresults dataframe
    bm_fitresults = []
    for param in fitresults['parameter']:
        if param == keyword_mtiac_g:
            bm_fitresults.append(round(tiac_g_rm,4))
        elif param == 'fitmodel':
            bm_fitresults.append(f'Bone Marrow projection: {bm_method} from organ: <{blood_key}>')
        else:
            bm_fitresults.append(None)

    fitresults = fitresults.T
    for bm_k, bm_v in bm_params.items():
        fitresults[bm_k] = None
        bm_fitresults.append(bm_v)
    fitresults = fitresults.T

    fitresults[keyword_projected_bm] = bm_fitresults
    columns_blood = []
    if type(blood_key) == list:
        columns_blood = list(blood_key).copy()
        columns_blood.extend([keyword_projected_bm])
    else:
        columns_blood = [blood_key,keyword_projected_bm]
    st.write(fitresults[columns_blood].T[[keyword_mtiac_g,'fitmodel']])
    fitresults['parameter'] = list(fitresults.index)
    return fitresults


def definition_boneMarrow_methods(bm_method):
    if bm_method == bm_sgourus:
        st.write(f'Method described by Sgourus (1993): https://jnm.snmjournals.org/content/34/4/689.long ')
        st.latex(r'''  
            [A]_{RM} = [A]_{BL} * \frac{RMECFF}{1 - HCT} 
                ''')
        st.latex(r'''[A]_{RM} = \text{Activity concentration in Red Marrow (Residence time in organ = } mTIAC_{org}\text{)}''')
        st.latex(r'''[A]_{BL} = \text{Activity concentration in Blood (Residence time in organ = } mTIAC_{org}\text{)}''')
        st.latex(r''' RMECFF = \text{Red Marrow Extracellular Fluid Fraction (physiological: 0.19; range 0.15-0.25)}''')
        st.latex(r''' HCT = \text{Hematocrit (physiological: 0.47; range 0.2-0.6)}''')    

    elif bm_method == bm_sgourus_greg:
        st.write(f'Method described by Sgourus (1993), adjusted by Greg: https://jnm.snmjournals.org/content/34/4/689.long ')
        st.latex(r'''  
            [A]_{RM} = [A]_{BL} * \frac{RMECFF}{1 - HCT} * \left(\frac{m_RM\left(patient\right)}{m_BL\left(patient\right)}\right)
                ''')
        st.latex(r'''[A]_{RM} = \text{Activity concentration in Red Marrow (Residence time in organ = } mTIAC_{org}\text{)}''')
        st.latex(r'''[A]_{BL} = \text{Activity concentration in Blood (Residence time in organ = } mTIAC_{org}\text{)}''')
        st.latex(r''' RMECFF = \text{Red Marrow Extracellular Fluid Fraction (physiological: 0.19; range 0.15-0.25)}''')
        st.latex(r''' HCT = \text{Hematocrit (physiological: 0.47; range 0.2-0.6)}''')    
        st.latex(r''' m_{RM\left(patient\right)} = \text{Mass of Red Marrow of patient (ICRP 89 male phantom: 1170g)}''')    
        st.latex(r''' m_{BL\left(patient\right)} = \text{Mass of Blood of patient (ICRP 89 male phantom: 6363g)}''')    
           

    elif bm_method == bm_traino:
        st.write(f' Method described by Traino et al. (2007): DOI 10.1088/0031-9155/52/17/009 (https://iopscience.iop.org/article/10.1088/0031-9155/52/17/009)')
        st.latex(r'''  
            [A]_{RM} = [A]_{BL} * \frac{RMECFF}{1 - HCT} * m_{RM} = [A]_{BL} * RMBLR * m_{RM} 
                ''')
        st.latex(r'''RMBLR = \frac{RMECFF}{1 - HCT} = \text{Red Marrow to Blood Ratio (RMBLR = 1)}''')
        st.latex(r'''[A]_{RM} = \text{Activity concentration in Red Marrow (Residence time in organ = } mTIAC_{org}\text{)}''')
        st.latex(r'''[A]_{BL} = \text{Activity concentration in Blood (Residence time in organ = } mTIAC_{org}\text{)}''')
        st.latex(r''' m_{RM\left(patient\right)} = \text{Mass of Red Marrow of patient (ICRP 89 male phantom: 1170g)}''')    

    elif bm_method == bm_grudzinski:
        st.write(f'Method proposed by Joe Grudzinski (consultant 01/2024)')
        st.latex(r'''  
            [A]_{RM} = [A]_{BL} * \frac{m_{total body, animal}}{m_{total body, human}} * m_{RM}
                ''')
        st.latex(r'''[A]_{RM} = \text{Activity concentration in Red Marrow (Residence time in organ = } mTIAC_{org}\text{)}''')
        st.latex(r'''[A]_{BL} = \text{Activity concentration in Blood (Residence time in organ = } mTIAC_{org}\text{)}''')
        st.latex(r''' m_{total body, animal} = \text{total body weight (animal)}''')
        st.latex(r''' m_{total body, human} = \text{total body weight (human)}''')
        st.latex(r''' m_{RM\left(patient\right)} = \text{Mass of Red Marrow of patient (ICRP 89 male phantom: 1170g)}''')    
           

    elif bm_method == bm_graves:
        st.write(f'Method proposed by Graves (consultant 01/2024) - based on Sgours 1993 with RMBLR = 1')
        st.latex(r'''  
            [A]_{RM} = [A]_{BL} * \frac{RMECFF}{1 - HCT} 
                ''')
        st.latex(r'''[A]_{RM} = \text{Activity concentration in Red Marrow (Residence time in organ = } mTIAC_{org}\text{)}''')
        st.latex(r'''[A]_{BL} = \text{Activity concentration in Blood (Residence time in organ = } mTIAC_{org}\text{)}''')
        st.latex(r'''RMBLR = \frac{RMECFF}{1 - HCT} = \text{Red Marrow to Blood Ratio (RMBLR = 1)}''')

    else:
        st.error('This method is not defined yet')         

def rename_olinda_tissues(olinda_input,source_organ, source_dose):
    if source_organ == 'Gallbladder Wall':
        olinda_input['Gallbladder Contents'] = [source_dose, source_organ]
    elif fnmatch(source_organ,'*LLI Wall*'):
        olinda_input['Left Colon'] = [source_dose, source_organ]
    elif source_organ == 'Stomach Wall':
        olinda_input['Stomach Contents'] = [source_dose, source_organ]
    elif source_organ == 'ULI Wall':
        olinda_input['Right Colon'] = [source_dose, source_organ]
    elif source_organ == 'Heart Wall':
        olinda_input['Heart Contents'] = [float(0.000), '-']
        olinda_input['Heart Wall'] = [source_dose, source_organ]
    elif source_organ == 'Osteogenic Cells':
        olinda_input['Cortical Bone'] = [source_dose, source_organ]
        olinda_input['Trabecular Bone'] = [float(0.000), '-']
    elif source_organ == 'Urinary Bladder Wall':
        olinda_input['Urinary Bladder Contents'] = [source_dose, source_organ]
    else:
        olinda_input[source_organ] = [source_dose, source_organ]
    return olinda_input

def dosimetry_from_hTIAC_org(rayz_id, batch_registration_id, results_scaling_df,molecule_batch_ID,df_sfactors,radioisotope2,doselimits_file):
    st.latex(r'''  absorbed Dose = hTIAC_g * m_{organ} * dosefactor_{organ,nucl} (T <- S)''')
    st.latex(r'''  MTA (MBq) = DoseLimit_{organ} (Gy) / dose_{organ} (Gy/MBq)''')

    targets_fromCSV_list = df_sfactors.columns.tolist() 
    if 'Unnamed: 0' in targets_fromCSV_list: targets_fromCSV_list.remove('Unnamed: 0') 
    targets_fromCSV_list.sort()
    
    sources_fromCSV_list = df_sfactors[df_sfactors.columns[0]].tolist()
    if 'Any bone activity on bone surfaces' in sources_fromCSV_list: sources_fromCSV_list.remove('Any bone activity on bone surfaces') 
    if 'OLINDA - Organ Level INternal Dose Assessment Code (Version 2.2; copyright Vanderbilt University - 2012)' in sources_fromCSV_list: sources_fromCSV_list.remove('OLINDA - Organ Level INternal Dose Assessment Code (Version 2.2; copyright Vanderbilt University - 2012)') 
    # if 'Total Dose Conversion Factors [mSv/MBq-s] - Nuclide:Lu-177 (6.65 d) for ICRP 89 Adult Male' in sources_fromCSV_list: sources_fromCSV_list.remove('Total Dose Conversion Factors [mSv/MBq-s] - Nuclide:Lu-177 (6.65 d) for ICRP 89 Adult Male') 
    total_dose_cf_col = [col_name for col_name in sources_fromCSV_list if fnmatch(str(col_name), '*Total Dose Conversion Factors*')]
    if len(total_dose_cf_col) > 0:
        total_dose_cf_col = total_dose_cf_col[0]
    if total_dose_cf_col in sources_fromCSV_list: sources_fromCSV_list.remove(total_dose_cf_col) 

    sources_fromCSV_list = [x for x in sources_fromCSV_list if str(x) != 'nan'] 

    sources_fromCSV_list.insert(0,'do not use for dosimetry')

    with st.form('Change Input for Olinda'):
        olinda_input_organ = dict()
        for ik, tissue_source in enumerate(results_scaling_df.tissue):
            # Get whole body weight human
            # wb_weight = results_scaling_df[results_scaling_df[key_human_wb_weight].keys()==tissue_source][key_human_wb_weight].iloc[0]

            # dose = calc_dosimetry(tissue_source=tissue_target, sourceslist=results_scaling_df, Sfactors=df_sfactors, hTIAC_org_list=results_scaling_df[htiac_org_key],key_human_organ_weight=key_human_organ_weight,wb_weight=wb_weight,radioisotope_proj_human=radioisotope2)
            if tissue_source == 'tumor':
                olinda_input_organ[tissue_source]   = 'Tumor'
            elif tissue_source == 'bone':
                preselect_olinda_organname_idx = sources_fromCSV_list.index('Total Body')
                olinda_source_organname = st.selectbox(f'Select Olinda input name corresponding to ****{tissue_source}****',options=sources_fromCSV_list, key=f'{tissue_source}-olinda', index=preselect_olinda_organname_idx)
                olinda_input_organ[tissue_source]   = olinda_source_organname
            elif tissue_source == 'tail':
                preselect_olinda_organname_idx = sources_fromCSV_list.index('do not use for dosimetry')
                olinda_source_organname = st.selectbox(f'Select Olinda input name corresponding to ****{tissue_source}****',options=sources_fromCSV_list, key=f'{tissue_source}-olinda', index=preselect_olinda_organname_idx)
                olinda_input_organ[tissue_source]   = olinda_source_organname
            else:
                target_organ_defined, organ_name_options = get_matching_organnames(tissue_source,sources_fromCSV_list)
                preselect_olinda_organname_idx = sources_fromCSV_list.index('Total Body')
                for organ_option in organ_name_options:
                    organ_option = organ_option.strip()
                    if organ_option in sources_fromCSV_list:
                        preselect_olinda_organname_idx = sources_fromCSV_list.index(organ_option)
                        break
                        
                olinda_source_organname = st.selectbox(f'Select Olinda input name corresponding to ****{tissue_source}****',options=sources_fromCSV_list, key=f'{tissue_source}-olinda', index=preselect_olinda_organname_idx)
                olinda_input_organ[tissue_source]   = olinda_source_organname
        st.form_submit_button(f'Confirm Change of Olinda Input')

    results_scaling_df[key_olindainput] = olinda_input_organ.values()
    mapped_olinda_input = results_scaling_df.pivot_table(index=key_olindainput, aggfunc = {htiac_org_key:sum}).reset_index()

    # Need to rename the Tissues from Dose Factor Files since Olinda Kinetics Input has different Organ Names
    dosimetry_input = dict()
    olinda_input = dict()
    for source_organ in sources_fromCSV_list:
        if source_organ in list(mapped_olinda_input[key_olindainput]):
            dosimetry_input[source_organ] = float(list(mapped_olinda_input[htiac_org_key][mapped_olinda_input[key_olindainput]==source_organ])[0])
            olinda_input = rename_olinda_tissues(olinda_input, source_organ = source_organ, source_dose = float(list(mapped_olinda_input[htiac_org_key][mapped_olinda_input[key_olindainput]==source_organ])[0]))
        else:
            dosimetry_input[source_organ] = 0.
            olinda_input = rename_olinda_tissues(olinda_input, source_organ = source_organ, source_dose = float(0.000))

    olinda_input_cleaned = olinda_input
    olinda_input_cleaned.pop('do not use for dosimetry',None)
    olinda_input_df = pd.DataFrame.from_dict(olinda_input_cleaned, orient='index',columns=['# of Disintegrations (MBq-hr/MBq)', 'source organ Olinda'])
    dosimetry_input.pop('do not use for dosimetry',None)
    dosimetry_input_df = pd.DataFrame.from_dict(dosimetry_input, orient='index',columns=['# of Disintegrations (MBq-hr/MBq)'])

    with st.expander('Organ name mapping to Olinda input'):
        st.write(results_scaling_df[[key_olindainput,htiac_org_key]])
        st.write('Resulting aggregated input')
        st.write(mapped_olinda_input)
        st.write('Olinda format')
        st.table(olinda_input_df)
        olinda_input_csv = convert_df_to_csv(olinda_input_df)
        st.download_button(
            label = "Download Olinda input",
            data  = olinda_input_csv,
            file_name =f"{rayz_id}-olinda-input.csv",
            mime = "text/csv",
            key = 'download-olinda-csv'
            )

    doselimit_expander = st.expander(f'Change dose limits to organ')
    results_dosimetry_alltissues = list()
    for ik, tissue_target in enumerate(targets_fromCSV_list):
        saf_vector= df_sfactors[tissue_target].iloc[0:25]
        dosimetry_input_df[tissue_target] = np.array(list(saf_vector), dtype=float)
        dosimetry_input_df[tissue_target+' Dose'] = dosimetry_input_df[tissue_target] * dosimetry_input_df['# of Disintegrations (MBq-hr/MBq)'] * 60 * 60

        # Store dosimetry calculations in dictionary results_dosimetry
        results_dosimetry=dict()
        results_dosimetry['RAYZ ID'] = rayz_id
        results_dosimetry['Batch ID'] = batch_registration_id
        results_dosimetry['Molecule-Batch-ID'] = molecule_batch_ID
        results_dosimetry['Target tissue'] = tissue_target 
        dose =  dosimetry_input_df[tissue_target+' Dose'].sum()
        results_dosimetry[keyword_dose] = dose
        doselimit_organ_names_found, doselimit_organ_names = get_matching_organnames(tissue_target,doselimits_file.Organ)
        if tissue_target == 'Red Mar.':
            doselimit_organ_names_found = True
            doselimit_organ_names = ['Bone Marrow']

        if doselimit_organ_names_found:
            doselimit_from_file = doselimits_file[doselimits_file.Organ == doselimit_organ_names[0]].iloc[0].iloc[1]    # in Gray
            doselimit_from_file = doselimit_expander.number_input(f'{doselimit_organ_names[0]}: MAX {doselimit_from_file} Gy ',value=doselimit_from_file, min_value=0.)#- max. injectable dose: {max_dose_toinject/1000:.2f} GBq')
            max_dose_toinject = doselimit_from_file / (dose*1e-3) #absorbed Dose in mGy/MBq 
            max_dose_toinject_mCi = max_dose_toinject*2.7027E-2     #conversion from MBq to mCi, see https://www.unitconverters.net/radiation-activity/megabecquerel-to-curie.htm
            results_dosimetry[doselimit_keyword] = round(doselimit_from_file,1)
            results_dosimetry[mta_organ_keyword] = round(max_dose_toinject/1000,2)
            doselimit_calculated = True
        results_dosimetry_alltissues.append(results_dosimetry)

    if 'tumor' in results_scaling_df.tissue or 'Tumor' in results_scaling_df.tissue:
        hTIAC_org_tumor = float(list(mapped_olinda_input[htiac_org_key][mapped_olinda_input[key_olindainput] == 'Tumor'])[0])
        tumor_input = True
    else:
        tumor_input = st.checkbox(f'Do you want to enter human residence time for tumor to project human doses?')
        if tumor_input:
            hTIAC_org_tumor = st.number_input(f'Enter residence time for human tumor [MBq*h/MBq]',value=1.0)

    if tumor_input:
        # tumor spheres model: https://jnm.snmjournals.org/content/jnumed/35/1/152.full.pdf
        # disintegrations to MeV conversion: http://www.sfu.ca/~mxchen/phys1021124/P102Lec29B.pdf
        df_tumor_sphere = tumor_spheres_file[tumor_spheres_file.radioisotope==radioisotope2]['dose (mGy/MBq) per disintegration in 310g spherical tumor'].iloc[0]
        dose_total_tumor = hTIAC_org_tumor * df_tumor_sphere
        calc_details = st.expander(f'dosimetry calculation (hTIAC per organ and dosefactors) - target: Tumor - in Olinda: Spheres model - 310g for {radioisotope2}')
        calc_details.write(f'***Tumor***')
        calc_details.write(f'dose: hTIAC_org * df_nucl = {hTIAC_org_tumor:.2e} * {df_tumor_sphere:.2e} = {dose_total_tumor:.2e}')  
        
        results_dosimetry=dict()
        results_dosimetry['RAYZ ID'] = rayz_id
        results_dosimetry['Batch ID'] = batch_registration_id
        results_dosimetry['Molecule-Batch-ID'] = molecule_batch_ID    
        results_dosimetry['Target tissue'] = 'Tumor' 
        results_dosimetry[keyword_dose] = dose_total_tumor 
        results_dosimetry_alltissues.append(results_dosimetry)

    results_doselimits_df = pd.DataFrame.from_dict(results_dosimetry_alltissues)
    if tumor_input:
        results_doselimits_df['Tumor dose at MTA per organ [Gy]'] = results_doselimits_df[mta_organ_keyword] * dose_total_tumor

    # Get dose limiting organ and MTA of dose limiting organ
    mta = results_doselimits_df[mta_organ_keyword].min()
    results_doselimits_df[mta_keyword] = mta
    tissue_DL = results_doselimits_df['Target tissue'][results_doselimits_df[mta_organ_keyword]==mta]
    tissue_DL = tissue_DL.values[0]
    results_doselimits_df['Dose limiting organ'] = tissue_DL

    st.write(f'**Dose limiting organ:** {tissue_DL}')
    st.write(f'MTA = **{mta:.2f} GBq** (dose limiting: {tissue_DL})')
    if tumor_input:
        tumor_exposure_at_DL = mta * dose_total_tumor   # GBq * mGy/MBq 
        results_doselimits_df[keyword_tumor_dose_mta] = tumor_exposure_at_DL
        st.write(f'Tumor exposure at MTA: **{tumor_exposure_at_DL:.2f} Gy** (dose limiting: {tissue_DL})')

    #Get second dose limiting organ
    mta_list = list(results_doselimits_df[mta_organ_keyword])
    mta_list = [x for x in mta_list if not np.isnan(x)]
    mta_list.sort()
    if len(mta_list) >= 1:
        mta_second = mta_list[1]
        tissue_DL_second = results_doselimits_df['Target tissue'][results_doselimits_df[mta_organ_keyword]==mta_second].values[0]
        results_doselimits_df['Second Dose limiting organ'] = tissue_DL_second

        st.write(f'**Second Dose limiting organ:** {tissue_DL_second}')
        st.write(f'Second dose limiting organ: MTA = **{mta_second:.2f} GBq** (dose limiting: {tissue_DL_second})')
        if tumor_input:
            tumor_exposure_at_DL_second = mta_second * dose_total_tumor   # GBq * mGy/MBq 
            results_doselimits_df['Second dose limiting organ: '+keyword_tumor_dose_mta] = tumor_exposure_at_DL_second
            st.write(f'Tumor exposure at second MTA: **{tumor_exposure_at_DL_second:.2f} Gy** (dose limiting: {tissue_DL_second})')

    results_doselimits_df['radioisotope 2'] = radioisotope2
    results_doselimits_df['Dosimetry method'] = irdc_version

    with st.expander('See dosimetry calculation result table (absorbed dose coefficients, MTA, tumor dose):'):
        st.write('****Results absorbed dose coefficients per tissue:****')
        if tumor_input:
            st.write(results_doselimits_df[['Target tissue',keyword_dose,mta_keyword,keyword_tumor_dose_mta,'Dose limiting organ',doselimit_keyword]])
        else:
            st.write(results_doselimits_df[['Target tissue',keyword_dose,mta_keyword,'Dose limiting organ',doselimit_keyword]])
    results_dosimetry_csv = convert_df_to_csv(results_doselimits_df)
    st.download_button(
        label = "Press to Download",
        data  = results_dosimetry_csv,
        file_name ="results-dosimetry.csv",
        mime = "text/csv",
        key = 'download-dosimetry-csv'
        )

    return [olinda_input_df,results_doselimits_df]

def scaling_mTIAC_hTIAC(scaling_method,fitresults,radioisotope1_study,x_fit,data_tissues,cdd_weights=[],blood_key=''):
    with st.expander('Model definition'):   
        definition_scaling_method(scaling_method, 'fitmodel varies throughout organs')

    with st.expander('Change tissue and body weights'):
        col_calc_1, col_calc_2 = st.columns(2) 
        endcolumn = st.container()  

    tissues = list(fitresults.T.index)
    # to exclude the columns parameter and Decay corrected in dataframe fitresults (returns only tissue headers)
    tissues = [element for element in tissues if element not in ['parameter','Decay corrected']]
    radioisotope1_decaycorr = fitresults['Decay corrected'].iloc[0]
    if radioisotope1_decaycorr == 'no decay correction':
        radioisotope1 = radioisotope1_study                           
    else:
        radioisotope1 = radioisotope1_decaycorr
    if radioisotope1 == 'n/a':
        radioisotope1 = st.radio(f"Which radioisotope should be used for scaling?",(isotopes_human_halflives.keys()))     
    lambda_radioisotope1 = np.log(2)/float(isotopes_BioD_halflives[radioisotope1])

    # Store extrapolation information in dictionaries
    results_scaling_alltissues = dict()
    alpha_dict = dict()

    # Iterate through all tissues for individual extrapolation
    for ik, tissue in enumerate(tissues):
        fitmodel = fitresults[tissue]['fitmodel']

        ###################################################################################
        ##### Enter Whole body weights for mouse and human
        ###################################################################################
        if ik < 1:  #to enter whole body weight only once
            if len(cdd_weights) != 0 and tissue is not keyword_projected_bm:
                wholebody_mouse = cdd_weights['organ weight'][cdd_weights['tissue'] =='total body']
            else:
                wholebody_mouse = 25. # gram
            wb_mouse = float(col_calc_1.number_input('body weight mouse [g]', min_value = 0., value = float(wholebody_mouse), step = 1.))
            wb_human = float(col_calc_2.number_input('body weight human [kg]', min_value = 0., value = float(wholebody_human), step = 1.))*1000


        ###################################################################################
        ##### Enter organ weights for mouse and human
        ###################################################################################
        # Get standard organ mass from excel file 'tissue_masses_mouse_file'
        org_weight_mouse_from_file = get_mass_from_file(tissue, tissue_masses_mouse_file, mouse_keyword_tissues, mouse_keyword_masses, 1e-3)
        if len(cdd_weights) != 0 and tissue is not keyword_projected_bm:
            org_weight_mouse_from_cdd = cdd_weights['organ weight'][cdd_weights['tissue'] ==tissue].item()

            # Calculate the difference of reference mouse organ mass and experimental weight from CDD, give an error if diff is >10%
            diff_cdd_standard = abs(org_weight_mouse_from_cdd-org_weight_mouse_from_file)/org_weight_mouse_from_file
            if int(org_weight_mouse_from_file) == 0.:
                org_weight_mouse_from_file = org_weight_mouse_from_cdd
            elif diff_cdd_standard > 0.2: #more than 10% difference of cdd data vs. standard model
                if tissue not in ['blood','muscle']:
                    endcolumn.error(f'Caution here: {tissue} weight in standard mouse is {org_weight_mouse_from_file:.2f}g vs. {org_weight_mouse_from_cdd:.2f}g from CDD Vault upload')
            else:
                org_weight_mouse_from_file = org_weight_mouse_from_cdd

            #     # always use standard weights for blood and muscle since BioD is done on parts of the organ only
            # if tissue == 'blood':   
            #     org_weight_mouse_from_file = get_mass_from_file(tissue, tissue_masses_mouse_file, mouse_keyword_tissues, mouse_keyword_masses, 1e-3)
            # elif tissue == 'muscle':
            #     org_weight_mouse_from_file = get_mass_from_file(tissue, tissue_masses_mouse_file, mouse_keyword_tissues, mouse_keyword_masses, 1e-3)

        org_weight_mouse = float(col_calc_1.number_input(f'{tissue} weight mouse [g]', min_value = 0., max_value = wb_mouse, value = org_weight_mouse_from_file, step = 0.1))
        org_weight_human_from_file = get_mass_from_file(tissue, tissue_masses_human_file, human_keyword_tissues, human_keyword_masses, 1.)
        org_weight_human = float(col_calc_2.number_input(f'{tissue} weight human [g]', min_value = 0., max_value = wb_human, value = org_weight_human_from_file, step = 100.))


        #####################################################################################################
        ##### Do scaling based on selected scaling method and store results in results_scaling dictionary
        #####################################################################################################
        results_scaling = dict()
        #previously fitted results
        tiac = fitresults._get_value(keyword_mtiac_g,tissue)
        datatis_df = pd.DataFrame.from_dict(data_tissues,columns=['x','y'],orient='index')
        mTIAC_organ = tiac*org_weight_mouse     #TIAC per organ
        if isinstance(data_tissues, pd.DataFrame) or isinstance(data_tissues, dict):
            if tissue is keyword_projected_bm:
                if type(blood_key) == list:
                    x_raw=data_tissues[blood_key[0]][0]
                    y_raw_lists = datatis_df.T[blood_key].T['y']
                    y_raw_values = pd.DataFrame(y_raw_lists.values.tolist()).T
                    y_raw = list(y_raw_values[blood_key].mean(axis=1))
                    fitmodel = fitresults[blood_key[0]]['fitmodel']
                else:
                    x_raw=data_tissues[blood_key][0]
                    y_raw=data_tissues[blood_key][1]
                    fitmodel = fitresults[blood_key]['fitmodel']
            else:
                x_raw=data_tissues[tissue][0]
                y_raw=data_tissues[tissue][1]

        # No Scaling
        if scaling_method == no_scaling:
            hTIAC_organ = mTIAC_organ
            hTIAC_g = hTIAC_organ / org_weight_human

            results_scaling[htiac_key] = float(hTIAC_g)
            results_scaling[htiac_org_key] = float(hTIAC_organ)
            results_scaling[keyword_mtiac_g] = float(tiac)
            results_scaling[keyword_mtiac_org] = float(mTIAC_organ)
            results_scaling[key_mouse_organ_weight] = float(org_weight_mouse)
            results_scaling[key_human_organ_weight] = float(org_weight_human)
            results_scaling[key_mouse_wb_weight] = float(wb_mouse)
            results_scaling[key_human_wb_weight] = float(wb_human)

        # Relative Mass Scaling
        elif scaling_method == rel_mass_scal_fda:
            hTIAC_organ = mTIAC_organ * (org_weight_human/wb_human)/(org_weight_mouse/wb_mouse)
            hTIAC_g = hTIAC_organ / org_weight_human

            results_scaling[htiac_key] = float(hTIAC_g)
            results_scaling[htiac_org_key] = float(hTIAC_organ)
            results_scaling[keyword_mtiac_g] = float(tiac)
            results_scaling[keyword_mtiac_org] = float(mTIAC_organ)
            results_scaling[key_mouse_organ_weight] = float(org_weight_mouse)
            results_scaling[key_human_organ_weight] = float(org_weight_human)
            results_scaling[key_mouse_wb_weight] = float(wb_mouse)
            results_scaling[key_human_wb_weight] = float(wb_human)

        elif scaling_method == rel_mass_scal:
            hTIAC_organ = mTIAC_organ * (org_weight_human/wb_human)/(org_weight_mouse/wb_mouse)
            hTIAC_g = hTIAC_organ / org_weight_human

            results_scaling[htiac_key] = float(hTIAC_g)
            results_scaling[htiac_org_key] = float(hTIAC_organ)
            results_scaling[keyword_mtiac_g] = float(tiac)
            results_scaling[keyword_mtiac_org] = float(mTIAC_organ)
            results_scaling[key_mouse_organ_weight] = float(org_weight_mouse)
            results_scaling[key_human_organ_weight] = float(org_weight_human)
            results_scaling[key_mouse_wb_weight] = float(wb_mouse)
            results_scaling[key_human_wb_weight] = float(wb_human)

        # Time Scaling from Beykan et al. (2018)
        elif scaling_method == time_scal:
            x_scal = x_raw*((wb_human/wb_mouse)**(0.25)) # x_scal are the new (extrapolated) timepoints
            if fitmodel == monofit:
                A1, lambda_fit, thalf, y_fit, AUC_calculated, AUC_numerical_extrapolated =  monofit_func(x_scal,y_raw,x_fit)
            elif fitmodel == monoexp_dec_fit:
                A1, lambda_fit, thalf, y_fit, AUC_calculated, AUC_numerical_extrapolated = monofit_dec_func(x_scal,y_raw,x_fit)
            elif fitmodel == bifit:
                A1, lambda1, A2, lambda2, y_fit, AUC_calculated, AUC_numerical_extrapolated = biexp_func(x_scal,y_raw,x_fit)
            elif fitmodel == trapfit:    
                xtrap, ytrap, AUC_calculated = trapezoidal_func(x_scal,y_raw)
                AUC_extrapolated = 0.
            elif fitmodel == linexpfit:
                xtrap, ytrap, AUC_calculated, AUC_extrapolated = trap_linextrapol_func(x_scal,y_raw, halflive = isotopes_BioD_halflives[radioisotope1])
            elif fitmodel == linphysdecay:
                xtrap, ytrap, AUC_calculated,AUC_extrapolated_ratio = trap_physdecayextrapol_func(x_scal,y_raw, halflive = isotopes_BioD_halflives[radioisotope1])
            elif fitmodel == biexpelim:
                ka, ke, cl, y_fit, AUC_calculated,AUC_numerical_extrapolated = biexpelim_func(x_scal,y_raw,x_fit)
            elif fitmodel == biexpelim2:
                lambda1,lambda2,A1, y_fit, AUC_calculated,AUC_numerical_extrapolated = biexpelim_func2(x_scal,y_raw,x_fit)
            elif fitmodel == biexp_abs_elim_two_components:
                lambda1,lambda2,A1,A2, y_fit, AUC_calculated,AUC_numerical_extrapolated = biexpelim_func3(x_scal,y_raw,x_fit)
            elif fitmodel == 'mTIAC_g as input':
                pass
            else:
                st.write('Fit model not implemented for time scaling')
                # results_scaling = dict()
                A1, lambda_fit, thalf, y_fit, AUC_calculated, AUC_numerical_extrapolated = [0.,0.,0.,0.,0.,0.]
                AUC_calculated = 0.

            tiac_adjusted = float(AUC_calculated)/100
            mTIAC_organ_adjusted = tiac_adjusted*org_weight_mouse   #TIAC per organ
            hTIAC_organ = mTIAC_organ_adjusted
            hTIAC_g = hTIAC_organ / org_weight_human

            results_scaling[htiac_key] = float(hTIAC_g)
            results_scaling[htiac_org_key] = float(hTIAC_organ)
            results_scaling[keyword_mtiac_g] = float(tiac)
            results_scaling[keyword_mtiac_org] = float(mTIAC_organ)
            results_scaling[keyword_mtiac_g+' scaled'] = float(tiac_adjusted)
            results_scaling[keyword_mtiac_org+' scaled'] = float(mTIAC_organ_adjusted)
            results_scaling[key_mouse_organ_weight] = float(org_weight_mouse)
            results_scaling[key_human_organ_weight] = float(org_weight_human)
            results_scaling[key_mouse_wb_weight] = float(wb_mouse)
            results_scaling[key_human_wb_weight] = float(wb_human)

        # Combined Time and Mass Scaling from Beykan et al. (2018)
        elif scaling_method == time_mass_scal:
            x_scal = x_raw*((wb_human/wb_mouse)**(1./4.))
            if fitmodel == monofit:
                A1, lambda_fit, thalf, y_fit, AUC_calculated, AUC_numerical_extrapolated =  monofit_func(x_scal,y_raw,x_fit)
            elif fitmodel == monoexp_dec_fit:
                A1, lambda_fit, thalf, y_fit, AUC_calculated, AUC_numerical_extrapolated = monofit_dec_func(x_scal,y_raw,x_fit)
            elif fitmodel == bifit:
                A1, lambda1, A2, lambda2, y_fit, AUC_calculated, AUC_numerical_extrapolated = biexp_func(x_scal,y_raw,x_fit)
            elif fitmodel == trapfit:    
                xtrap, ytrap, AUC_calculated = trapezoidal_func(x_scal,y_raw)
                AUC_extrapolated = 0.
            elif fitmodel == linexpfit:
                xtrap, ytrap, AUC_calculated,AUC_extrapolated = trap_linextrapol_func(x_scal,y_raw, halflive = isotopes_BioD_halflives[radioisotope1])
            elif fitmodel == linphysdecay:
                xtrap, ytrap, AUC_calculated = trap_physdecayextrapol_func(x_scal,y_raw, halflive = isotopes_BioD_halflives[radioisotope1])
            elif fitmodel == biexpelim:
                ka, ke, cl, y_fit, AUC_calculated,AUC_numerical_extrapolated = biexpelim_func(x_scal,y_raw,x_fit)
            elif fitmodel == biexpelim2:
                lambda1,lambda2,A1, y_fit, AUC_calculated,AUC_numerical_extrapolated = biexpelim_func2(x_scal,y_raw,x_fit)
            elif fitmodel == biexp_abs_elim_two_components:
                lambda1,lambda2,A1,A2, y_fit, AUC_calculated,AUC_numerical_extrapolated = biexpelim_func3(x_scal,y_raw,x_fit)

            elif fitmodel == 'mTIAC_g as input':
                pass
            else:
                st.write('Fit model not implemented for time scaling')
                # results_scaling = dict()
                A1, lambda_fit, thalf, y_fit, AUC_calculated, AUC_numerical_extrapolated = [0.,0.,0.,0.,0.,0.]
                AUC_calculated = 0.

            tiac_adjusted = float(AUC_calculated)/100
            mTIAC_organ_adjusted = tiac_adjusted*org_weight_mouse   #TIAC per organ

            hTIAC_organ = mTIAC_organ_adjusted * (org_weight_human/wb_human)/(org_weight_mouse/wb_mouse)
            hTIAC_g = hTIAC_organ / org_weight_human

            results_scaling[htiac_key] = float(hTIAC_g)
            results_scaling[htiac_org_key] = float(hTIAC_organ)
            results_scaling[keyword_mtiac_g] = float(tiac)
            results_scaling[keyword_mtiac_org] = float(mTIAC_organ)
            results_scaling[keyword_mtiac_g+' scaled'] = float(tiac_adjusted)
            results_scaling[keyword_mtiac_org+' scaled'] = float(mTIAC_organ_adjusted)
            results_scaling[key_mouse_organ_weight] = float(org_weight_mouse)
            results_scaling[key_human_organ_weight] = float(org_weight_human)
            results_scaling[key_mouse_wb_weight] = float(wb_mouse)
            results_scaling[key_human_wb_weight] = float(wb_human)

        # Combined Time and Mass Scaling from FDA guideline
        elif scaling_method == time_mass_fda:
            if fitmodel == monofit:
                A1, lambda_fit, thalf, y_fit, AUC_calculated,AUC_numerical_extrapolated =  monofit_func(x_raw,y_raw,x_fit)
                lambda_fit = lambda_fit * (org_weight_human/wb_human)/(org_weight_mouse/wb_mouse)
                AUC_scaled = (A1/lambda_fit)    # from book M.Stabin: Basic Principles of Internal Dosimetry Calculations, page 7
            elif fitmodel == monoexp_dec_fit:
                A1, lambda_fit, thalf, y_fit, AUC_calculated, AUC_numerical_extrapolated = monofit_dec_func(x_raw,y_raw,x_fit)
                lambda_fit = lambda_fit * (org_weight_human/wb_human)/(org_weight_mouse/wb_mouse)
                AUC_scaled = (A1/lambda_fit)    # from book M.Stabin: Basic Principles of Internal Dosimetry Calculations, page 7                
            elif fitmodel == bifit:
                A1, lambda1, A2, lambda2, y_fit, AUC_calculated, AUC_numerical_extrapolated = biexp_func(x_raw,y_raw,x_fit)
                lambda1 = lambda1 * (org_weight_human/wb_human)/(org_weight_mouse/wb_mouse)
                lambda2 = lambda2 * (org_weight_human/wb_human)/(org_weight_mouse/wb_mouse)
                if lambda1 < 1e-4:
                    AUC_scaled = (A2/lambda2)
                elif lambda2 < 1e-4:
                    AUC_scaled = (A1/lambda1)
                else:
                    AUC_scaled = (A1/lambda1) + (A2/lambda2)
            elif fitmodel == biexpelim_func2:
                lambda1, lambda2, A1, y_fit, AUC_calculated,AUC_numerical_extrapolated = biexpelim_func2(x_raw,y_raw,x_fit)
                lambda1 = lambda1 * (org_weight_human/wb_human)/(org_weight_mouse/wb_mouse)
                lambda2 = lambda2 * (org_weight_human/wb_human)/(org_weight_mouse/wb_mouse)
                if lambda1 < 1e-4:
                    AUC_scaled = (A2/lambda2)
                elif lambda2 < 1e-4:
                    AUC_scaled = (A1/lambda1)
                else:
                    AUC_scaled = (A1/lambda1) + (A2/lambda2)
            elif fitmodel == biexp_abs_elim_two_components:
                lambda1,lambda2,A1,A2, y_fit, AUC_calculated,AUC_numerical_extrapolated = biexpelim_func3(x_scal,y_raw,x_fit)
                lambda1 = lambda1 * (org_weight_human/wb_human)/(org_weight_mouse/wb_mouse)
                lambda2 = lambda2 * (org_weight_human/wb_human)/(org_weight_mouse/wb_mouse)
                if lambda1 < 1e-4:
                    AUC_scaled = (A2/lambda2)
                elif lambda2 < 1e-4:
                    AUC_scaled = (A1/lambda1)
                else:
                    AUC_scaled = (A1/lambda1) + (A2/lambda2)
            elif fitmodel == 'mTIAC_g as input':
                pass
            else:
                st.write(f'Fit model not implemented for time scaling. (Please select between {monofit},{monoexp_dec_fit},{bifit} and {biexpelim_func2})')
                A1, lambda_fit, thalf, y_fit, AUC_calculated, AUC_numerical_extrapolated = [0.,0.,0.,0.,0.,0.]
                AUC_calculated = 0.

            hTIAC_organ = mTIAC_organ * (org_weight_human/wb_human)/(org_weight_mouse/wb_mouse)
            hTIAC_g = hTIAC_organ / org_weight_human

            results_scaling[htiac_key] = float(hTIAC_g)
            results_scaling[htiac_org_key] = float(hTIAC_organ)
            results_scaling[keyword_mtiac_g] = float(tiac)
            results_scaling[keyword_mtiac_org] = float(mTIAC_organ)
            results_scaling[key_mouse_organ_weight] = float(org_weight_mouse)
            results_scaling[key_human_organ_weight] = float(org_weight_human)
            results_scaling[key_mouse_wb_weight] = float(wb_mouse)
            results_scaling[key_human_wb_weight] = float(wb_human)

        #'Method #5 in https://doi.org/10.1155/2019/6438196 (BUT misreferenced in this publication - has to be revisited)'
        elif scaling_method == allom_scal:
            if tissue == 'kidney':
                b_scaling_beykan = 0.85
            elif tissue == 'liver':
                b_scaling_beykan = 0.92
            else:
                b_scaling_beykan = 2.

            hTIAC_organ = mTIAC_organ * ((wb_human/wb_mouse)**(b_scaling_beykan-1))
            hTIAC_g = hTIAC_organ / org_weight_human

            results_scaling[htiac_key] = float(hTIAC_g)
            results_scaling[htiac_org_key] = float(hTIAC_organ)
            results_scaling[keyword_mtiac_g] = float(tiac)
            results_scaling[keyword_mtiac_org] = float(mTIAC_organ)
            results_scaling[key_mouse_organ_weight] = float(org_weight_mouse)
            results_scaling[key_human_organ_weight] = float(org_weight_human)
            results_scaling[key_mouse_wb_weight] = float(wb_mouse)
            results_scaling[key_human_wb_weight] = float(wb_human) 
            results_scaling['allometric scaling coefficient'] = float(b_scaling_beykan)-1. 

        elif scaling_method == metabol_scal:
            hTIAC_organ = mTIAC_organ * ((wb_human/wb_mouse)**(1/4)) * (org_weight_human/wb_human)/(org_weight_mouse/wb_mouse)
            hTIAC_g = hTIAC_organ / org_weight_human

            results_scaling[htiac_key] = float(hTIAC_g)
            results_scaling[htiac_org_key] = float(hTIAC_organ)
            results_scaling[keyword_mtiac_g] = float(tiac)
            results_scaling[keyword_mtiac_org] = float(mTIAC_organ)
            results_scaling[key_mouse_organ_weight] = float(org_weight_mouse)
            results_scaling[key_human_organ_weight] = float(org_weight_human)
            results_scaling[key_mouse_wb_weight] = float(wb_mouse)
            results_scaling[key_human_wb_weight] = float(wb_human)

        elif scaling_method == alpha_scal: 
            tissue_results_key = tissue
            key_tissue_preset = tissue
            if tissue is keyword_projected_bm:
                if type(blood_key) == list:
                    tissue_results_key = blood_key[0]
                    key_tissue_preset = tissue_results_key
                else:
                    tissue_results_key = blood_key

            # key_tissue_preset = tissue
            if tissue == 'tumor':
                tumor_alpha_fromBlood = st.checkbox(f'Do you want to use alpha from ****blood for scaling of tumor****?',value=False)
                if tumor_alpha_fromBlood:
                    if 'blood calculated' in alpha_dict.keys():
                        key_tissue_preset = tissues.index('blood')
                        key_blood = 'blood'
                    else:
                        key_tissue_preset = 0
                        key_blood = st.selectbox(f'Which tissue is surrogate for blood?', options=tissues,index=key_tissue_preset)

            elif tissue == keyword_projected_bm:
                bm_fit = fitresults[keyword_projected_bm].T['fitmodel']
                bm_keyword_fit = bm_fit.split('<')[1]
                bm_keyword_fit = bm_keyword_fit.split('>')[0]
                if fnmatch(bm_keyword_fit,'''['*'''):
                    bm_keyword_fit = bm_keyword_fit.split(',')[0].split("'")[1]
                if bm_keyword_fit+' calculated' in alpha_dict.keys():
                    key_tissue_preset = tissues.index(bm_keyword_fit)
                    key_blood = bm_keyword_fit
                elif 'blood calculated' in alpha_dict.keys():
                    key_tissue_preset = tissues.index('blood')
                    key_blood = 'blood'
                else:
                    key_tissue_preset = 0
                key_blood = st.selectbox(f'Which tissue is surrogate for blood?', options=tissues,index=key_tissue_preset)

                fitmodel = str(fitresults[key_blood].T['fitmodel'])

            alpha_expander = st.expander(f'Alpha calculated for ****{tissue}**** from {tissue_results_key} fit results - change here if you want')  

            if fitmodel == monofit:
                lambda_fit = fitresults._get_value('lambda_fit',tissue_results_key)

            elif fitmodel == monoexp_dec_fit:
                lambda_fit = fitresults._get_value('lambda_fit',tissue_results_key)

            elif fitmodel == bifit or fitmodel == biexpelim2 or fitmodel == biexp_abs_elim_two_components:
                lambda1 = fitresults._get_value('lambda1',tissue_results_key)
                lambda2 = fitresults._get_value('lambda2',tissue_results_key)
                # with st.expander("Three principles about effective half-time"):
                #     ''' Three principles always apply to the effective half-time:'''
                #     '''1. The units for both half-times have to be the same.'''
                #     '''2. The effective half-time must always be shorter than the smaller of the biological and physical half-times.'''
                #     '''3. As one half-time becomes very large relative to the other, the effective half-time converges to the smaller of the two.'''
                #     '''from Stabin, The Practice of Internal Dosimetry in Nuclear Medicine, page 7
                #     '''
                if lambda1 < lambda2:
                    lambda_fit = lambda1
                else:
                    lambda_fit = lambda2
            elif fitmodel == biexpelim:
                ka_m = fitresults._get_value('ka',tissue_results_key)
                ke_m = fitresults._get_value('ke',tissue_results_key)
                cl = fitresults._get_value('cl',tissue_results_key)
                ke_h = lambda_radioisotope1 - (ke_m-lambda_radioisotope1)/4.7  # o_lambda = lambda_radioisotope1, p_lambda = lambda_radioisotope2
                ka_h = (ka_m+lambda_radioisotope1)/4.7+lambda_radioisotope1
            elif fitmodel == trapfit or fitmodel == linexpfit or fitmodel == linphysdecay:    
                st.error('Please select an exponential fit model for alpha scaling')
                AUC_calculated = 0.
                lambda_fit = 0.
                results_scaling = dict()
            elif fitmodel == 'mTIAC_g as input':
                pass
            else:
                st.write(f'Fit model ***{fitmodel}*** not implemented for alpha scaling')
                results_scaling = dict()
                AUC_calculated = 0.
            
            if fitmodel == biexpelim:
                alpha_calculated = ((ke_m*ka_m)/(ka_m-ke_m))/((ke_h*ka_h)/(ka_h-ke_h)) # from David Huang, reference?
                alpha_expander.latex(r'''  \alpha = \frac{\frac{ke_m*ka_m}{ka_m-ke_m}}{\frac{ke_h*ka_h}{ka_h-ke_h}} ''')

            else:
                if fitmodel == trapfit or fitmodel == linexpfit or fitmodel == linphysdecay or fitmodel == 'mTIAC_g as input':
                    alpha_calculated = 3.7
                elif radioisotope1_decaycorr == 'no decay correction':
                    alpha_calculated = lambda_fit/(lambda_radioisotope1+(lambda_fit-lambda_radioisotope1)/4.7)
                    alpha_expander.error(f'alpha is calculated with isotope {radioisotope1} used in study')
                else:
                    alpha_calculated = lambda_fit/(lambda_radioisotope1+(lambda_fit-lambda_radioisotope1)/4.7)

            alpha_dict[tissue+' calculated'] = alpha_calculated
            alpha_expander.latex(f'{tissue}: '+r'''  \alpha = '''+f'{alpha_calculated:.2f}')
            if tissue == 'tumor':
                if tumor_alpha_fromBlood:
                    alpha_tumor = alpha_dict[key_blood+' calculated']
                    alpha_expander.write(f'alpha for blood ({key_blood}) is used for tumor scaling: alpha = {alpha_tumor:.2f}')
                else:
                    alpha_tumor = alpha_calculated
                    alpha_expander.write(f'alpha(tumor) is used for scaling: alpha = {alpha_tumor:.2f}')

            # 'You can change alpha (in the range 0-20) if you want:'
            placeholder_key = 5555
            if tissue == 'tumor':
                alpha_calculated = alpha_tumor
            try:
                new_alpha = alpha_expander.number_input(f'alpha for {tissue}:', min_value = 0., max_value = 20., value = alpha_calculated, step= 1., key = f'{placeholder_key}-{tissue}')
            except:
                new_alpha = alpha_expander.number_input(f'alpha for {tissue}:', min_value = 0., max_value = 20., value = 3.7, step= 1., key = f'{placeholder_key}-{tissue}-standard')
            alpha_dict[tissue+' changed'] = new_alpha
            if alpha_expander.button(f'Reset alpha for {tissue}'):
                new_alpha = alpha_dict[tissue+' calculated']
                placeholder_key += 1
                new_alpha = alpha_expander.number_input(f'alpha for {tissue}:', min_value = 0., max_value = 20., value = alpha_dict[tissue+' calculated'], step= 1.,  key = f'{placeholder_key}-{tissue}')

            alpha_expander.write(f'****Scaling for {tissue} will be performed with alpha = {new_alpha:.2f}****')

            if tissue == 'kidney':
                hTIAC_organ = new_alpha * mTIAC_organ  # Stephen Graves
            else:
                hTIAC_organ = new_alpha * mTIAC_organ * (org_weight_human/wb_human)/(org_weight_mouse/wb_mouse)


            hTIAC_g = hTIAC_organ / org_weight_human
            results_scaling[htiac_key] = float(hTIAC_g)
            results_scaling[htiac_org_key] = float(hTIAC_organ)
            results_scaling['alpha used for scaling'] = float(new_alpha)
            results_scaling['alpha calculated'] = float(alpha_dict[tissue+' calculated'])
            results_scaling[keyword_mtiac_g] = float(tiac)
            results_scaling[keyword_mtiac_org] = float(mTIAC_organ)
            results_scaling[key_mouse_organ_weight] = float(org_weight_mouse)                    
            results_scaling[key_human_organ_weight] = float(org_weight_human)
            results_scaling[key_mouse_wb_weight] = float(wb_mouse)
            results_scaling[key_human_wb_weight] = float(wb_human)
            results_scaling[keyword_Radioisotope1] = radioisotope1

        else:
            st.write('Not defined yet')

        results_scaling['decay fit method'] = fitmodel
        results_scaling['scaling method'] = scaling_method
        results_scaling['tissue'] = tissue

        results_scaling_alltissues[tissue] = results_scaling

    results_scaling_df = pd.DataFrame.from_dict(results_scaling_alltissues, orient='index')
    with st.expander('**See results**'):
        st.table(results_scaling_df)
    results_scaling_csv = convert_df_to_csv(results_scaling_df)
    st.download_button(
        label = "Press to Download",
        data  = results_scaling_csv,
        file_name ="results-scaling.csv",
        mime = "text/csv",
        key = 'download-scaling-csv'
        )

    return results_scaling_df

def plot_rawdata_bioD(keywords_list,tissues,data_input):
    with st.expander("plot rawdata"):
        col4, col5 = st.columns(2)
        fig_rawdata = go.Figure()
        div = 260/len(keywords_list)
        for ik, tissue in enumerate(tissues): #https://plotly.com/python/colorscales/
            ik_rgb = ik*div
            injdose_keyword_tissue = tissue
            fig_rawdata.add_trace(go.Bar(
                x=data_input[time_keyword],
                y=data_input[injdose_keyword_tissue],
                name=tissue,
                marker={'color': f'rgb({ik_rgb},{ik_rgb},{ik_rgb})'}
            ))
            fig_rawdata.add_trace(go.Scatter(
                x=data_input[time_keyword],
                y=data_input[injdose_keyword_tissue],
                name=tissue
            ))     
        fig_rawdata.update_layout(
            yaxis_title=f"[%ID/g]",
            xaxis_title=f"{time_keyword} [{time_unit}]",
        )
        col4.plotly_chart(fig_rawdata)

def select_scaling(calculate_scaling, scaling_options=scaling_options):
    if st.checkbox('Show comparison of extrapolation/scaling methods'):
        st.image(Image.open('./pics/comparison-scaling-methods.png'), caption='Comparison with rawdata from Beykan et al. - yellow indicates real measured hTIAC')  

    with st.form(key='options_scaling'):
        scaling_method = st.radio(
                f"Scaling method for mTIAC -> hTIAC",
                (scaling_options))

        submitted = st.form_submit_button(f'click to calculate')
        if submitted:
            calculate_scaling = True
    return [scaling_method, calculate_scaling]

def fit_quality_calc(y_raw,y_fit_quality):                                
    d = np.array(list(y_raw)) - np.array(list(y_fit_quality))    #difference between arrays
    mse_f = np.mean(d**2)                               #mean squared error
    mae_f = np.mean(abs(d))                             #mean absolute error
    rmse_f = np.sqrt(mse_f)                             #root mean squared error
    r2_f = 1-(sum(d**2)/sum((np.array(list(y_raw))-np.mean(np.array(list(y_raw))))**2))         # r squared/Coefficient of determination
    return [mse_f,mae_f,rmse_f,r2_f]

def fit_decay_fitmodel(data_input,tissues,time_keyword,injdose_keyword,radioisotope1,decay_corrected):
    same_fitmodel_all = st.checkbox('Same fitmodel for all tissues?',value=False, key='same_fitmodel_checkbox')
    if same_fitmodel_all:
        fitmodel_presubmit = st.selectbox('Which fitmodel for all?',options= fitmodel_options)
    else:
        fitmodel_presubmit='free to adjust each individual organ'

    AUC_calculated = 'not calculated'

    data_tissues=dict()      
    fitresults_all = dict()
    all_figures_fit= dict()

    # iterate through all tissues to fit each tissue individually
    for ik, injdose_keyword_tissue in enumerate(tissues):
        tissue = injdose_keyword_tissue.replace(injdose_keyword, '')
        with st.form(f'Fit Time Activity Curve for {tissue}'):
            # get y values for selected tissue
            y_raw = data_input[injdose_keyword_tissue]
            # remove any empty cells (and adapt the timepoints/x values accordingly)
            if y_raw.isnull().values.any():
                new_Df = pd.concat([data_input[time_keyword],y_raw],axis=1)
                new_Df=new_Df.dropna()
                x_raw = new_Df[new_Df.keys()[0]]
                y_raw = new_Df[new_Df.keys()[1]]
            else:
                x_raw = data_input[time_keyword]

            # create x values in the range from min(x) to max(x)
            x_fit = np.linspace(np.min(x_raw),np.max(x_raw), 100)
            data_tissues[tissue] = [x_raw,y_raw]
            fitresults_tissue = dict()

            # create the figure that will show x and y input (and add the fit later)
            ik_rgb = 0 #ik*div
            fig_fit = go.Figure()
            fig_fit.add_trace(go.Scatter(
                x=x_raw,
                y=y_raw,
                name=tissue,
                mode='markers',
                marker={'color': f'rgb({ik_rgb},{ik_rgb},{ik_rgb})'}
            ))

            # Create a header with anchor, so you can create a hyperlink later in the app if you have problems with the fit 
            st.header(f'{tissue}', anchor = str('anchor-fitdecay-')+str(tissue))

                # Try all fits and select the fit with R squared closest to 1.
            if fitmodel_presubmit == 'free to adjust each individual organ':
                R_squared_testfits = {}
                trap_integrated = False
                fit_successful_test = False
                for fitmodel_test in fitmodel_options:
                    fit_successful_test = False
                    if fitmodel_test == monofit:
                        try: 
                            A1, lambda_fit, thalf, y_fit, AUC_calculated, AUC_extrapolated =  monofit_func(x_raw,y_raw,x_fit)   
                            A1, lambda_fit, thalf, y_fit_quality, AUC_calculated, AUC_extrapolated =  monofit_func(x_raw,y_raw,x_raw)   
                            fit_variables = 2
                            fit_successful_test = True
                        except:
                            r2_f = -1.

                    elif fitmodel_test == monoexp_dec_fit:
                        try:
                            A1, lambda_fit, thalf, y_fit, AUC_calculated, AUC_extrapolated = monofit_dec_func(x_raw,y_raw,x_fit)             
                            A1, lambda_fit, thalf, y_fit_quality, AUC_calculated, AUC_extrapolated = monofit_dec_func(x_raw,y_raw,x_raw) 
                            fit_variables = 2
                            fit_successful_test = True
                        except:
                            r2_f = -1.

                    elif fitmodel_test == bifit:
                        fit_variables = 4
                        if len(x_raw) < fit_variables:
                            fit_successful_test = False
                        else:
                            try:
                                A1, lambda1, A2, lambda2, y_fit, AUC_calculated, AUC_extrapolated = biexp_func(x_raw,y_raw,x_fit)
                                A1, lambda1, A2, lambda2, y_fit_quality, AUC_calculated, AUC_extrapolated = biexp_func(x_raw,y_raw,x_raw)
                                fit_successful_test = True
                            except:
                                r2_f = -1.

                    elif fitmodel_test == biexpelim:
                        fit_variables = 3
                        if len(x_raw) < fit_variables:
                            fit_successful_test = False
                        else:
                            try:
                                ka, ke, cl, y_fit, AUC_calculated,AUC_extrapolated = biexpelim_func(x_raw,y_raw,x_fit)
                                ka, ke, cl, y_fit_quality, AUC_calculated,AUC_extrapolated = biexpelim_func(x_raw,y_raw,x_raw)
                                fit_successful_test = True
                            except:
                                r2_f = -1.

                    elif fitmodel_test == biexpelim2:
                        fit_variables = 3
                        if len(x_raw) < fit_variables:
                            fit_successful_test = False
                        else:
                            try:
                                lambda1, lambda2, A1, y_fit, AUC_calculated,AUC_extrapolated = biexpelim_func2(x_raw,y_raw,x_fit)
                                lambda1, lambda2, A1, y_fit_quality, AUC_calculated,AUC_extrapolated = biexpelim_func2(x_raw,y_raw,x_raw)
                                fit_successful_test = True
                            except:
                                r2_f = -1.

                    elif fitmodel_test == biexp_abs_elim_two_components:
                        fit_variables = 4
                        if len(x_raw) < fit_variables:
                            fit_successful_test = False
                        else:
                            try:
                                lambda1, lambda2, A1, A2, y_fit, AUC_calculated,AUC_extrapolated = biexpelim_func3(x_raw,y_raw,x_fit)
                                lambda1, lambda2, A1, A2, y_fit_quality, AUC_calculated,AUC_extrapolated = biexpelim_func3(x_raw,y_raw,x_raw)
                                fit_successful_test = True
                            except:
                                r2_f = -1.

                    elif fitmodel_test == trapfit:    #https://personal.math.ubc.ca/~pwalls/math-python/integration/trapezoid-rule/
                        xtrap, ytrap, AUC_calculated = trapezoidal_func(x_raw,y_raw)
                        AUC_extrapolated = 0.
                        trap_integrated = True

                    elif fitmodel_test == linexpfit:
                        xtrap, ytrap, AUC_calculated,AUC_extrapolated = trap_linextrapol_func(x_raw,y_raw, halflive = isotopes_BioD_halflives[radioisotope1], write_error_warning = False)
                        trap_integrated = True

                    elif fitmodel_test == linphysdecay:
                        xtrap, ytrap, AUC_calculated,AUC_extrapolated = trap_physdecayextrapol_func(x_raw,y_raw, halflive = isotopes_BioD_halflives[radioisotope1])
                        trap_integrated = True

                    if fit_successful_test:
                        mse_f,mae_f,rmse_f,r2_f = fit_quality_calc(y_raw,y_fit_quality)  # get mean squared error, mean absolute error, root mean squared error and r squared
                        n_obs = len(y_raw)  #number of observations
                        p_var = fit_variables  #number of independent variables in model
                        try:
                            r2_adj = 1 - (1- r2_f) * ((n_obs-1)/(n_obs - p_var - 1)) # R squared adjusted is normalized with number of independent variables of model
                        except:
                            r2_adj = 0.
                    elif trap_integrated:
                        r2_adj,r2_f,n_obs,p_var = [-1.,-1.,-1.,-1.]
                    else: 
                        r2_adj,r2_f,n_obs,p_var = [-1.,-1.,-1.,-1.]
                        AUC_calculated = None
                        AUC_extrapolated = None

                    R_squared_testfits[fitmodel_test] = [round(r2_adj,4),round(r2_f,4),p_var,n_obs,AUC_calculated,AUC_extrapolated]

                R_squared_testfits_df = pd.DataFrame.from_dict(R_squared_testfits, orient='index', columns=['R squared adjusted','R squared', 'ind. var. in fit','n# observations','AUC','AUC extrapolated'])
                
                # to avoid chosing R squared adjusted if only 3 datapoints are passed (in that case, R squared adjusted can't be determined for monoexponential fits)
                if len(x_raw) < 4:
                    best_test_fit = list(R_squared_testfits_df[R_squared_testfits_df['R squared']==max(R_squared_testfits_df['R squared'])].index)[0]
                else:
                    best_test_fit = list(R_squared_testfits_df[R_squared_testfits_df['R squared adjusted']==max(R_squared_testfits_df['R squared adjusted'])].index)[0]

                with st.expander(f'best fit: {best_test_fit}; see all fit results'):
                    st.write('R squared adjusted is normalized with number of independent variables of model')
                    st.latex(r'''R^2_{adj} = 1 - (1 - R^2) * (n-1)/(n - p - 1)''')
                    st.write(f'n = number of observations')
                    st.write(f'p = number of independent variables in model')
                    st.write(R_squared_testfits_df)

                # select the best fit as fitmodel (index = best test fit)
                fitmodel = st.selectbox(f'What model do you want to use for fitting {injdose_keyword_tissue}?', options=fitmodel_options,key=ik, index = fitmodel_options.index(best_test_fit))
            else:
                # if the same fitmodel is used for all tissues: presubmitted fitmodel will be selected for all tissues (index = select_index)
                select_index = fitmodel_options.index(fitmodel_presubmit)
                fitmodel = st.selectbox(f'What model do you want to use for fitting {injdose_keyword_tissue}?', options=fitmodel_options,key=ik, index = select_index)


            st.write('Fit: %s '%(fitmodel))
            col_fit_1, col_fit_3 = st.columns(2)
            plot_fitresult = True
            if fitmodel == monofit:
                if len(x_raw) < 3:
                    st.error('At least 3 datapoints are needed for fitting! ')
                    plot_fitresult = False
                    fit_successful = False
                else:
                    try:
                        A1, lambda_fit, thalf, y_fit, AUC_calculated, AUC_extrapolated =  monofit_func(x_raw,y_raw,x_fit)   
                        A1, lambda_fit, thalf, y_fit_quality, AUC_calculated, AUC_extrapolated =  monofit_func(x_raw,y_raw,x_raw)   
                        fit_successful = True
                    except:
                        fit_successful = False
                        st.error('Fit not converging')
                    col_fit_1.latex(r'''  A(t) = A_1 * e^{-\lambda_{eff,1}  * t}''')
                    if fit_successful:
                        col_fit_1.latex(r'''fit result:  c(t) = '''+f'{A1:.2f}'+r''' * e^{-'''+f'{lambda_fit:.2f}'+r''' * t}''')
                        col_fit_1.latex(       r''' A_1 = ''' + f'{A1:.2f}')
                        col_fit_1.latex(       r''' \lambda_{eff,1} = ''' + f'{lambda_fit:.2f}   [1/h]')
                        col_fit_1.latex(       r''' T_{1/2} = ''' + f'{thalf:.2f}')
                        plot_fitresult = True
                        fig_fit.add_trace(go.Scatter(x=x_fit,y=y_fit,fill='tozeroy',mode='none')) #mode = 'markers+lines' or delete 

                        fitresults_tissue['A1'] = round(A1,4)
                        fitresults_tissue['lambda_fit'] = round(lambda_fit,4)
                        fitresults_tissue['thalf'] = round(thalf,4)

            elif fitmodel == monoexp_dec_fit:
                if len(x_raw) < 3:
                    st.error('At least 3 datapoints are needed for fitting! ')
                    plot_fitresult = False
                    fit_successful = False
                else:
                    try:
                        A1, lambda_fit, thalf, y_fit, AUC_calculated, AUC_extrapolated = monofit_dec_func(x_raw,y_raw,x_fit)             
                        A1, lambda_fit, thalf, y_fit_quality, AUC_calculated, AUC_extrapolated = monofit_dec_func(x_raw,y_raw,x_raw)             
                        fit_successful = True
                    except:
                        fit_successful = False
                        st.error('Fit not converging')

                    col_fit_1.latex(r'''  A(t) = A_1 * e^{-\lambda_{eff,1}  * t}''')
                    if fit_successful:
                        col_fit_1.latex(r'''fit result:  c(t) = '''+f'{A1:.2f}'+r''' * e^{-'''+f'{lambda_fit:.2f}'+r''' * t}''')
                        col_fit_1.latex(       r''' A_1 = ''' + f'{A1:.2f}')
                        col_fit_1.latex(       r''' \lambda_{eff,1} = ''' + f'{lambda_fit:.2f}   [1/h]')
                        col_fit_1.latex(       r''' T_{1/2} = ''' + f'{thalf:.2f}')
                        plot_fitresult = True
                        fig_fit.add_trace(go.Scatter(x=x_fit,y=y_fit,fill='tozeroy',mode='none')) #mode = 'markers+lines' or delete 

                        fitresults_tissue['A1'] = round(A1,4)
                        fitresults_tissue['lambda_fit'] = round(lambda_fit,4)
                        fitresults_tissue['thalf'] = round(thalf,4)

            elif fitmodel == bifit:
                if len(x_raw) < 4:
                    st.error('At least 4 datapoints are needed for fitting! ')
                    plot_fitresult = False
                    fit_successful = False
                else:
                    try:
                        A1, lambda1, A2, lambda2, y_fit, AUC_calculated, AUC_extrapolated = biexp_func(x_raw,y_raw,x_fit)
                        A1, lambda1, A2, lambda2, y_fit_quality, AUC_calculated, AUC_extrapolated = biexp_func(x_raw,y_raw,x_raw)
                        fit_successful = True
                    except:
                        st.error('Fit not converging')
                        fit_successful = False
                        plot_fitresult = False
                    col_fit_1.latex(r'''  
                        A(t) = A_1  * e^{-\lambda_{eff,1} * t} + A_2  * e^{-\lambda_{eff,2} * t}
                        ''')
                    if fit_successful:
                        col_fit_1.latex(       r''' A_1 = ''' + f'{A1:.2f}')
                        col_fit_1.latex(       r''' \lambda_{eff,1} = ''' + f'{lambda1:.2f}')
                        col_fit_1.latex(       r''' A_2 = ''' + f'{A2:.2f}')
                        col_fit_1.latex(       r''' \lambda_{eff,2} = ''' + f'{lambda2:.2f}')
                        plot_fitresult = True
                        fig_fit.add_trace(go.Scatter(x=x_fit,y=y_fit,fill='tozeroy',mode='none')) #mode = 'markers+lines' or delete 
                        
                        fitresults_tissue['A1'] = round(A1,4)
                        fitresults_tissue['lambda1'] = round(lambda1,4)
                        fitresults_tissue['A2'] = round(A2,4)
                        fitresults_tissue['lambda2'] = round(lambda2,4)

            elif fitmodel == trapfit:    #https://personal.math.ubc.ca/~pwalls/math-python/integration/trapezoid-rule/
                xtrap, ytrap, AUC_calculated = trapezoidal_func(x_raw,y_raw)
                y_fit_quality = y_raw
                plot_fitresult = False
                fig_fit.add_trace(go.Scatter(x=xtrap,y=ytrap,fill='tozeroy',mode='none')) #mode = 'markers+lines' or delete 
                fit_successful = True
                AUC_extrapolated = 0.

            elif fitmodel == linexpfit:
                xtrap, ytrap, AUC_calculated,AUC_extrapolated = trap_linextrapol_func(x_raw,y_raw, halflive = isotopes_BioD_halflives[radioisotope1])
                y_fit_quality = y_raw
                plot_fitresult = False
                fig_fit.add_trace(go.Scatter(x=xtrap,y=ytrap,fill='tozeroy',mode='none')) #mode = 'markers+lines' or delete 
                fit_successful = True

            elif fitmodel == linphysdecay:
                xtrap, ytrap, AUC_calculated,AUC_extrapolated = trap_physdecayextrapol_func(x_raw,y_raw, halflive = isotopes_BioD_halflives[radioisotope1])
                y_fit_quality = y_raw
                plot_fitresult = False
                fig_fit.add_trace(go.Scatter(x=xtrap,y=ytrap,fill='tozeroy',mode='none')) #mode = 'markers+lines' or delete 
                x_decay = np.linspace(list(x_raw)[-1],list(x_raw)[-1]*20,num=20)
                x_decay_fit = x_decay - list(x_raw)[-1]
                y_decay = list(y_raw)[-1] * np.exp(-x_decay_fit/(isotopes_BioD_halflives[radioisotope1]/np.log(2)))
                fig_fit.add_trace(go.Scatter(x=x_decay,y=y_decay,fill='tozeroy',mode='none')) #mode = 'markers+lines' or delete 
                fit_successful = True

            elif fitmodel == biexpelim:
                if len(x_raw) < 3:
                    st.error('At least 3 datapoints are needed for fitting! ')
                    plot_fitresult = False
                    fit_successful = False
                else:
                    try:
                        ka, ke, cl, y_fit, AUC_calculated,AUC_extrapolated = biexpelim_func(x_raw,y_raw,x_fit)
                        ka, ke, cl, y_fit_quality, AUC_calculated,AUC_extrapolated = biexpelim_func(x_raw,y_raw,x_raw)
                        fit_successful = True
                    except:
                        fit_successful = False
                        st.error('Fit not converging')

                    col_fit_1.latex(r'''  
                        A(t) = \frac{k_e * k_a}{cl * (k_a - k_e)}  * (e^{-k_e * t} - e^{-k_a * t})
                        ''')
                    if fit_successful:
                        col_fit_1.latex(       r''' k_a = ''' + f'{ka:.2f}')
                        col_fit_1.latex(       r''' k_e = ''' + f'{ke:.2f}')
                        col_fit_1.latex(       r''' cl = ''' + f'{cl:.2f}')

                        plot_fitresult = True
                        fig_fit.add_trace(go.Scatter(x=x_fit,y=y_fit,fill='tozeroy',mode='none')) #mode = 'markers+lines' or delete 
                        
                        fitresults_tissue['ka'] = round(ka,4)
                        fitresults_tissue['ke'] = round(ke,4)
                        fitresults_tissue['cl'] = round(cl,4)

            elif fitmodel == biexpelim2:
                if len(x_raw) < 3:
                    st.error('At least 3 datapoints are needed for fitting! ')
                    plot_fitresult = False
                    fit_successful = False
                else:
                    try:
                        lambda1, lambda2, A1, y_fit, AUC_calculated,AUC_extrapolated = biexpelim_func2(x_raw,y_raw,x_fit)
                        lambda1, lambda2, A1, y_fit_quality, AUC_calculated,AUC_extrapolated = biexpelim_func2(x_raw,y_raw,x_raw)
                        fit_successful = True
                    except:
                        fit_successful = False
                        st.error('Fit not converging')

                    col_fit_1.latex(r'''  
                        A(t) = A1 * [ e^{-lambda_1 * t} - e^{- lambda_2 * t})
                        ''')
                    if fit_successful:
                        col_fit_1.latex(       r''' lambda_1 = ''' + f'{lambda1:.2f}')
                        col_fit_1.latex(       r''' lambda_2 = ''' + f'{lambda2:.2f}')
                        col_fit_1.latex(       r''' A1 = ''' + f'{A1:.2f}')

                        plot_fitresult = True
                        fig_fit.add_trace(go.Scatter(x=x_fit,y=y_fit,fill='tozeroy',mode='none')) #mode = 'markers+lines' or delete 
                        
                        fitresults_tissue['lambda1'] = round(lambda1,4)
                        fitresults_tissue['lambda2'] = round(lambda2,4)
                        fitresults_tissue['A1'] = round(A1,4)

            elif fitmodel == biexp_abs_elim_two_components:
                if len(x_raw) < 4:
                    st.error('At least 4 datapoints are needed for fitting! ')
                    plot_fitresult = False
                    fit_successful = False
                else:
                    try:
                        lambda1, lambda2, A1, A2, y_fit, AUC_calculated,AUC_extrapolated = biexpelim_func3(x_raw,y_raw,x_fit)
                        lambda1, lambda2, A1, A2, y_fit_quality, AUC_calculated,AUC_extrapolated = biexpelim_func3(x_raw,y_raw,x_raw)
                        fit_successful = True
                    except:
                        fit_successful = False
                        st.error('Fit not converging')

                    col_fit_1.latex(r'''  
                        A(t) = A1 * e^{-lambda_1 * t} - A2 * e^{- lambda_2 * t}
                        ''')
                    if fit_successful:
                        col_fit_1.latex(       r''' lambda_1 = ''' + f'{lambda1:.2f}')
                        col_fit_1.latex(       r''' lambda_2 = ''' + f'{lambda2:.2f}')
                        col_fit_1.latex(       r''' A1 = ''' + f'{A1:.2f}')
                        col_fit_1.latex(       r''' A2 = ''' + f'{A2:.2f}')

                        plot_fitresult = True
                        fig_fit.add_trace(go.Scatter(x=x_fit,y=y_fit,fill='tozeroy',mode='none')) #mode = 'markers+lines' or delete 
                        
                        fitresults_tissue['lambda1'] = round(lambda1,4)
                        fitresults_tissue['lambda2'] = round(lambda2,4)
                        fitresults_tissue['A1'] = round(A1,4)
                        fitresults_tissue['A2'] = round(A2,4)

            else:
                st.error('not defined yet')
                plot_fitresult = False


            fitresults_tissue['AUC calculated'] = round(AUC_calculated,4)
            try:
                fitresults_tissue['AUC extrapolated'] = round(AUC_extrapolated,4)
                AUC_ratio = 100* AUC_extrapolated/AUC_calculated
                if AUC_ratio > 15:
                    st.error(f'Warning! {AUC_ratio:.2f}% of total AUC is extrapolated')
                # else:
                #     st.write(f'{AUC_ratio:.4f}% of total AUC is extrapolated')
                fitresults_tissue['Percentage of extrapolated TIAC [%]'] = round(AUC_ratio,2)
            except:
                fitresults_tissue['Percentage of extrapolated TIAC [%]'] = 0.

            if fitmodel in [monofit, monoexp_dec_fit, bifit, biexpelim, biexpelim2,biexp_abs_elim_two_components]:
                # Fit quality parameters:
                if fit_successful:
                    mse_f,mae_f,rmse_f,r2_f = fit_quality_calc(y_raw,y_fit_quality)  # get mean squared error, mean absolute error, root mean squared error and r squared
                    fitresults_tissue['MSE'] = round(mse_f,4)
                    fitresults_tissue['MAE'] = round(mae_f,4)
                    fitresults_tissue['RMSE'] = round(rmse_f,4)
                    fitresults_tissue['R squared'] = round(r2_f,4)
                    col_fit_3.markdown(f'**R squared: {r2_f:.2f}**')
                else:
                    st.error('data not fitted')
            with col_fit_3.expander('fit results'):
                st.table(fitresults_tissue)
                if fit_successful:
                    st.write(f'MSE: Mean Squared Error')
                    st.write(f'MAE: Mean Absolute Error')
                    st.write(f'RMSE: Root Mean Squared Error')
                    st.write(f'R squared: Coefficient of Determination')

            with col_fit_3:
                yaxis_type = st.radio(
                    f"{tissue} [%ID/g] axis representation",
                    ('normal','logarithmic'),key=f'{tissue}'+'_axiskey_'+str(ik))
            if yaxis_type == 'logarithmic':
                plotLog = True
            else:
                plotLog = False
            if plotLog:
                fig_fit.update_yaxes(type="log")

            if plot_fitresult:
                fig_fit.add_trace(
                    go.Scatter(
                        x=x_fit,
                        y=y_fit,
                        mode="lines",
                        line=go.scatter.Line(color="red"),
                        showlegend=False)
                        )
                fig_fit.add_annotation(
                    xref="x domain",
                    yref="y domain",
                    # The arrow head will be 25% along the x axis, starting from the left
                    x=0.5,
                    # The arrow head will be 40% along the y axis, starting from the bottom
                    y=1.1,
                    text = f'{tissue}: {fitmodel}',
                    # text="AUC from 0 to infinity",
                    # arrowhead=1,
                    showarrow = False,
                    font=dict(
                        family="Courier New, monospace",
                        size=14,
                        color="Black"
                        ),
                        )
                fig_fit.add_annotation(
                    xref="x domain",
                    yref="y domain",
                    # The arrow head will be 25% along the x axis, starting from the left
                    x=0.5,
                    # The arrow head will be 40% along the y axis, starting from the bottom
                    y=1.17,
                    text = f'R squared: {r2_f:.2f}',
                    showarrow = False,
                    font=dict(
                        family="Courier New, monospace",
                        size=14,
                        color="Red"
                        ),
                        )
                
            fig_fit.update_layout(
                yaxis_title=f"[%ID/g]",
                xaxis_title=f"{time_keyword} [{time_unit}]",
            )
            st.plotly_chart(fig_fit)

           
            #Store the png bytes object in a dict
            try:
                all_figures_fit[tissue] = fig_fit.to_image(format='png')
            except:
                pass

            with col_fit_3:
                try:
                    tiac = float(AUC_calculated)/100
                    if fit_successful:
                        tiac_calculated = True
                    else:
                        tiac_calculated = False
                    if tiac <= 0:
                        tiac_calculated = False
                    if np.isnan(tiac):
                        tiac_calculated = False

                except:
                    tiac = 0.
                    tiac_calculated = False
                fitresults_tissue['fitted'] = tiac_calculated

            fitresults_tissue['fitmodel'] = fitmodel
            fitresults_tissue[keyword_mtiac_g] = round(tiac,4)

            fitresults_all[tissue] = fitresults_tissue

            fit_submitted = st.form_submit_button(f'Re-run fit for {tissue}')
                
        # Add a download button to download the figure (traces with fit result)
        if st.button(f'Download Chart {tissue}',key=f'{tissue}-dwnld'):
            # Convert the figure to an image
            try:    
                image = fig_fit.to_image(format='png')
            except:
                pass


            # Prompt the user to enter a new file name
            new_file_name = st.text_input('Enter a new file name', f'{tissue}-fit.png')

            # Download the image with the new file name
            st.download_button(label='Download', data=image, file_name=new_file_name, mime='image/png')

    fitresults = pd.DataFrame.from_dict(fitresults_all)
    if decay_corrected:
        fitresults.insert(loc=0, column = 'Decay corrected',value = radioisotope1)
    else:
        fitresults.insert(loc=0, column = 'Decay corrected',value = 'no decay correction')
    fitresults.insert(loc=0, column = 'parameter',value = fitresults.index)
    fitresults = pd.DataFrame.from_dict(fitresults)

    # Get a list of all organs that could not be fitted
    tiac_not_fitted = [key_organ for key_organ,organ_fitted in fitresults.loc['fitted'].items() if not organ_fitted] 
    tiac_calculated = True
    if len(tiac_not_fitted) > 0:
        tiac_calculated = False
        for not_fitted in tiac_not_fitted:
            st_error_columns = st.columns(2)
            st_error_columns[0].error(f'Fit {not_fitted} was not converging')
            st_error_columns[1].error(f'Please change selected fit model')
            st_error_columns[1].markdown("Go to [%s](#%s)"%(str(not_fitted),str('anchor-fitdecay-')+str(not_fitted)), unsafe_allow_html=True)

    return [tiac_calculated, data_tissues, fitresults, x_fit, all_figures_fit]

def tiac_results_download(fitresults, data_input, rawdata=[], all_figures_fit=[], rayz_id=''):
    st.write('**See results**')
    st.write(fitresults.astype(str))

    output_data_tiac = io.BytesIO()
    with pd.ExcelWriter(output_data_tiac) as writer:  
        fitresults.to_excel(writer, sheet_name = "Decay fit", index = False)
        if len(data_input) != 0:
            data_input.to_excel(writer, sheet_name = "Data input decay corrected", index = False)
        if len(rawdata) != 0:
            rawdata.to_excel(writer, sheet_name = "Data input", index = False)
        workbook  = writer.book
        for fig_nr,figure_fit in enumerate(all_figures_fit):
            img_stream = pd.io.common.BytesIO(all_figures_fit[figure_fit])
            img = workbook.add_worksheet(f'{figure_fit}-{fig_nr}')
            img.insert_image('D2', 'image.png', {'image_data': img_stream})
            img.set_column('A:A', 35)
            img.set_column('B:B', 25)
            fitresults[['parameter',figure_fit]].to_excel(writer, sheet_name = f'{figure_fit}-{fig_nr}', index = False)

        
    st.download_button(label='Download Fit Results',
                            data=output_data_tiac.getvalue(),
                            file_name= f"{rayz_id}-results-decayfit.xlsx")
    

def all_results_download(results_dosimetry_df,results_scaling_df,fitresults, data_input, rawdata, results_filtered, all_figures_fit, rayz_id, olinda_input_df = []):
    output_data = io.BytesIO()

    #Create one dataframe results_total with all dosimetry calculation information (Olinda output and calculator input)
    id_cols = ['Molecule-Batch-ID', 'RAYZ ID', 'Batch ID']
    results_total = results_dosimetry_df[id_cols]
    results_total['tissue Olinda output'] = results_dosimetry_df['Target tissue']
    if 'tumor' in list(results_total['tissue Olinda output']) or 'Tumor' in list(results_total['tissue Olinda output']):
        headers_col = [keyword_dose,'Dose limiting organ',keyword_tumor_dose_mta,'radioisotope 2','Dosimetry method',doselimit_keyword,mta_organ_keyword,mta_keyword]
    else:
        headers_col = [keyword_dose,'Dose limiting organ','radioisotope 2','Dosimetry method',doselimit_keyword,mta_organ_keyword,mta_keyword]
    results_total[headers_col] = results_dosimetry_df[headers_col]

    scaling_df_upload = results_scaling_df
    scaling_df_upload[id_cols] = [results_total['Molecule-Batch-ID'].iloc[0],results_total['RAYZ ID'].iloc[0],results_total['Batch ID'].iloc[0]]
    scaling_df_upload['Dosimetry method'] = results_total['Dosimetry method'].iloc[0]

    if isinstance(fitresults, pd.DataFrame):
        if keyword_projected_bm in list(fitresults.columns):
            fitmodel = fitresults[keyword_projected_bm][fitresults['parameter']=='fitmodel'].values[0]
            tissue_list = list(scaling_df_upload['tissue'])
            com = []
            for tis in tissue_list:
                if tis == keyword_projected_bm:
                    com.append(fitmodel)
                else:
                    com.append(None)
            scaling_df_upload['Comment'] = com
    results_total = pd.concat([results_total,scaling_df_upload])

    # Create a dataframe in a format that can be uploaded to CDD Vault including only the tissues that were submitted to the calculator
    upload_cdd = scaling_df_upload
    upload_cdd_combined = dict()

    for organ in list(upload_cdd['tissue']):
        organ_spellings = get_all_possible_tissue_spellings(organ)
        match_organ_names = list(set(organ_spellings).intersection(list(results_total['tissue Olinda output'])))
        if len(match_organ_names)== 0:
            match_organ_names = ['Tot Body']
        for organ_name in match_organ_names:
            olinda_output = results_total[results_total['tissue Olinda output'] == organ_name]
            olinda_output['tissue'] = organ
            olinda_output['combine'] = 'yes'
            olinda_output.set_index('combine',inplace = True)
            break
        scaling_output = scaling_df_upload[upload_cdd['tissue']==organ]
        scaling_output['tissue Olinda output'] = organ_name
        scaling_output['combine'] = 'yes'
        scaling_output.set_index('combine',inplace = True)
        upload_cdd_combined_organ = scaling_output.combine_first(olinda_output[headers_col])
        upload_cdd_combined_organ_values = upload_cdd_combined_organ.T['yes'].tolist()
        upload_cdd_combined[organ] = upload_cdd_combined_organ_values

    upload_cdd_combined = pd.DataFrame.from_dict(upload_cdd_combined)
    upload_cdd_combined['parameter'] = list(upload_cdd_combined_organ.T.index)
    upload_cdd_combined = upload_cdd_combined.set_index('parameter', drop=True)
    upload_cdd_combined = upload_cdd_combined.T
    
    # This is to re-arrange the order of columns in the upload sheet
    upload_cols = [   "Molecule-Batch-ID",  "RAYZ ID",  "Batch ID",  
        "tissue", "Olinda Input Organ Name", "tissue Olinda output",  "absorbed dose [mSv/MBq = mGy/MBq]",
        "Dose limit [Gy]", mta_organ_keyword, "Dose limiting organ",   "MTA [GBq]",  "tumor exposure at MTA [Gy]",
        "hTIAC_g [MBq * h/MBq]",  "hTIAC_org [MBq * h/MBq]", "mTIAC_org [MBq * h/MBq]", "mTIAC_g [MBq * h/MBq]",
        "mouse organ weight [g]","human organ weight [g]","wb_mouse [g]","wb_human [g]",
        "radioisotope 1 (used in BioD)","radioisotope 2",
        "alpha used for scaling",  "alpha calculated",
        "decay fit method",  "scaling method", "Dosimetry method", "Comment"]
    
    upload_cols = [element for element in upload_cols if element in list(upload_cdd_combined.columns)]
    upload_cdd_combined = upload_cdd_combined[upload_cols]

    st.write(upload_cdd_combined)

    with pd.ExcelWriter(output_data) as writer:  
        upload_cdd_combined.to_excel(writer, sheet_name = "uploadsheetCDD", index = False)
        results_total.to_excel(writer, sheet_name = "All results Olinda", index = False)
        if isinstance(olinda_input_df, pd.DataFrame):
            olinda_input_df = olinda_input_df.reset_index()
            olinda_input_df.to_excel(writer, sheet_name = "Olinda Kinetics Input", index = False)
        results_dosimetry_df.to_excel(writer, sheet_name = "Dosimetry", index = False)
        results_scaling_df.to_excel(writer, sheet_name = "Scaling mTIAC -> hTIAC", index = False)
        if isinstance(fitresults, pd.DataFrame):
            fitresults.to_excel(writer, sheet_name = "Decay fit", index = False)
        data_input.to_excel(writer, sheet_name = "Data input decay corrected", index = False)
        if isinstance(rawdata, pd.DataFrame):
            rawdata.to_excel(writer, sheet_name = "Data input", index = False)
        if isinstance(results_filtered, pd.DataFrame):
            results_filtered.to_excel(writer, sheet_name = "All data from CDD", index = False)

        # To return the fitresults and plot for each organ individually on a separate worksheet
        workbook  = writer.book
        for fig_nr,figure_fit in enumerate(all_figures_fit):
            img_stream = pd.io.common.BytesIO(all_figures_fit[figure_fit])
            img = workbook.add_worksheet(f'{figure_fit}-{fig_nr}')
            img.insert_image('D2', 'image.png', {'image_data': img_stream})
            img.set_column('A:A', 35)
            img.set_column('B:B', 25)
            fitresults[['parameter',figure_fit]].to_excel(writer, sheet_name = f'{figure_fit}-{fig_nr}', index = False)


    st.download_button(label='Download Results',
                            data=output_data.getvalue(),
                            file_name= f'{rayz_id}-dose_calculations.xlsx')

   
