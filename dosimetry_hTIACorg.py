import io
import sys
import streamlit as st
import pandas as pd

from scipy.optimize import curve_fit
from scipy.stats import linregress
from fnmatch import fnmatch
from BioDfunctions import *

def _get_state():
    try:
        session = st.session_state()
    except:
        session = sys
    if not hasattr(session, "submitted"):
        session.submitted = False
    if not hasattr(session, "calc_scaling"):
        session.calc_scaling = False
    if not hasattr(session, "tiac_calculated"):
        session.tiac_calculated = False
    return session


def app(CDD_TOKEN='None'):
    ########################################################
    # Initialization to create a memory for streamlit app
    ########################################################
    st.session_state = _get_state()

    st.title(irdc_version)
    markdowntext = "[%s](#%s)"%(str('Input'),str('anchor-input'))
    st.sidebar.markdown(markdowntext, unsafe_allow_html=True)
    st.header(f'Input', anchor=str('anchor-input'))
    
    # Allow user to enter information about molecule 
    molecule_batch_ID = str(st.text_input('Molecule Batch ID from CDD Vault','RB-0000000-000',key='molbatchID-dosimetry-idg')).lstrip()
    rayz_id = st.text_input('RAYZ ID','RAYZ-')
    batch_registration_id = st.text_input('Batch-ID','') 

    # get all available tissues 
    tissues = get_selected_tissues_list(tissue_masses_mouse_file)

    # Allow user to enter hTIAC_org
    with st.form(key='hTIAC_org'):
        st.write('Enter human residence time per organ [MBq*h/MBq] (= input for Olinda)')
        data_input = {}
        for ie,tissue_ie in enumerate(tissues):
            data_input[tissue_ie] = list()
        for j, (k, v) in enumerate(data_input.items()):
            v.append((st.text_input(f'{k} {htiac_org_key}',key=str(j)+str(k))))
        submitted = st.form_submit_button("Submit")
        if submitted:
            st.session_state.submitted = True
            st.session_state.tiac_calculated=True
    

    ########################################################################################################################################################
    ##### Next section: Extrapolation from animal to human residence time is skipped, input is converted to align to other dosimetry calculations
    ########################################################################################################################################################
            
    if st.session_state.tiac_calculated:   
        sex_key= st.selectbox('Please select phantom gender',['male','female'])


        # Allow user to change body weight and organ weight for mouse and human  
        results_scaling_alltissues = dict()
        col_calc_1,col_calc_2 = st.columns(2)
        for ik, tissue in enumerate(tissues):
            try:
                htiac_org = float(data_input[tissue][0])
            except:
                st.error(f'Please enter hTIAC for {tissue} (hTIAC per organ)')
                htiac_org = 0.
            results_scaling = dict()
            if ik > 0:
                wb_entered = True   
            else:
                wb_entered = False             

            if not wb_entered:  #to enter whole body weight only once
                wb_mouse = float(col_calc_1.number_input('body weight mouse [g]', min_value = 0., max_value = 100., value = float(wholebody_mouse), step = 1.))
                wb_human = float(col_calc_2.number_input('body weight human [kg]', min_value = 0., max_value = 200., value = float(wholebody_human), step = 1.))*1000

            org_weight_mouse_from_file = get_mass_from_file(tissue, tissue_masses_mouse_file, mouse_keyword_tissues, mouse_keyword_masses, 1e-3)
            org_weight_mouse = float(col_calc_1.number_input(f'{tissue} weight mouse [g]', min_value = 0., max_value = wb_mouse, value = org_weight_mouse_from_file, step = 0.1))
            org_weight_human_from_file = get_mass_from_file(tissue, tissue_masses_human_file, human_keyword_tissues, human_keyword_masses, 1.)
            org_weight_human = float(col_calc_2.number_input(f'{tissue} weight human [g]', min_value = 0., max_value = wb_human, value = org_weight_human_from_file, step = 100.))


            # Create results_scaling dictionaries to stay consistent with other scripts (dosimetry_CDD and dosimetry_idg) and re-use dosimetry functions later

            htiac = float(htiac_org) / org_weight_human
            results_scaling[htiac_key] = float(htiac)
            results_scaling[htiac_org_key] = float(htiac_org)
            results_scaling['mouse organ weight [g]'] = float(org_weight_mouse)
            results_scaling[key_human_organ_weight] = float(org_weight_human)
            results_scaling['wb_mouse [g]'] = float(wb_mouse)
            results_scaling[key_human_wb_weight] = float(wb_human)
            results_scaling['tissue'] = tissue

            results_scaling_alltissues[tissue] = results_scaling

        results_scaling_df = pd.DataFrame.from_dict(results_scaling_alltissues, orient='index')
        with st.expander('See results:'):
            st.table(results_scaling_df)
        data_input_df = pd.DataFrame.from_dict(data_input)
        tiac_results_download(fitresults=results_scaling_df, data_input=data_input_df, all_figures_fit=[], rayz_id=rayz_id)

        if st.button('Continue with Dosimetry'):
            st.session_state.calc_scaling = True

        if st.session_state.calc_scaling:
            ############################################################################
            ##### Next section: Dosimetry from human residence time
            ############################################################################
            markdowntext = "[%s](#%s)"%(str('Dosimetry'),str('anchor-dosimetry-')+str(rayz_id))
            st.sidebar.markdown(markdowntext, unsafe_allow_html=True)
            st.header(f'Dosimetry', anchor=str('anchor-dosimetry-')+str(rayz_id))    

            # if not scaling_method == alpha_scal:
            radioisotope2 = st.radio(
                f"Radioisotope 2: projected use in human study",
                options=(isotopes_BioD_halflives.keys()), index=0)                                
            st.write(f'You selected {radioisotope2}')
            

            olinda_input_df, results_doselimits_df = dosimetry_from_hTIAC_org(rayz_id, batch_registration_id, results_scaling_df,molecule_batch_ID,sex_key,radioisotope2,doselimits_file)


            #############################################################################
            # Dosimetry finished, download results button
            #############################################################################
            markdowntext = "[%s](#%s)"%(str('Results'),str('anchor-')+str(rayz_id))
            st.sidebar.markdown(markdowntext, unsafe_allow_html=True)
            st.header(f'Congratulations: Dosimetry is done!', anchor=str('anchor-')+str(rayz_id))

            st.write(f'Dosimetry was calculated from hTIAC_org')
            st.write(f'Download your results here (and upload them to CDD Vault!!! )')

            all_results_download(results_dosimetry_df=results_doselimits_df,results_scaling_df=results_scaling_df,fitresults=None, data_input=data_input_df, rawdata=None, results_filtered=[], all_figures_fit=[], rayz_id=rayz_id, olinda_input_df=olinda_input_df)




    st.write('----')
    if st.button('Clear all calculations: Click button and reload homepage ( "Ctrl + R" or F5 )'):
        st.cache_data.clear()
        st.cache_resource.clear()
        st.session_state = _get_state()
        st.session_state.submitted = False
        st.session_state.calc_scaling = False
        st.session_state.tiac_calculated = False
        st.stop()


if __name__ == "__main__":
    app()