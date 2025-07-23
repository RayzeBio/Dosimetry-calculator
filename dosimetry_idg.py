import streamlit as st
import pandas as pd
from fnmatch import fnmatch
from BioDfunctions import *

def _get_state(session):
    if not hasattr(session, "submitted"):
        session.submitted = False
    if not hasattr(session, "dragdropinput"):
        session.dragdropinput = False        
    if not hasattr(session, "calc_scaling"):
        session.calc_scaling = False
    if not hasattr(session, "decay_correction"):
        session.decay_correction = False
    if not hasattr(session, 'decay_corr_alpha_scal'):
        session.decay_corr_alpha_scal = False
    return session


def app(CDD_TOKEN='None'):
    ########################################################
    # Initialization to create a memory for streamlit app
    ########################################################
    st.session_state = _get_state(st.session_state)
    if not hasattr(st.session_state, "dragdropinput"):
        st.session_state.dragdropinput = False    
    tiac_calculated = False 

    st.title(irdc_version)
    markdowntext = "[%s](#%s)"%(str('Input'),str('anchor-input'))
    st.sidebar.markdown(markdowntext, unsafe_allow_html=True)
    st.header(f'Input', anchor=str('anchor-input'))

    # Allow user to enter information about molecule 
    molecule_batch_ID = str(st.text_input('Molecule Batch ID from CDD Vault','RB-0000000-000',key='molbatchID-dosimetry-idg')).lstrip()
    rayz_id = st.text_input('RAYZ ID','RAYZ-')
    batch_registration_id = st.text_input('Batch-ID','') 

    # Data input from file (csv)
    data_input_dragdrop = st.file_uploader("Select a dataset with %ID/g from .csv file with columns 'Time [hr]', 'Organ names' ", type = ["csv"],on_change=reset_session_state_dragndropInput)
    if data_input_dragdrop is not None:
        df_input = pd.read_csv(data_input_dragdrop)
        data_input,tissues,keywords_list,st.session_state.dragdropinput = biod_input_file(df_input)

    # Manual data input 
    if not st.session_state.dragdropinput:
        data_input,tissues,keywords_list,st.session_state.submitted = biod_input_manually(st.session_state.submitted)

    try:
        data_input = pd.DataFrame.from_dict(data_input)
    except:
        st.session_state.submitted = False
        st.session_state.dragdropinput = False


    if st.session_state.submitted or st.session_state.dragdropinput:
        # sort input data by time
        data_input.sort_values([time_keyword], inplace=True)
        # make %ID/g input numeric
        for datakey in data_input.keys():
            data_input[datakey]= pd.to_numeric(data_input[datakey], errors='coerce')

        # checkbox option if input data should be decay corrected
        st.session_state.decay_correction = select_decay_correction_input()

        # select isotope for decay correction of input
        radioisotope1 = select_radioisotope1(rayz_id)
                        
        # decay correction of input data
        rawdata, data_input = decay_corr_input(data_input,radioisotope1)


    if st.session_state.submitted or st.session_state.dragdropinput:
        sex_key= st.selectbox('Please select phantom gender',['male','female'])

        with st.expander('show input'):
            st.write(data_input)

        # plot rawdata
        plot_rawdata_bioD(keywords_list,tissues,data_input)


        ###########################################################################
            # Decay FIT
        ###########################################################################
        markdowntext = "[%s](#%s)"%(str('Fit'),str('anchor-fit-')+str(rayz_id))
        st.sidebar.markdown(markdowntext, unsafe_allow_html=True)
        st.header(f'Decay fit', anchor=str('anchor-fit-')+str(rayz_id))

        tiac_calculated, data_tissues, fitresults, x_fit, all_figures_fit = fit_decay_fitmodel(data_input,keywords_list,time_keyword,injdose_keyword=injdose_keyword, radioisotope1=radioisotope1, decay_corrected=st.session_state.decay_correction)

    if tiac_calculated:
        # Download button for decay fit results and plots
        tiac_results_download(fitresults, data_input, rawdata, all_figures_fit, rayz_id)


        ############################################################################
        ##### Next section: Bone Marrow Estimation from blood residence time
        ############################################################################

        markdowntext = "[%s](#%s)"%(str('BoneMarrowEstimation'),str('anchor-bonemarrow-')+str(rayz_id))
        st.sidebar.markdown(markdowntext, unsafe_allow_html=True)
        st.header(f'Bone Marrow Estimation', anchor=str('anchor-bonemarrow-')+str(rayz_id))

        st.session_state.calculate_bm = st.checkbox('Do you want to project Bone Marrow uptake from blood?', value=False)

        if st.session_state.calculate_bm:
            bm_method,blood_key = select_bonemarrow_projection(fitresults)
            fitresults = project_bonemarrow(fitresults,blood_key,bm_method)
        else:
            blood_key = None

        ############################################################################
        ##### Next section: Extrapolation from animal to human residence time
        ############################################################################

        markdowntext = "[%s](#%s)"%(str('Extrapolation'),str('anchor-extrapolation-')+str(rayz_id))
        st.sidebar.markdown(markdowntext, unsafe_allow_html=True)
        st.header(f'Human extrapolation: Mouse to Human', anchor=str('anchor-extrapolation-')+str(rayz_id))

        scaling_method, st.session_state.calc_scaling = select_scaling(st.session_state.calc_scaling)
      
        if st.session_state.calc_scaling:
            st.write(f'You selected: {scaling_method}')

            results_scaling_df = scaling_mTIAC_hTIAC(scaling_method,fitresults,radioisotope1,x_fit=x_fit,data_tissues=data_tissues,cdd_weights=[],blood_key=blood_key)

            ############################################################################
            ##### Next section: Dosimetry from human residence time
            ############################################################################
            markdowntext = "[%s](#%s)"%(str('Dosimetry'),str('anchor-dosimetry-')+str(rayz_id))
            st.sidebar.markdown(markdowntext, unsafe_allow_html=True)
            st.header(f'Dosimetry', anchor=str('anchor-dosimetry-')+str(rayz_id))          

            with st.expander('Select radioisotope for absorbed dose calculation'):
                radioisotope2 = st.radio(
                    f"Radioisotope 2: projected use in human study",
                    options=(isotopes_BioD_halflives.keys()), index=list(isotopes_BioD_halflives.keys()).index(radioisotope1))                                
            st.write(f'You selected {radioisotope2}')

            # results_total, results_dosimetry_df = dosimetry_from_hTIAC_org(rayz_id, batch_registration_id, results_scaling_df,molecule_batch_ID,key_human_wb_weight,key_human_organ_weight,df_sfactors,radioisotope2,doselimits_file)
            olinda_input_df, results_doselimits_df = dosimetry_from_hTIAC_org(rayz_id, batch_registration_id, results_scaling_df,molecule_batch_ID,sex_key,radioisotope2,doselimits_file)


            #############################################################################
            # Dosimetry finished, download results button
            #############################################################################
            markdowntext = "[%s](#%s)"%(str('Results'),str('anchor-')+str(rayz_id))
            st.sidebar.markdown(markdowntext, unsafe_allow_html=True)
            st.header(f'Congratulations: Dosimetry is done!', anchor=str('anchor-')+str(rayz_id))

            st.write(f'Dosimetry was calculated using **{scaling_method}** as scaling method')
            st.write(f'Download your results here (and upload them to CDD Vault!!! )')

            all_results_download(results_doselimits_df,results_scaling_df,fitresults, data_input, rawdata,[],all_figures_fit, rayz_id,olinda_input_df)


    st.write('----')
    if st.button('Clear all calculations: Click button and reload homepage ( "Ctrl + R" or F5 )'):
        st.cache_data.clear()
        st.cache_resource.clear()
        st.session_state = _get_state(st.session_state)
        st.session_state.submitted = False
        st.session_state.dragdropinput = False
        st.session_state.calc_scaling = False
        st.session_state.decay_correction = False
        st.stop()


if __name__ == "__main__":
    app()