import sys

import streamlit as st
import pandas as pd
from BioDfunctions import *
from glob import glob


def append_attribute(session, attribute_name, attribute_state):
    if not hasattr(session, attribute_name):
        session.attribute_name = attribute_state    
    return session

@st.cache_resource
def _get_state():
    try:
        session = st.session_state()
    except:
        session = sys
    session = append_attribute(session, 'searchBioD', False)
    session = append_attribute(session, 'tissues_updated', False)
    session = append_attribute(session, 'token', '')
    session = append_attribute(session, 'token_entered', False)
    session = append_attribute(session, 'BioD_bigDataSet', False)
    session = append_attribute(session, 'calc_scaling', False)
    session = append_attribute(session, 'continueDosimetry', False)
    session = append_attribute(session, 'continue_with_dosimetry', False)
    session = append_attribute(session, 'first_molecule_requested', True)
    session = append_attribute(session, 'process_id', None)
    session = append_attribute(session, 'continueCDDsearch', False)
    session = append_attribute(session, 'decay_correction', False)
    session = append_attribute(session, 'start_tissues', ["tumor", "kidney"])
    session = append_attribute(session, 'moleculefound', False)
    return session


def app(CDD_TOKEN='None'):
    ########################################################
    # Initialization to create a memory for streamlit app
    ########################################################
    st.session_state = _get_state()
    if not hasattr(st.session_state, 'searchBioD'):
        st.session_state.searchBioD = False
    if not hasattr(st.session_state, 'moleculefound'):
        st.session_state.moleculefound = False
    if not hasattr(st.session_state, 'continueCDDsearch'):
        st.session_state.continueCDDsearch = False
    if not hasattr(st.session_state, 'start_tissues'):
        st.session_state.start_tissues=["tumor", "kidney"]
    if not hasattr(st.session_state, 'continue_with_dosimetry'):
        st.session_state.continue_with_dosimetry= False
    if not hasattr(st.session_state, 'calculate_bm'):
        st.session_state.calculate_bm = False

    st.title(f"""**{irdc_version}**
                Download BioD data from CDD Vault
                """)
    st.markdown('-'*17)

    ########################################################
    # Check for proper input of CDD token
    ########################################################
    if CDD_TOKEN == 'None': # If None is entered a text input field to ask for CDD token will be shown, otherwise this section can be skipped (CDD Vault does check if the token is valid)
        token = st.text_input('enter your personal token for CDD Vault', CDD_TOKEN)
        if len(token) > 10:
            st.session_state.token_entered = True
            st.session_state.token = token
        else:
            st.session_state.token_entered = False
            st.session_state.continueCDDsearch = False 
    else:
        token = CDD_TOKEN
        st.session_state.token = CDD_TOKEN
        st.session_state.token_entered = True
        

    ########################################################
        # Start the app after the token was entered
    ########################################################
    if st.session_state.token_entered:
        # If moleculefound from session state was set to true the search in CDD Vault can be started
        if st.session_state.moleculefound:
            st.session_state.continueCDDsearch = True 
        else:
            st.session_state.continueCDDsearch = False 
        token = st.session_state.token


        #############################################################
        # Search for Molecule Batch ID in CDD Vault
        #############################################################
        col1, col2 = st.columns([1,1])
        with col1:
            markdowntext = "[%s](#%s)"%(str('Search'),str('anchor-search'))
            st.sidebar.markdown(markdowntext, unsafe_allow_html=True)
            st.header(f'CDD search', anchor=str('anchor-search'))
            molecule_batch_ID,st.session_state.continueCDDsearch,st.session_state.moleculefound = search_moleculebatchID_CDD(st.session_state.continueCDDsearch,st.session_state.moleculefound)
                
        if st.session_state.continueCDDsearch:
            st.write('You chose: ', molecule_batch_ID)
            if st.session_state.first_molecule_requested: #only run the CDD Vault query if the molecule was entered for the first time, otherwise get information from previous API request through export ID
                st.session_state.export_id_molecule_request, molecule_request,st.session_state.process_id = find_molecule_batch_id(molecule_batch_ID,token,bigDataSet=True) 
                st.session_state.first_molecule_requested = False
            else:
                molecule_request = get_export(st.session_state.export_id_molecule_request, token)

            ######################################################################################
            # Check that molecule details are correct, show RAYZ ID and synonyms in expander
            ######################################################################################   
            with st.expander("Molecule details"):
                try:
                    nr_molecules = molecule_request['count'] 
                    if nr_molecules == 0:
                        st.session_state.moleculefound = False
                except:
                    st.success("Please 'Search in CDD Vault'")
                    nr_molecules = 0
                    st.session_state.moleculefound = False
                    molecule_batch_ID = '000'
                    rayz_id = '000'

                if nr_molecules == 0:
                    st.error('This Molecule-Batch ID does not exist in our Vault')
                    st.session_state.moleculefound = False
                else:
                    st.session_state.moleculefound = True
                    st.write(f'There is {nr_molecules} entry in CDD Vault')
                    col1, col2 = st.columns([1,1])
                    for element in molecule_request['objects']:
                            batch_id = element['id']
                            molecule_id = element['molecule']['id']
                    batch_information = get_batch_information(batch_id,token)
                    rayz_id = batch_information['molecule_fields']['RAYZ Internal ID']
                    batch_registration_id = element['batch_fields']['Batch ID']

                    df_batch_information = pd.DataFrame({
                        'name': pd.Series(batch_information['name']),
                        'RAYZ Internal ID': pd.Series(rayz_id),
                        'synonyms': pd.Series(batch_information['synonyms'])
                        })

                    st.table(data = df_batch_information)

            #####################################################################################
            # Button to confirm to get BioD data uploaded to CDD Vault
            #####################################################################################
            if st.session_state.moleculefound:
                if st.button(f'Get BioD data for {molecule_batch_ID} / {rayz_id}'):
                    st.session_state.searchBioD = True

            if st.session_state.searchBioD and st.session_state.moleculefound:
                BioD_protocol = BioD_protocol_request(st.session_state.token)

                st.title('Uploaded Conditions')
                readout_definitions = BioD_protocol['readout_definitions']

                # only query CDD Vault once for BioD data, get uploaded data from export ID at second run
                if not st.session_state.BioD_bigDataSet:
                    export_id, bioD_data = get_molecules_from_BioD_bigDataset(molecule_id,batch_id,token)
                    st.session_state.export_id = export_id
                    st.session_state.BioD_bigDataSet = True
                else:
                    identified_bioD = get_export(st.session_state.export_id, token)
                    bioD_data = list(filter(lambda item: item['molecule'] == molecule_id and item['batch'] == batch_id, identified_bioD['objects']))

                if len(bioD_data) == 0:
                        st.error(f'There are no BioD data uploaded in CDD Vault for {molecule_batch_ID}')
                        st.session_state.tissues_updated = False
                        st.session_state.searchBioD = False
                else:    
                    # create a dataframe from all readouts uploaded                
                    readouts_list = []
                    for readout_nr in enumerate(bioD_data):
                        readouts_list.append(readout_nr[1]['readouts'])   
                    bioD_df = pd.DataFrame.from_dict(readouts_list)

                    # split readouts into value/note/outlier/outlier-type (CDD logic), exclude outliers from bioD_df
                    outlier_list = []
                    with st.expander('Data uploaded to CDD Vault'):
                        for element in readout_definitions:
                            st.write('%10i: %s'%(element['id'], element['name']))
                        for column_name in bioD_df:
                            if type(bioD_df.loc[:,column_name].values.tolist()[0]) == dict:
                                bioD_df[[f'{column_name}-value', f'{column_name}-note', f'{column_name}-outlier', f'{column_name}-outlier_type']] = bioD_df[column_name].apply(lambda x: pd.Series(extract_values_dict(x),dtype='object'))
                                outlier_list.append(f'{column_name}-outlier')                            
                        for column_name_outlier in outlier_list:
                            bioD_df = bioD_df[bioD_df[column_name_outlier] != True] #only include rows where all outlier fields are marked as False (= no outlier)
                        st.write(bioD_df)

                    # pre-defined standard conditions in BioD protocols that can be individually selected by user to include in dosimetry calculations
                    parameters_experiment_BioD = [  'Data Source',
                                                    'Tumor Type',
                                                    'Species',
                                                    'Analysis Method',
                                                    'Time point']


                    ###### Show a table where specific parameters can be selected
                    ###### UPLOADED CONDITIONS
                    ################################################################
                    selected_conditions_all = select_conditions_BioD(parameters_experiment_BioD,readout_definitions,bioD_data)
               
                    st.title('Uploaded tissues')
                    tissue_cdd_id = get_readout_name_id('Tissue',readout_definitions)
                    all_tissues = get_parameters_fromBioD(tissue_cdd_id,bioD_data)

                    ###### To select all tissues or clear selection
                    button_columns = st.columns(2)
                    if button_columns[0].button('Select all tissues'):
                        st.session_state.start_tissues = all_tissues
                    if button_columns[1].button('Clear selection'):
                        st.session_state.start_tissues = []
                        st.session_state.tissues_updated = False
                        st.session_state.continue_with_dosimetry = False
                    ####### To select individual tissues
                    st.session_state.tissues_updated,tissues = change_tissue_input(all_tissues,st.session_state.start_tissues)

                if st.session_state.tissues_updated:
                    ###### Filter the downloaded results for specific previously selected conditions
                    tissue_list_selected = {}
                    for selected_tissue in tissues:
                        tissue_list_selected[selected_tissue] = True
                    selected_conditions_all['Tissue'] = tissue_list_selected

                    # Filter BioD output only for selected conditions
                    bioD_data_filtered = filter_BioD_conditions(bioD_df,selected_conditions_all,readout_definitions)
                        
                    timepoints = [k for k, v in selected_conditions_all['Time point'].items() if v]
                    timepoints.sort()

                    if len(bioD_data_filtered) == 0:
                        st.error(f'There are no BioD data uploaded in CDD Vault for {molecule_batch_ID}')
                        st.session_state.tissues_updated = False
                        st.session_state.searchBioD = False
                    else:
                        st.write(f'There are BioD data uploaded in CDD Vault for {molecule_batch_ID}')
                        tissues,results_for_dosimetry_Calc_df,results_filtered,condition = get_condition_raw_aver(tissues,timepoints,bioD_data_filtered,readout_definitions)

                        averaged_weight_all_timepoints_df = get_masses_CDD(results_filtered)
                        
                        with st.expander(f'{condition} filtered and averaged'):
                            st.write(results_for_dosimetry_Calc_df)

                        ### Create plot of input data
                        plot_input_data_biod(results_for_dosimetry_Calc_df,condition,rayz_id)

                        # checkbox option if input data should be decay corrected
                        st.session_state.decay_correction = select_decay_correction_input()

                        # select isotope for decay correction of input
                        radioisotope1 = select_radioisotope1(rayz_id)

                        # decay correction of input data
                        rawdata, results_for_dosimetry_Calc_df = decay_corr_input(results_for_dosimetry_Calc_df,radioisotope1)

                        if st.session_state.decay_correction:
                            with st.expander('Decay corrected rawdata'):
                                st.write('**From CDD Vault: Original data**')
                                st.table(rawdata)
                                st.write(f'**Decay corrected for {radioisotope1}**')
                                st.table(results_for_dosimetry_Calc_df)


                        #########################################################################
                        # Download button for CDD rawdata (incl. decay correction if selected)
                        #########################################################################
                        download_cdd_rawdata(rayz_id,radioisotope1,results_for_dosimetry_Calc_df,rawdata,averaged_weight_all_timepoints_df,tissues,results_filtered)

                    if st.button('Continue with Dosimetry'):
                        st.session_state.continue_with_dosimetry = True

            ####################################################################
            # Show BioD data input and change order of submitted organs
            ####################################################################
            if st.session_state.continue_with_dosimetry:
                st.title("BioD data Input")
                try:
                    data_input = results_for_dosimetry_Calc_df
                    st.session_state.continueDosimetry = True
                except:
                    st.error('Please press button')

                if st.session_state.continueDosimetry:                    
                    with st.expander('show input'):
                        st.write(data_input)      
                  
                    keywords_list = []
                    for ie,tissue_ie in enumerate(tissues):
                        keyword = injdose_keyword+tissue_ie
                        keywords_list.append(keyword)

                    plot_rawdata_bioD(keywords_list,tissues,data_input)


                    ###########################################################################
                        # Decay FIT
                    ###########################################################################
                    markdowntext = "[%s](#%s)"%(str('Fit'),str('anchor-fit-')+str(rayz_id))
                    st.sidebar.markdown(markdowntext, unsafe_allow_html=True)
                    st.header(f'Decay fit', anchor=str('anchor-fit-')+str(rayz_id))

                    tiac_calculated, data_tissues, fitresults, x_fit,all_figures_fit = fit_decay_fitmodel(data_input,tissues,time_keyword=time_keyword,injdose_keyword=injdose_keyword, radioisotope1= radioisotope1, decay_corrected=st.session_state.decay_correction)
              
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
                        st.header(f'Human extrapolation', anchor=str('anchor-extrapolation-')+str(rayz_id))

                        scaling_method, st.session_state.calc_scaling = select_scaling(st.session_state.calc_scaling)
                    
                        if st.session_state.calc_scaling:
                            st.write(f'You selected: {scaling_method}')
                            results_scaling_df = scaling_mTIAC_hTIAC(scaling_method,fitresults,radioisotope1,x_fit=x_fit,data_tissues=data_tissues,cdd_weights=averaged_weight_all_timepoints_df,blood_key=blood_key)

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

                            try:
                                df_sfactors = pd.read_csv(f"ICRP89_DF_{radioisotope2}.csv")
                            except:
                                st.error(f'Dose factors for {radioisotope2} are not defined, please ask admin for future implementation or select different isotope for dosimetry')
                            olinda_input_df, results_doselimits_df = dosimetry_from_hTIAC_org(rayz_id, batch_registration_id, results_scaling_df,molecule_batch_ID,df_sfactors,radioisotope2,doselimits_file)


                            #############################################################################
                            # Dosimetry finished, download results button
                            #############################################################################
                            markdowntext = "[%s](#%s)"%(str('Results'),str('anchor-')+str(rayz_id))
                            st.sidebar.markdown(markdowntext, unsafe_allow_html=True)
                            st.header(f'Congratulations: Dosimetry is done!', anchor=str('anchor-')+str(rayz_id))

                            st.write(f'Dosimetry was calculated using **{scaling_method}** as scaling method')
                            st.write(f'Download your results here (and upload them to CDD Vault!!! )')

                            all_results_download(results_doselimits_df,results_scaling_df,fitresults, data_input, rawdata,results_filtered,all_figures_fit, rayz_id, olinda_input_df=olinda_input_df)




    
if __name__ == '__main__':
    st.set_page_config(layout="wide")

    app(CDD_TOKEN='None')