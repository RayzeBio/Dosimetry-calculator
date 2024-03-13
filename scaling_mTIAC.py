import io
import sys
import streamlit as st
import pandas as pd
from fnmatch import fnmatch
from BioDfunctions import *

def _get_state():
    try:
        session = st.session_state()
    except:
        session = sys
    if not hasattr(session, "submitted"):
        session.submitted = False
    if not hasattr(session, "tiac_calculated"):
        session.tiac_calculated = False
    if not hasattr(session, "data_input"):
        st.session_state.data_input = dict()
    if not hasattr(session, 'decay_corr_alpha_scal'):
        st.session_state.decay_corr_alpha_scal = False
    if not hasattr(session, 'calculate_bm'):
        st.session_state.calculate_bm = False
    return session


def app(CDD_TOKEN='None'):
    # Initialization
    radioisotope1 = 'n/a'
    st.session_state = _get_state()

    st.title(irdc_version)
    markdowntext = "[%s](#%s)"%(str('Input'),str('anchor-input'))
    st.sidebar.markdown(markdowntext, unsafe_allow_html=True)
    st.header(f'Input', anchor=str('anchor-input'))

    # Allow user to enter information about molecule 
    molecule_batch_ID = str(st.text_input('Molecule Batch ID from CDD Vault','RB-0000000-000',key='molbatchID-dosimetry-idg')).lstrip()
    rayz_id = st.text_input('RAYZ ID','RAYZ-')
    batch_registration_id = st.text_input('Batch-ID','') 

    tissues = get_selected_tissues_list(tissue_masses_mouse_file)

    fitmodel = f'input is {keyword_mtiac_g}'

    with st.form(key='mTIAC_g'):
        st.write(f'Please enter the TIAC from mouse per g (for each tissue)')
        data_input = {}
        for ie,tissue_ie in enumerate(tissues):
            data_input[tissue_ie] = float(st.text_input(f'{tissue_ie} {keyword_mtiac_g}',value = '1.0',key=keyword_mtiac_g+str(tissue_ie)))
        submitted = st.form_submit_button("Submit")
        if submitted:
            st.session_state.submitted = True
            st.session_state.tiac_calculated=True



    if st.session_state.tiac_calculated:
        data_input = pd.DataFrame.from_dict(data_input, orient = 'index', columns=[keyword_mtiac_g])
        st.session_state.data_input = data_input

        # Create fitresults dataframe for following dosimetry
        fitresults = data_input.T
        fitresults['parameter'] = [keyword_mtiac_g]
        fitresults = fitresults.T
        fit_list = len(data_input.T.columns)* ['mTIAC_g as input']
        fit_list.append('fitmodel')
        fitresults['fitmodel'] = fit_list
        fitresults = fitresults.T
        fitresults['Decay corrected'] = 'n/a'
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

        scaling_method, st.session_state.calc_scaling = select_scaling(st.session_state.calc_scaling,scaling_options=[no_scaling,rel_mass_scal,alpha_scal])
      
        if st.session_state.calc_scaling:
            st.write(f'You selected: {scaling_method}')
            results_scaling_df = scaling_mTIAC_hTIAC(scaling_method,fitresults,radioisotope1,x_fit=[],data_tissues=[],cdd_weights=[],blood_key=blood_key)

            ############################################################################
            ##### Next section: Dosimetry from human residence time
            ############################################################################
            markdowntext = "[%s](#%s)"%(str('Dosimetry'),str('anchor-dosimetry-')+str(rayz_id))
            st.sidebar.markdown(markdowntext, unsafe_allow_html=True)
            st.header(f'Dosimetry', anchor=str('anchor-dosimetry-')+str(rayz_id))          

            with st.expander('Select radioisotope for absorbed dose calculation'):
                radioisotope2 = st.radio(f"Radioisotope 2: projected use in human study",(isotopes_human_halflives.keys()))                                
            st.write(f'You selected {radioisotope2}')

            try:
                df_sfactors = pd.read_csv(f"ICRP89_DF_{radioisotope2}.csv")
            except:
                st.error(f'Dose factors for {radioisotope2} are not defined, please ask admin for future implementation or select different isotope for dosimetry')
            # results_total, results_dosimetry_df = dosimetry_from_hTIAC_org(rayz_id, batch_registration_id, results_scaling_df,molecule_batch_ID,key_human_wb_weight,key_human_organ_weight,df_sfactors,radioisotope2,doselimits_file)
            olinda_input_df, results_doselimits_df = dosimetry_from_hTIAC_org(rayz_id, batch_registration_id, results_scaling_df,molecule_batch_ID,df_sfactors,radioisotope2,doselimits_file)

            #############################################################################
            # Dosimetry finished, download results button
            #############################################################################
            markdowntext = "[%s](#%s)"%(str('Results'),str('anchor-')+str(rayz_id))
            st.sidebar.markdown(markdowntext, unsafe_allow_html=True)
            st.header(f'Congratulations: Dosimetry is done!', anchor=str('anchor-')+str(rayz_id))

            st.write(f'Dosimetry was calculated using **{scaling_method}** as scaling method')
            st.write(f'Download your results here (and upload them to CDD Vault!!! )')

            all_results_download(results_dosimetry_df=results_doselimits_df,results_scaling_df=results_scaling_df,fitresults=fitresults, data_input=data_input, rawdata=[],results_filtered=[],all_figures_fit=[], rayz_id=rayz_id, olinda_input_df=olinda_input_df)

    # if st.session_state.tiac_calculated:
    #     st.title('Human extrapolation')
    #     if st.checkbox('Show comparison of extrapolation methods'):
    #         st.image(Image.open('./pics/comparison-scaling-methods.png'), caption='Comparison with rawdata from Beykan et al. - yellow indicates real measured hTIAC')  
    #     scaling_method = st.radio(
    #             f"Scaling method for mTIAC -> hTIAC",
    #             ([no_scaling,rel_mass_scal,alpha_scal]))

      
    #     if True:
    #         st.write(f'You selected: {scaling_method}')
    #         if scaling_method == alpha_scal:
    #             st.session_state.decay_corr_alpha_scal = st.checkbox(f'Calculate alpha with decay correction? (Check if Yes)',value=True,key="decay_correction_alpha_scaling"+tissue_ie)

    #         wb_entered, results_scaling_df, radioisotope1, radioisotope2 = scaling_mTIAC_hTIAC(scaling_method,fitmodel,tissues,st.session_state.data_input[keyword_mtiac_g],keyword_mtiac_g,radioisotope1,decay_corr_alpha_scal = st.session_state.decay_corr_alpha_scal)
          
    #         st.title('Dosimetry')

    #         if not scaling_method == alpha_scal:
    #             radioisotope2 = st.radio(
    #                 f"Radioisotope 2: projected use in human study",
    #                 (isotopes_human_halflives.keys()))                                
    #             st.write(f'You selected {radioisotope2}')

    #         df_sfactors = pd.read_csv(f"ICRP89_DF_{radioisotope2}.csv")

    #         results_dosimetry_df = dosimetry_from_hTIAC_g(rayz_id,batch_registration_id,results_scaling_df,molecule_batch_ID,key_human_wb_weight,key_human_organ_weight,df_sfactors,radioisotope2,doselimits_file)

    #         st.title('Congratulations: Dosimetry is done!')
    #         st.write(f'Dosimetry was calculated using **{fitmodel}** fit method and **{scaling_method}** as scaling method')
    #         st.write(f'Download your results here (and upload them to CDD Vault!!! )')

    #         try:
    #             if htiac_org_key in results_dosimetry_df.keys():
    #                 results_dosimetry_df=results_dosimetry_df.drop(columns = [htiac_org_key])
    #             if htiac_key in results_dosimetry_df.keys():
    #                 results_dosimetry_df=results_dosimetry_df.drop(columns = [htiac_key])


    #             results_total = results_dosimetry_df.join(results_scaling_df)

    #             output_data = io.BytesIO()
    #             with pd.ExcelWriter(output_data) as writer:  
    #                 results_total.to_excel(writer, sheet_name = "uploadsheetCDD", index = False)
    #                 results_dosimetry_df.to_excel(writer, sheet_name = "Dosimetry", index = False)
    #                 results_scaling_df.to_excel(writer, sheet_name = "Scaling mTIAC -> hTIAC", index = False)
    #                 # fitresults.to_excel(writer, sheet_name = "Decay fit", index = False)
    #                 st.session_state.data_input.to_excel(writer, sheet_name = "Data input", index = False)
    #                 writer.save()
                    
    #             st.download_button(label='Download Results',
    #                                     data=output_data.getvalue(),
    #                                     file_name= 'dose_calculations.xlsx')
    #         except:
    #             pass

    st.write('----')
    if st.button('Clear all calculations: Click button and reload homepage ( "Ctrl + R" or F5 )'):
        st.experimental_memo.clear()
        st.experimental_singleton.clear()
        st.stop()


if __name__ == "__main__":
    app()