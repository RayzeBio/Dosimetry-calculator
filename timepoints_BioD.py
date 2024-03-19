import numpy as np
import streamlit as st
import matplotlib.pyplot as plt
import pandas as pd

from PIL import Image
from BioDfunctions import *



def app(CDD_TOKEN = 'None'):
    st.header("Best sampling timepoints")

    #Select isotope for planned study:
    radioisotope1 = st.radio(f"Which radioisotope should be used in your study?",(isotopes_human_halflives.keys()))     

    half_life_isotope = float(isotopes_BioD_halflives[radioisotope1])
    lambda_radioisotope1 = np.log(2)/half_life_isotope
    st.write(f'You selected {radioisotope1}: Thalf = {half_life_isotope:.2f} h')

    half_life_biol = st.number_input('Enter the biological half-life of your compound in hours. Recommendation: Use halflife of tumor or the longest halflife that was observed across all tissues',value=1., step=0.01)

    st.write(f'Biological half life of your compound: {half_life_biol:.2f} hours')

    with st.expander('Definition of Effective half life'):
        try:
            st.image(Image.open('./pics/T_Effective.PNG'), caption='Image 1: Definition of Effective Half Life')  
        except:
            pass
        st.write('https://www.studysmarter.co.uk/explanations/physics/medical-physics/effective-half-life/#:~:text=multiple%20choice%20flashcards-,Effective%20half%2Dlife%20is%20the%20time%20it%20takes%20for%20the,will%20be%20discussed%20later%20on.')


    t_effective = (half_life_isotope*half_life_biol)/(half_life_isotope+half_life_biol)
    st.write(f'Effective half life of your compound with {radioisotope1}: {t_effective:.2f} hours')


    with st.expander('Recommendation for Sampling Timepoints'):
        try:    
            st.image(Image.open('./pics/T_Effective_timepoints.PNG'), caption='Image 2: Timepoints')  
        except:
            pass
        st.write('https://link.springer.com/article/10.1007/s11307-023-01868-9')
        try:
            st.image(Image.open('./pics/T_Effective_integration.PNG'), caption='Image 3: Integration')  
        except:
            pass
        st.write('https://doi.org/10.1007/s00259-022-05727-7')

    results_timepoints = dict()
    results_timepoints['Te/3']      = t_effective/3.
    results_timepoints['2Te/3']     = 2 * t_effective/3.
    results_timepoints['3Te/2']     = 3 * t_effective/2.
    results_timepoints['3Te']       = 3 * t_effective
    results_timepoints['5Te']       = 5 * t_effective

    results_timepoints_df = pd.DataFrame.from_dict(results_timepoints, orient='index', columns=['timepoints [h]'])
    results_timepoints_df['timepoints [min]'] = results_timepoints_df['timepoints [h]']*60.
    results_timepoints_df['timepoints [days]'] = results_timepoints_df['timepoints [h]']/24.
    st.write(results_timepoints_df)







        

    


if __name__ == '__main__':
    app(CDD_TOKEN='None')