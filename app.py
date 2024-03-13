import streamlit as st
# import gr_calculator
import decay_calculator
import dosimetry_idg
import dosimetry_hTIACg
import dosimetry_hTIACorg
import dosimetry_CDD
import dosimetry_CDD_imaging
import scaling_mTIAC

# dictionary of all available apps that can be accessed on the main streamlit page
PAGES = {
    "dosimetry: from %ID/g": dosimetry_idg,
    "dosimetry: from hTIAC(g)": dosimetry_hTIACg,
    "dosimetry: from hTIAC(org)": dosimetry_hTIACorg,
    "dosimetry: from CDD BioD": dosimetry_CDD,
    "dosimetry: from CDD Imaging": dosimetry_CDD_imaging,
    "scaling: mTIAC -> hTIAC": scaling_mTIAC,
    # "GR":   gr_calculator,
    "Radioisotope decay": decay_calculator,    
}

st.sidebar.title('Navigation')

# Check if URL has parameters (selection/CDD_TOKEN): e.g. URL/?CDD_TOKEN=None&selection=GR
params = st.query_params
if "selection" in params:
    try:    
        index_selection = list(PAGES.keys()).index(params["selection"])
    except:
        index_selection=0
else:
    index_selection = 0

if "CDD_TOKEN" in params:
    try:    
        CDD_TOKEN = params["CDD_TOKEN"]
    except:
        CDD_TOKEN='None'
else:
    CDD_TOKEN = 'None'

# radio menu to select from all available PAGES
selection = st.sidebar.radio("Apps", list(PAGES.keys()),index=index_selection)
page = PAGES[selection]

# run app() function of selected app
page.app(CDD_TOKEN)