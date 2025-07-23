# calculator
IRDC (Inhouse Rayzebio Dosimetry Calculator)
Written by Charlotte Lorenz, RayzeBio Inc. 

Version 2.0: 2025-07-23
updated:
- Ac225 with all daughters in decay chain included for dose deposition
- female and male phantoms integrated for scaling and dose deposition

Version 0.7: 2023-01-16
updated:
- Scaling methods validated with Beykan et al. publication

Version 0.6: 2023-01-09
updated:
- Data input from file: corrected readout from unassigned column headers
- added rawdata columns in case of decay correction
- added calculation of alpha for biexponential with absorption and elimination

Version 0.5: 2022-12-20
updated:
- removed scaling options 'no scaling' and 'relative mass and time - FDA'
- removed decay correction for alpha scaling

Version 0.4: 2022-11-11
updated:
- colors in CDD input rawdata
- bone marrow calculation
- blood, muscle tissue weights from standard 25g mouse model (not CDD input)
- %ID/g from file input
- from CDD Vault: select tissues new st.session state
- muscle weight changed to 28000g
- dosimetry results download file was changed to be used as CDD Vault upload file (22.11.22)

Version 0.3: 2022-11-04
updated:
- Dose calculation: hTIAC(org) * 60 * 60 to adjust for df in MBq-s
- Decay correction for input possible 
- From CDD Vault: Tissue and total body weight mouse from CDD Vault 

Version 0.2: 2022-10-26
updated:
- other dose limiting organs added
- tumor dose from 310g sphere (Olinda)


Version 0.1: 2022-10-04
presented in dosimetry company talk.
Features:
1) BioD request from CDD Vault
2) Uploaded Conditions/Uploaded tissues
3) BioD data Input
4) Fitting decay
- monoexponential direct method
- monoexponential lin. reg. on log(y)
- biexponential
- biexponential with absorption and elimination
- trapezoidal
- trapezoidal with linear extrapolation
5) Human extrapolation
- mass scaling
- time scaling
- combined time and mass scaling
- no scaling
- alpha scaling
6) Dosimetry
7) Results

General To Dos:
- calculate bone marrow dose based on blood
- write version number on title and store it in excel output
- create CDD Vault protocol for upload of results
- Make request based on RAYZ ID
- Store masses per tissue in excel
- cross validate with publications
- calculate total body dose
- show if ICPMS data are available?
- calculate dose from all scaling methods and give range of results (?)
