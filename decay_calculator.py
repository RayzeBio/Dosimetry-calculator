import numpy as np
import streamlit as st
import matplotlib.pyplot as plt
import pandas as pd

from PIL import Image


# see also https://www.nuclear-power.com/nuclear-engineering/radiation-protection/units-of-radioactivity/becquerel-unit-of-radioactivity/
def calc_enddose(thalf, initialdose, endtime):
    rateconstant = -np.log(1/2)/thalf
    enddose = np.exp(np.log(initialdose)-rateconstant*endtime)
    return enddose

def get_survivors(N_initial, q):
    r = np.random.random(N_initial)
    survive = np.sum(r<q)       # No. of atoms survived
    return survive

def initialize_nucleide(thalf,stable=False):
    lambda_nucleide = np.log(2)/thalf
    if stable:
        p = 0
        q = 1
    else:
        if lambda_nucleide > 1:         # this means that the thalf in this time unit is too small, so probability of decay = 1. To get an accurate estimate, calculate in smaller unit (preferably seconds, but computationally demanding)
            p = 1
            q = 1-p
        else:
            p = lambda_nucleide         # Decay probability, dt =1 t_half_unit
            q = 1 - p                   # Survival probability
    return [0,list(),lambda_nucleide,p,q,0]

# @st.experimental_memo
def get_decay_chain_Ac225(t_half_unit,N,T_end,my_bar):
# def get_decay_chain_Ac225(t_half_unit,N,T_end):
    thalf_actinium225,thalf_fr221,thalf_at217,thalf_bi213,thalf_po213,thalf_pb209 = get_halflives(t_half_unit=t_half_unit,printresults = False) # in minutes

    N_pop1,population_1,lam,p,q,pop1_survive = initialize_nucleide(thalf=thalf_actinium225)

    N_pop2,population_2,lam_pop2,p_pop2,q_pop2,pop2_survive = initialize_nucleide(thalf=thalf_fr221)
    N_pop3,population_3,lam_pop3,p_pop3,q_pop3,pop3_survive = initialize_nucleide(thalf=thalf_at217)
    N_pop4,population_4,lam_pop4,p_pop4,q_pop4,pop4_survive = initialize_nucleide(thalf=thalf_bi213)
    N_pop5,population_5,lam_pop5,p_pop5,q_pop5,pop5_survive = initialize_nucleide(thalf=thalf_po213)
    N_pop6,population_6,lam_pop6,p_pop6,q_pop6,pop6_survive = initialize_nucleide(thalf=thalf_pb209)

    N_pop7,population_7,lam_pop7,p_pop7,q_pop7,pop7_survive = initialize_nucleide(thalf=1.0,stable = True) #Bi209 is stable

    for t in range(T_end):
        survive         = get_survivors(N,q)
        N_pop1_decayed  = N-survive
        population_1.append(survive)

        N_pop2          += N_pop1_decayed
        N_pop2_survive  = get_survivors(N_pop2, q_pop2)
        N_pop2_decayed  = N_pop2 - N_pop2_survive
        N_pop2          -= N_pop2_decayed
        population_2.append(N_pop2_survive)

        N_pop3          += N_pop2_decayed
        N_pop3_survive  = get_survivors(N_pop3, q_pop3)
        N_pop3_decayed  = N_pop3 - N_pop3_survive
        N_pop3          -= N_pop3_decayed
        population_3.append(N_pop3_survive)

        N_pop4          += N_pop3_decayed
        N_pop4_survive  = get_survivors(N_pop4, q_pop4)
        N_pop4_decayed  = N_pop4 - N_pop4_survive
        N_pop4          -= N_pop4_decayed
        population_4.append(N_pop4_survive)

        N_pop5          += N_pop4_decayed
        N_pop5_survive  = get_survivors(N_pop5, q_pop5)
        N_pop5_decayed  = N_pop5 - N_pop5_survive
        N_pop5          -= N_pop5_decayed
        population_5.append(N_pop5_survive)

        N_pop6          += N_pop5_decayed
        N_pop6_survive  = get_survivors(N_pop6, q_pop6)
        N_pop6_decayed  = N_pop6 - N_pop6_survive
        N_pop6          -= N_pop6_decayed
        population_6.append(N_pop6_survive)

        N_pop7          += N_pop6_decayed
        N_pop7_survive  = get_survivors(N_pop7, q_pop7)
        N_pop7_decayed  = N_pop7 - N_pop7_survive
        population_7.append(N_pop7_survive)

        N = survive

        my_bar.progress(t/T_end)

    results = {
        f'time': range(T_end),
        f'Ac-225': population_1,
        f'Fr-221': population_2,
        f'At-217': population_3,
        f'Bi-213': population_4,
        f'Po-213': population_5,
        f'Pb-213': population_6,
        f'Bi-209': population_7
    }

    df_results = pd.DataFrame(results)

    return df_results

#final product: 209Bi is stable
def get_halflives(t_half_unit, printresults=False):    
    thalf_actinium225 = 9.92*24*60
    thalf_fr221 = 4.9 
    thalf_at217 = 32.3e-3/60
    thalf_bi213 = 45.6
    thalf_po213 = 4.2e-6/60
    thalf_pb209 = 3.3 *60

    if t_half_unit == 'days':
        thalf_actinium225 /= 24*60
        thalf_fr221 /= 24*60
        thalf_at217 /= 24*60
        thalf_bi213 /= 24*60
        thalf_po213 /= 24*60
        thalf_pb209 /= 24*60
    elif t_half_unit == 'hours':
        thalf_actinium225 /= 60
        thalf_fr221 /= 60
        thalf_at217 /= 60
        thalf_bi213 /= 60
        thalf_po213 /= 60
        thalf_pb209 /= 60
    elif t_half_unit == 'minutes':
        pass

    if printresults:        
        print(f'half life Act 225: {thalf_actinium225} {t_half_unit}')
        print(f'half life Fr 221: {thalf_fr221} {t_half_unit}')
        print(f'half life At 217: {thalf_at217} {t_half_unit}')
        print(f'half life Bi 213: {thalf_bi213} {t_half_unit}')
        print(f'half life Po 213: {thalf_po213} {t_half_unit}')
        print(f'half life Pb 209: {thalf_pb209} {t_half_unit}')
    return [thalf_actinium225,thalf_fr221,thalf_at217,thalf_bi213,thalf_po213,thalf_pb209]

def plot_Ac225_decay(decay_chain_actinium_plot,t_half_unit, N):
    percentages_decay_products = decay_chain_actinium_plot.values[-1]
    percentages_decay_products = 100*percentages_decay_products/N
    percentages_decay_products = percentages_decay_products.tolist()[1:]      


    labels = [
        '{:44}'.format('Ac-225, thalf = 9.92 days'),
        '{:47}'.format('Fr-221, thalf = 4.9 min'),
        '{:45}'.format('At-217, thalf = 32.3 ms'),
        '{:45}'.format('Bi-213, thalf = 45.6 min'),
        '{:47}'.format('Po-213, thalf = 4.2 us'),
        '{:47}'.format('Pb-213, thalf = 3.3 h'),
        '{:52}'.format('Bi-209, stable')
        ]

    for nucl,perc_nucl in enumerate(percentages_decay_products):
        labels[nucl] = f'{labels[nucl]}, {perc_nucl:>5.2f}%'    

    timelist = np.asarray(decay_chain_actinium_plot.iloc[:,0])
    population_1 = np.asarray(decay_chain_actinium_plot.iloc[:,1])
    population_2 = np.asarray(decay_chain_actinium_plot.iloc[:,2])
    population_3 = np.asarray(decay_chain_actinium_plot.iloc[:,3])
    population_4 = np.asarray(decay_chain_actinium_plot.iloc[:,4])
    population_5 = np.asarray(decay_chain_actinium_plot.iloc[:,5])
    population_6 = np.asarray(decay_chain_actinium_plot.iloc[:,6])
    population_7 = np.asarray(decay_chain_actinium_plot.iloc[:,7])

    timelistdata = np.asarray(timelist)


    colors = plt.cm.BuPu(np.linspace(0.3, 0.8, len(labels)))
    width = 0.8

    fig, ax = plt.subplots()

    timelist = range(len(timelistdata))


    ax.bar(timelist, population_1, width, label=labels[0], color = colors[0])
    ax.bar(timelist, population_2, width, bottom = population_1, label=labels[1], color = colors[1])

    bottom3 = np.array(population_1) + np.array(population_2)
    ax.bar(timelist, population_3, width, bottom = bottom3, label=labels[2], color = colors[2])

    bottom4 = bottom3 + np.array(population_3)
    ax.bar(timelist, population_4, width, bottom = bottom4, label=labels[3], color = colors[3])
    bottom5 = bottom4 + np.array(population_4)
    ax.bar(timelist, population_5, width, bottom = bottom5, label=labels[4], color = colors[4])
    bottom6 = bottom5 + np.array(population_5)
    ax.bar(timelist, population_6, width, bottom = bottom6, label=labels[5], color = colors[5])
    bottom7 = bottom6 + np.array(population_6)
    ax.bar(timelist, population_7, width, bottom = bottom7, label=labels[6], color = colors[6])


    ax.set_ylabel('Number of atoms')
    ax.set_title('Decay chain')
    ax.set_xticks(np.arange(len(timelistdata)))
    if t_half_unit == 'minutes':
        timelistdata = ['{:.2f}'.format(x) for x in (np.asarray(timelistdata)/60/24)]
        t_half_unit = 'days'
    ax.set_xticklabels(timelistdata)
    ax.set_xlabel(f'time [{t_half_unit}]')
    plt.legend(bbox_to_anchor =(1, 1.25))

    # plt.show()
    st.pyplot(fig)

def app(CDD_TOKEN = 'None'):
    st.header("Ac-225 decay calculation")

    with st.expander("Sources"):
        st.write("based on: https://digitalshowcase.lynchburg.edu/jaupli-b/vol3/iss1/2/")
        st.write("see also online calculator: http://www.radprocalculator.com/Decay.aspx or http://www.radprocalculator.com/TimedDecay.aspx")

    with st.expander("Show Ac-225 decay scheme"):
        pic_Ac225decay = Image.open('./pics/1_Ac225-decay.png')
        st.image(pic_Ac225decay, caption='Image 1: Ac-225 with decay half-lives from Radchenko et al. (2015), DOI: 10.1016/j.chroma.2014.12.045')


    #all half lives in minutes !

    thalf_actinium225 = 9.92*24*60
    #10 days from http://nucleardata.nuclear.lu.se/toi/nuclide.asp?iZA=890225&sortG=D&sortA=E
    #9.92 days from Radchenko et al. 2015

    thalf_fr221 = 4.9 
    #4.9min from Radchenko et al. 2015

    thalf_at217 = 32.3e-3/60
    #32.3ms from Radchenko et al. 2015

    thalf_bi213 = 45.6
    #45.6min from Radchenko et al. 2015

    thalf_po213 = 4.2e-6/60
    #4.2 us from Radchenko et al. 2015

    thalf_pb209 = 3.3 *60
    #3.3h from Radchenko et al. 2015

    #final product: 209Bi is stable

    with st.expander("Endpoint dose calculation"):
        col1,col2 = st.columns(2)
        with col1:
            initialdose = st.text_input("Enter initial Ac-225 concentration", 10)
            initialdose_unit = st.text_input("Enter unit of initial Ac-225 concentration", 'mCi')
            submitted = st.button("Calculate")
        with col2:
            endtime = float(st.text_input("Enter time endpoint", 10))
            time_unit = st.radio(label= "Unit of time for endpoint calculation", options=['min','h','d'])

        if submitted:
            if time_unit == 'h':
                endtimecalc = endtime* 60
            elif time_unit == 'd':
                endtimecalc = endtime* 60*24
            else:
                endtimecalc = endtime
            enddose = calc_enddose(thalf=float(thalf_actinium225), initialdose=float(initialdose), endtime=float(endtimecalc))

            # st.write(f't(1/2) for Actinium-225: {thalf_actinium225:.2f} min')
            thalf_actinium225_d = thalf_actinium225 / 60 / 24
            st.write(f't(1/2) for Actinium-225: {thalf_actinium225_d:.2f} d')
            st.write(f'Your final dose after {endtime} {time_unit}: **{enddose:.2f} {initialdose_unit}**')

    if True:
        col1,col2 = st.columns(2)
        with col1:
            N = int(st.text_input("Enter initial number of atoms Ac-225", 1000))
            #future work: implement conversion from Curie to number of atoms
            # N = int(st.text_input("Enter initial activity of Ac-225 [mCi]", 10))
            
            submitted_Ac_decay = st.button("Calculate Ac decay scheme")
        with col2:
            T_end = int(st.text_input("Enter time endpoint", 4))
            t_end_unit = st.radio(label= "Unit of time for endpoint calculation", options=['days','hours','minutes','seconds'])

        # if True:  
        if submitted_Ac_decay:  
            my_bar = st.progress(0)  
            if t_end_unit == 'minutes':
                t_half_unit ='minutes'
            elif t_end_unit == 'hours':
                T_end *= 60
                t_half_unit ='minutes'
            elif t_end_unit == 'days':
                T_end *= 60*24
                t_half_unit ='minutes'
            elif t_end_unit == 'seconds':
                t_half_unit ='seconds'     

            
            # decay_chain_actinium =  get_decay_chain_Ac225(t_half_unit,N,T_end)
            decay_chain_actinium =  get_decay_chain_Ac225(t_half_unit,N,T_end,my_bar)

            # bin_size = int(st.number_input('select bin size', 5, 50, 5, 5))
            bin_size = 5
            bins_to_select = int(len(decay_chain_actinium)/bin_size)
            decay_chain_actinium_plot = decay_chain_actinium[::bins_to_select]
            last_row = decay_chain_actinium.tail(1)
            # try:
            #     decay_chain_actinium_plot = decay_chain_actinium_plot.append(last_row, ignore_index=True)
            # except:
            decay_chain_actinium_plot = pd.concat([decay_chain_actinium_plot,last_row])


            plot_Ac225_decay(decay_chain_actinium_plot,t_half_unit, N)

            # future work: implement average over 1000 runs




        

    


if __name__ == '__main__':
    app()