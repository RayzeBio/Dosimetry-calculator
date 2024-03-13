import streamlit as st
import pandas as pd
import math
import io
import time




# from https://github.com/sorgerlab/gr50_tools

def compute_gr(data):
    """Compute Growth Response value for an entire dataset.

    The input dataframe must contain at least the following numeric fields:

    * cell_count: Number of cells detected per sample.
    * cell_count__time0: Number of cells in the treatment_duration=0 control for each sample.
    * cell_count__ctrl: Number of cells in the no-perturbation control.

    The input must not already contain a column named 'GRvalue'.

    A new dataframe will be returned with the GR values stored in a new
    'GRvalue' column.

    Parameters
    ----------
    data : pandas.DataFrame
        Input data on which to compute the metrics.

    Returns
    -------
    pandas.DataFrame
        Shallow copy of input data with a 'GRvalue' column appended.

    Example
    -------

    """
    if 'GRvalue' in data:
        raise ValueError("Data already contains a 'GRvalue' column; aborting")
    result = data.copy(deep=False)
    result['GRvalue'] = data.apply(compute_gr_single, axis=1)
    return result

def compute_gr_single(record):
    """Compute Growth Response value for a single sample.

    The input is a namedtuple or pandas Series with at least the following
    numeric fields:

    * cell_count: Number of cells detected in this sample.
    * cell_count__time0: Number of cells in the treatment_duration=0 control for this sample.
    * cell_count__ctrl: Number of cells in the no-perturbation control.

    Parameters
    ----------
    record : Union[namedtuple, pandas.Series]
        Input data on which to compute the GR value.

    Returns
    -------
    float
        The computed GR value.

    Example
    -------
    >>> from collections import namedtuple
    >>> Record = namedtuple('Record',
    ...              ['cell_count', 'cell_count__ctrl', 'cell_count__time0'])
    >>> rec = Record(cell_count=1710, cell_count__ctrl=1766.0,
    ...              cell_count__time0=492.8)
    >>> print compute_gr_single(rec)
    0.965305500206
    """
    cc_t0 = float(record.cell_count__time0)
    log2nn = _normalize_log2(float(record.cell_count), cc_t0)
    log2nn_ctrl = _normalize_log2(float(record.cell_count__ctrl), cc_t0)
    gr = 2 ** (log2nn / log2nn_ctrl) - 1
    return gr
    
def _normalize_log2(n, n_0_0):
    normalized = max(n / n_0_0, 1e-6) # avoiding negative and null values
    return math.log(normalized, 2)

def convert_timelist_to_GRdataframe(cell_line,agent,concentration,cell_count,timepoint,cell_count__time0,cell_count__ctrl):
    datalist = list()
    if not len(cell_count) == len(cell_count__ctrl):
        print('Error in input for list conversion')
    else:
        for i,t in enumerate(timepoint):
            sublist = [str(cell_line)]
            sublist.append(str(agent))
            sublist.append(str(concentration))
            sublist.append(str(cell_count[i]))
            sublist.append(str(t))
            sublist.append(str(cell_count__time0))
            sublist.append(str(cell_count__ctrl[i]))

            datalist.append(sublist)
        GRdataframe = pd.DataFrame(datalist,
        columns=['cell_line','agent', 'concentration', 'cell_count','timepoint','cell_count__time0','cell_count__ctrl']
        )
        return GRdataframe
  
def df_to_excel(df,filename):
    output = io.BytesIO()
    with pd.ExcelWriter(output,mode ='w',engine='xlsxwriter') as writer:
        workbook  = writer.book
        df.to_excel(writer)
        writer.save()

        st.download_button(
            label='download results',
            data= output,
            file_name=filename
            )



def app(CDD_TOKEN='None'):
    st.header("GR calculation")
    st.write("Hafner, M., Niepel, M. Chung, M. and Sorger, P.K., Metrics of drug sensitivity and resistance based on growth rate inhibition correct for the confounding effects of variable division rates, (2016) Nature Methods, doi:10.1038/nmeth.3853")

    columns_input = st.columns(4)
    with columns_input[1]:
        cell_line = st.text_input("Cell line", "T84")
        agent = st.text_input("Compound", '7228-Lu177')
        agent_conc = st.text_input("Concentration of compound [nM] / radioactivity [uCi]", '0.1')
        cell_count_time0 = st.text_input("tumor volume [mm3] before experiment (time 0)", '400')
        timepoints_length = st.slider("How many timepoints?", 1, 20, 5)

    example_times = [2,5,14,27,38]
    example_control = [420,450,500,750,1000]
    example_cells = [410,420,440,580,700]


    input_columns = st.columns(4)
    timepoints = list()
    cell_count_list = list()
    cell_count_control_list = list()
    for i in range(timepoints_length):
        i+=1
        with input_columns[0]:
            st.markdown('***')
            time_name = st.text_input('time [days] %i'%(i), str(example_times[i-1]))
            timepoints.append(time_name)
        with input_columns[1]:
            st.markdown('***')
            cell_count = st.text_input(f'tumor volume [mm3]  {i}', str(example_cells[i-1]))
            cell_count_list.append(cell_count)
        with input_columns[2]:
            st.markdown('***')
            cell_count_control = st.text_input(f'tumor volume [mm3] of control {i}', str(example_control[i-1]))
            cell_count_control_list.append(cell_count_control)          


    if st.button('Calculate GR'):
        gr_dataframe = convert_timelist_to_GRdataframe(cell_line=cell_line, agent=agent, concentration=agent_conc, cell_count=cell_count_list, timepoint=timepoints, cell_count__time0=cell_count_time0,cell_count__ctrl=cell_count_control_list)
        results_table = compute_gr(gr_dataframe)

        st.dataframe(results_table.style.format({'GRvalue':'{:.2f}'}))

        date_save = time.localtime()
        date_save = str(date_save.tm_year)+str(date_save.tm_mon)+str(date_save.tm_mday)
        filenameGR = date_save+'_GR-results-'+agent+'.xlsx'
        df_to_excel(results_table,filenameGR)



if __name__ == '__main__':
    app()