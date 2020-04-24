"""
This file is derived from the eel example page. It sets up a basic eel environment that also
works with React and weppacked frontend.
This allows us to use npm modules in the frontend (i.e. React and D3) and at the same time gives
us TypeScript support in the frontend.

src:
https://github.com/samuelhwilliams/Eel/tree/master/examples/07%20-%20CreateReactApp

"""

import os
import platform
import random
import sys
import lightkurve as lk
import eel

import math
import matplotlib.pyplot as plt
import altair as alt
import numpy as np
import pandas as pd
from scipy import signal

# Astropy, Mast, LightKurve libraries
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
from scipy import stats
from astroquery.mast import Mast
from astroquery.mast import Observations
from astropy import units as u
from astropy.timeseries import BoxLeastSquares

#for sax conversion
from saxpy.znorm import znorm
from saxpy.paa import paa
from saxpy.sax import ts_to_string
from saxpy.alphabet import cuts_for_asize
from saxpy.sax import sax_via_window

from astropy.timeseries import TimeSeries
from time import process_time
import time

import period_matrix as pm

#Lating Hypercube imports
from pyDOE import *
from scipy.stats.distributions import norm

#Contains tuples with parameter combination of selected cells
selected_cells = []

def run_exo_vis():
    print("Building period matrix ...")

    global selected_cells
    #Download or use downloaded lightcurve files
    time_flux_tuple_arr = pm.get_lightcurve_data()
    #get ground truth values for all lc's with autocorrelation
    ground_truth_arr= pm.get_ground_truth_values(time_flux_tuple_arr)
    #transform durations from exoplanet archieve from hours to days
    actual_duration_arr=[3.88216/24, 2.36386/24, 3.98235/24 , 4.56904/24 ,3.60111/24, 5.16165/24, 3.19843/24 ] ##kepler-2,3,4,5,6,7,8 https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=cumulative
    #mean array of periods
    mean_period_arr = []
    #Initialize 12*14 matrix to store data being send to frontend
    matrix = np.array([[-1.0 for j in range(14)] for _ in range(12)])
    #Matrix to show progression of each cell in percentage
    progression_matrix = matrix.copy()

    #Calculate matrix values for all lighcurves
    for i in range(len (time_flux_tuple_arr)):
        # get flux, time, duration and ground thruth for i'th tuple
        ground_truth_period = ground_truth_arr[i]
        actual_duration = actual_duration_arr[i]
        time_flux_tuple = time_flux_tuple_arr[i]
        time = time_flux_tuple[0]
        norm_fluxes = time_flux_tuple[1]
        dat_size = norm_fluxes.size
        #Variable to keep count of number of cells finished in matrix
        cell_counter = 0
        cells = []

        #Create list with parameter combinations alphabets_size/PAA_division_interger for SAX
        for alphabet_size in range(3, 15):
            for paa_division_integer in range(1, 15):
                cells += [(alphabet_size, paa_division_integer)]

        # Do Latin-Hypercube sampling
        lhs_array = lhs(1, samples=168, criterion='center')
        #Sort Latin-Hypercube array
        lhs_array_sorted = np.argsort(lhs_array.flatten())
        lhs_final = []
        #Add values of the cells array to lhs_final by using idexes from sorted Latin-Hypercube array
        for index in range(168) :
            lhs_final.append(cells[lhs_array_sorted[index]])
        cells = lhs_final

        c = []
        #Find cells selected and put in front of array
        for selected_cell in selected_cells:
            #Remove and append to front of cells array
            for cell in cells:
                if cell[0] == selected_cell[0] and cell[1] == selected_cell[1]:
                    c += [cell]

        for selected_cell in c:
            cells.remove(selected_cell)
            cells[:0] = [selected_cell]

        #Find Period for eaxh parameter combination alphabets_size/PAA_division_interger of SAX
        for cell in range(len(cells)):
            alphabet_size = cells[cell][0]
            paa_division_integer = cells[cell][1]

            #PAA transformation procedure
            #Determine number of PAA points from the datasize devided by the paa_division_integer(number of points per segment)
            paa_points = int(dat_size/paa_division_integer)
            ## PAA transformation of data
            PAA_array = paa(norm_fluxes, paa_points)
            PAA_array = np.asarray(PAA_array, dtype=np.float32)

            #SAX conversion
            # Get breakpoints to convert segments into SAX string
            breakPointsArray = pm.getBreakPointsArray(PAA_array, alphabet_size)
            sax_output = ts_to_string(PAA_array, breakPointsArray)

            # Convert to numeric SAX representation
            numericSaxConversionArray = pm.getNumericSaxArray(breakPointsArray)
            numeric_SAX_flux = []
            for symbol_index in range(len(sax_output)):
                letter_represented_as_int = pm.getAlfabetToNumericConverter(sax_output[symbol_index], numericSaxConversionArray)
                numeric_SAX_flux.append(letter_represented_as_int)
            numeric_SAX_flux= np.asarray(numeric_SAX_flux,dtype=np.float32)
            numeric_SAX_time = time
            # Repeat each element in array x times, where x is the number of PAA points
            repeated_x_array= np.repeat(numeric_SAX_time,paa_points)
            # How many elements each list should have
            n = int(len(repeated_x_array)/paa_points)
            final_x_array=[]
            lists = list(pm.divide_array_in_chunks(repeated_x_array, n))
            #take mean of all chunks
            for l in lists:
                final_x_array.append(np.mean(l))
            numeric_SAX_time= final_x_array

            #BoxLeastSquares applied to numeric SAX representation
            BLS = BoxLeastSquares(numeric_SAX_time, numeric_SAX_flux)
            periodogram = BLS.autopower(actual_duration)
            #Find period with highest power in periodogram
            best_period = np.argmax(periodogram.power)
            period = periodogram.period[best_period]

            #Add error in percentage between best peiord and ground truth to array with periods
            ground_truth_error = (abs(period - ground_truth_period) / ground_truth_period)*100
            #Update mean periods array
            if i == 0:
                matrix[alphabet_size - 3][paa_division_integer - 1] = ground_truth_error
            else:
                #Update mean of particualr parameter combination
                current_value = matrix[alphabet_size - 3][paa_division_integer - 1]
                matrix[alphabet_size - 3][paa_division_integer - 1] = (current_value * i + ground_truth_error) / (i+1)
            #Update progression of cell in percentage
            progression_matrix[alphabet_size - 3][paa_division_integer - 1] += (1/len(ground_truth_arr))
            #Send mean periods array data to server every 2th iteration with loadData()
            if(cell_counter % 2 ==0):
                eel.sleep(0.01)
                loadData(matrix.flatten(), np.mean(progression_matrix, axis=0) , np.mean(progression_matrix, axis=1),progression_matrix ) #load mean period array if full
            cell_counter +=1
        loadData(matrix.flatten(), np.mean(progression_matrix, axis=0) , np.mean(progression_matrix, axis=1) ,progression_matrix) #load mean period array if full
        print("matrix ", i, " finished")


def loadData(data_arr, avg_progression_columns, avg_progression_rows,progression_matrix):
#def loadData(data_arr):
    data_dic = {}
    data_dic["matrix" ]= data_arr.tolist()
    data_dic["progression_column" ]= avg_progression_columns.tolist()
    data_dic["progression_row" ]= avg_progression_rows.tolist()
    data_dic["progression_matrix" ]= progression_matrix.tolist()
    #eel.send_data(data.tolist())
    eel.send_data(data_dic)

@eel.expose
def send_selected_cells(new_selected_cells):
    global selected_cells
    print("received new cells")
    selected_cells = new_selected_cells

@eel.expose
def register_client(message):
    print(message)

    starttime = time.time()
    interval = 0.25

    eel.spawn(run_exo_vis)

def start_eel(develop):
    """Start Eel with either production or development configuration."""
    print("Launching backend ...")

    if develop:
        directory = '../frontend/'
        app = None
        page = {'port': 3000}
    else:
        directory = 'build'
        app = 'chrome-app'
        page = 'index_lightcurve.html'

    eel.init(directory, ['.tsx', '.ts', '.jsx', '.js', '.html'])

    print('Backend launched successfully. Waiting for requests ...')

    # These will be queued until the first connection is made, but won't be repeated on a page reload
    # eel.say_hello_js('Python World!')   # Call a JavaScript function (must be after `eel.init()`)

    eel_kwargs = dict(
        host='localhost',
        port=8080,
        size=(1280, 800),
    )
    try:
        eel.start(page, mode=None, **eel_kwargs)
    except EnvironmentError:
        # If Chrome isn't found, fallback to Microsoft Edge on Win10 or greater
        if sys.platform in ['win32', 'win64'] and int(platform.release()) >= 10:
            eel.start(page, mode='edge', **eel_kwargs)
        else:
            raise

if __name__ == '__main__':
    import sys

    # Uses the production version in the "build" directory if passed a second argument
    start_eel(develop=len(sys.argv) == 1)

