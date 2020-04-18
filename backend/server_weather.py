
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

import lin_regression_matrix as pm

#Lating Hypercube imports
from pyDOE import *
from scipy.stats.distributions import norm

def run_exo_vis():
    print("Building period matrix ...")
    
    #Download or use downloaded weather data files
    time_avg_temp_arr = pm.get_weather_data_from_files()
    #get ground truth for the slope computed by the lin. regression model for all weather data files.
    ground_truth_arr= pm.get_ground_truth_values(time_avg_temp_arr)
    #mean array of periods
    mean_period_arr = []
    #Initialize 12*12 matrix to store data being send to frontend
    matrix = np.array([[0.0 for j in range(12)] for _ in range(12)])

    #Calculate matrix values for all lighcurves
    for i in range(len (time_avg_temp_arr)):
        # get flux, time, duration and ground thruth for i'th tuple
        ground_truth_slope = ground_truth_arr[i]
        time_avg_temp_tuple = time_avg_temp_arr[i]
        time = time_avg_temp_tuple[0]
        norm_avg_temp = time_avg_temp_tuple[1]
        dat_size = norm_avg_temp.size
        #Variable to keep count of number of cells finished in matrix
        cell_counter = 0
        cells = []

        #Create list with parameter combinations alphabets_size/PAA_division_interger for SAX
        for alphabet_size in range(3, 15):
            for paa_division_integer in range(1, 13):
                cells += [(alphabet_size, paa_division_integer)]

        #Contains tuples with parameter combination of selected cells
        list_of_selected_cells = []
        
        # Do Latin-Hypercube sampling
        lhs_array = lhs(1, samples=144, criterion='center')
        #Sort Latin-Hypercube array 
        lhs_array_sorted = np.argsort(lhs_array.flatten())
        lhs_final = []
        #Add values of the cells array to lhs_final by using idexes from sorted Latin-Hypercube array
        for index in range(144) :
            lhs_final.append(cells[lhs_array_sorted[index]])
        cells = lhs_final
        
        #Find cells selected and put in front of array
        for selected_cell in list_of_selected_cells:
            #Remove and append to front of cells array
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
            PAA_array = paa(norm_avg_temp, paa_points)
            PAA_array = np.asarray(PAA_array, dtype=np.float32)
            
            #SAX conversion
            # Get breakpoints to convert segments into SAX string
            breakPointsArray = pm.getBreakPointsArray(PAA_array, alphabet_size)
            sax_output = ts_to_string(PAA_array, breakPointsArray)
            
            # Convert to numeric SAX representation
            numericSaxConversionArray = pm.getNumericSaxArray(breakPointsArray)
            numeric_SAX_temp = []
            for symbol_index in range(len(sax_output)):
                letter_represented_as_int = pm.getAlfabetToNumericConverter(sax_output[symbol_index], numericSaxConversionArray)
                numeric_SAX_temp.append(letter_represented_as_int)
            numeric_SAX_temp= np.asarray(numeric_SAX_temp,dtype=np.float32)
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

            #Compute linear regression model
            slope, intercept, r_value, p_value, std_err = stats.linregress(numeric_SAX_time,numeric_SAX_temp)
            #Add error in percentage between best peiord and ground truth to array with periods
            ground_truth_error = (abs(slope - ground_truth_slope) / ground_truth_slope)*100
            #Update mean periods array
            if i == 0:
                matrix[alphabet_size - 3][paa_division_integer - 1] = ground_truth_error
            else:
                #Update mean of particualr parameter combination
                current_value = matrix[alphabet_size - 3][paa_division_integer - 1]
                matrix[alphabet_size - 3][paa_division_integer - 1] = (current_value * i + ground_truth_error) / (i+1)
            #Send mean periods array data to server every 2th iteration with loadData()
            if(cell_counter % 2 ==0):
                loadData(matrix.flatten()) #load mean period array if full
            cell_counter +=1
        loadData(matrix.flatten())
        print("matrix ", i, " finished")
    print(matrix)

def get_cell_order(columns, rows):
    cell_order = []

    while len(cell_order) < (columns * rows):
        next_col = int(random.uniform(0, columns))
        next_row = int(random.uniform(0, rows))
        next_cell = (next_col, next_row)

        if next_cell not in cell_order:
            cell_order += [next_cell]

    return cell_order


def get_updated_value(cell):
    return random.uniform(0, 1)


def send_progressive_updates(cell_order, counter):
    steps = 5
    updated_data = []

    for _ in range(0, steps):
        next_cell = cell_order[counter]
        preliminary_value = get_updated_value(next_cell)

        updated_data += [{
            'cell': next_cell,
            'value': preliminary_value
        }]

        counter += 1
        if counter == len(cell_order):
            counter = 0

    eel.send_data_to_frontend(updated_data)
    return counter


def loadData(data_arr):
    data = np.asarray(data_arr) 
    eel.send_data(data.tolist())

@eel.expose
def register_client(message):
    print(message)

    starttime = time.time()
    interval = 0.25

    run_exo_vis()
    #loadData()

    # cell_order = get_cell_order(10, 10)
    # progressiveness_counter = 0

    # while True:
    #     print("tick")
    #     time.sleep(interval - ((time.time() - starttime) % interval))
    #     progressiveness_counter = send_progressive_updates(cell_order, progressiveness_counter)



def start_eel(develop):
    """Start Eel with either production or development configuration."""
    print("hello there")
    if develop:
        directory = '../frontend/'
        app = None
        page = {'port': 3000}
    else:
        directory = 'build'
        app = 'chrome-app'
        page = 'index_weather.html'

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

