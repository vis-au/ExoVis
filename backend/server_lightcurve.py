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

from pyDOE import *
from scipy.stats.distributions import norm

MIN_PAA = 1
MAX_PAA = 14
MIN_SAX = 3
MAX_SAX = 14

MAX_ITERATIONS = 6

#Contains tuples with parameter combination of selected cells
selected_cells = []

def process_cell(matrix, cell, progress_indicator):
    alphabet_size = cell[0]
    paa_division_integer = cell[1]

    progression = 0 if progress_indicator == -1 else progress_indicator

    # Download or use downloaded lightcurve files
    time_flux_tuple_arr = pm.get_lightcurve_data()

    # get ground truth values for all lc's with autocorrelation
    ground_truth_arr= pm.get_ground_truth_values(time_flux_tuple_arr)

    # transform durations from exoplanet archieve from hours to days
    actual_duration_arr=[3.88216/24, 2.36386/24, 3.98235/24 , 4.56904/24 ,3.60111/24, 5.16165/24, 3.19843/24 ] ##kepler-2,3,4,5,6,7,8 https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=cumulative

    #mean array of periods
    mean_period_arr = []

    # Calculate matrix values for all lighcurves
    # get flux, time, duration and ground thruth for i'th tuple
    ground_truth_period = ground_truth_arr[progression]
    actual_duration = actual_duration_arr[progression]

    time_flux_tuple = time_flux_tuple_arr[progression]
    time = time_flux_tuple[0]
    norm_fluxes = time_flux_tuple[1]

    dat_size = norm_fluxes.size

    # Find Period for eaxh parameter combination alphabets_size/PAA_division_interger of SAX

    # PAA transformation procedure
    # Determine number of PAA points from the datasize devided by the paa_division_integer(number of points per segment)
    paa_points = int(dat_size/paa_division_integer)

    # PAA transformation of data
    PAA_array = paa(norm_fluxes, paa_points)
    PAA_array = np.asarray(PAA_array, dtype=np.float32)

    # SAX conversion
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

    # take mean of all chunks
    for l in lists:
        final_x_array.append(np.mean(l))
    numeric_SAX_time= final_x_array

    # BoxLeastSquares applied to numeric SAX representation
    BLS = BoxLeastSquares(numeric_SAX_time, numeric_SAX_flux)
    periodogram = BLS.autopower(actual_duration)

    # Find period with highest power in periodogram
    best_period = np.argmax(periodogram.power)
    period = periodogram.period[best_period]

    # Add error in percentage between best peiord and ground truth to array with periods
    ground_truth_error = (abs(period - ground_truth_period) / ground_truth_period)*100

    # Update matrix
    if progression == 0:
        matrix[alphabet_size - MIN_SAX][paa_division_integer - MIN_PAA] = ground_truth_error
    else:
        #Update mean of particualr parameter combination
        current_value = matrix[alphabet_size - MIN_SAX][paa_division_integer - MIN_PAA]
        matrix[alphabet_size - MIN_SAX][paa_division_integer - MIN_PAA] = (current_value * progression + ground_truth_error) / (progression+1)


def get_cell_progression(cell, progress_matrix):
    global MIN_PAA, MIN_SAX

    omega = cell[1] - MIN_PAA
    alpha = cell[0] - MIN_SAX

    return int(progress_matrix[alpha][omega])


def get_next_orderly_cell(matrix, progress_matrix, cell_order, index_cell_order):
    global MAX_ITERATIONS

    has_found_valid_cell = False
    candidate_cell = -1

    while not has_found_valid_cell:
        index_cell_order = index_cell_order % len(cell_order)
        candidate_cell = cell_order[index_cell_order]
        index_cell_order += 1
        has_found_valid_cell = get_cell_progression(candidate_cell, progress_matrix) < MAX_ITERATIONS

    return candidate_cell, index_cell_order


def get_next_selected_cell(matrix, progress_matrix, index_selected_cells):
    global MAX_ITERATIONS

    has_found_valid_cell = False
    tried = 0
    candidate_cell = -1

    while (not has_found_valid_cell) and tried < len(selected_cells):
        index_selected_cells = index_selected_cells % len(selected_cells)
        candidate_cell = selected_cells[index_selected_cells]
        index_selected_cells += 1
        has_found_valid_cell = get_cell_progression(candidate_cell, progress_matrix) < MAX_ITERATIONS
        tried += 1

    if has_found_valid_cell:
        return candidate_cell, index_selected_cells
    else:
        return -1, index_selected_cells


def process_next_cell(matrix, progress_matrix, cell_order, index_cell_order, index_selected_cells):
    next_cell = -1

    if len(selected_cells) > 0:
        next_cell, index_selected_cells = get_next_selected_cell(matrix, progress_matrix, index_selected_cells)

    if next_cell == -1:
        next_cell, index_cell_order = get_next_orderly_cell(matrix, progress_matrix, cell_order, index_cell_order)

    progression = get_cell_progression(next_cell, progress_matrix)
    omega = next_cell[1] - MIN_PAA
    alpha = next_cell[0] - MIN_SAX

    process_cell(matrix, next_cell, progression)
    progress_matrix[alpha][omega] = 1 if progression == -1 else progression + 1

    return index_cell_order, index_selected_cells


def get_cell_order():
    global MIN_PAA, MAX_PAA, MIN_SAX, MAX_SAX

    cells = []

    #Create list with parameter combinations alphabets_size/PAA_division_interger for SAX
    for alphabet_size in range(MIN_SAX, MAX_SAX + 1):
        for paa_division_integer in range(MIN_PAA, MAX_PAA + 1):
            cells += [(alphabet_size, paa_division_integer)]

    no_of_samples = (MAX_PAA - MIN_PAA + 1) * (MAX_SAX - MIN_SAX + 1)

    # Do Latin-Hypercube sampling
    lhs_array = lhs(1, samples=no_of_samples, criterion='center')

    #Sort Latin-Hypercube array
    lhs_array_sorted = np.argsort(lhs_array.flatten())
    lhs_final = []

    #Add values of the cells array to lhs_final by using idexes from sorted Latin-Hypercube array
    for index in range(no_of_samples) :
        lhs_final.append(cells[lhs_array_sorted[index]])

    return lhs_final

def update_frontend(matrix, progression_matrix):
    global MAX_ITERATIONS

    eel.sleep(0.01)
    non_negative = progression_matrix.copy()
    non_negative[non_negative == -1] = 0
    column_progression = np.mean(non_negative, axis=0) / MAX_ITERATIONS
    row_progression = np.mean(non_negative, axis=1) / MAX_ITERATIONS
    loadData(matrix.flatten(), column_progression, row_progression, non_negative)


def build_matrix():
    global MIN_PAA, MAX_PAA, MIN_SAX, MAX_SAX, MAX_ITERATIONS, selected_cells

    selected_cells = []
    processed_cell_counter = 0
    cells_to_process = (MAX_PAA - MIN_PAA + 1) * (MAX_SAX - MIN_SAX + 1)

    matrix = np.array([[-1.0 for j in range(MAX_PAA - MIN_PAA + 1)] for _ in range(MAX_SAX - MIN_SAX + 1)])
    progress_matrix = matrix.copy()

    index_selected_cells = 0
    index_cell_order = 0

    while processed_cell_counter < cells_to_process * MAX_ITERATIONS:
        cell_order = get_cell_order()
        index_cell_order = 0

        for _ in range(len(cell_order)):
            index_cell_order, index_selected_cells = process_next_cell(matrix, progress_matrix, cell_order, index_cell_order, index_selected_cells)
            processed_cell_counter += 1

            update_frontend(matrix, progress_matrix)


def run_exo_vis():
    print("Building period matrix ...")
    build_matrix()
    print("Done!")


def loadData(data_arr, avg_progression_columns, avg_progression_rows, progression_matrix):
    global MAX_ITERATIONS

    data_dic = {}
    data_dic["matrix"] = data_arr.tolist()
    data_dic["progression_column"] = avg_progression_columns.tolist()
    data_dic["progression_row"] = avg_progression_rows.tolist()
    data_dic["progression_matrix"] = (progression_matrix / MAX_ITERATIONS).tolist()

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

