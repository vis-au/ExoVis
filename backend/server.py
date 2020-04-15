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


def run_exo_vis():
    print("building period matrix ...")

    #array with periods for each parameter combination for PAA/SAX
    period_array = np.ones(238, dtype=int)
    period_array = [element * 5 for element in period_array]

    ## Search for lightcurve file of Exoplanet with LightKurve library - choose the corrected PDCSAP_FLUX and remove NaNs

    lc = lk.search_lightcurvefile('kepler-8', quarter=0).download().PDCSAP_FLUX.remove_nans()
    fluxes = lc.flux
    time = lc.time

    actual_period = 3.5224 ## found in NASA Database: https://exoplanetarchive.ipac.caltech.edu/cgi-bin/DisplayOverview/nph-DisplayOverview?objname=K00010.01&type=KEPLER_CANDIDATE
    actual_duration_hours =3.1984  ## found in NASA Database: https://exoplanetarchive.ipac.caltech.edu/cgi-bin/DisplayOverview/nph-DisplayOverview?objname=K00010.01&type=KEPLER_CANDIDATE
    actual_duration_days=  0.13327 #3.19843/24

    dat_size= fluxes.size
    #Z-Normalize data
    dat_znorm = stats.zscore(fluxes)
    cell_counter = 0

    for alphabet_size in range(3, 20):
        for segment_size in range(1, 15):

            ###PAA transformation procedure
            #Determine number of PAA points from the datasize devided by the segment_size(number of points per segment)
            paa_points = int(dat_size/segment_size)

            ## PAA transformation of data
            PAA_array = paa(dat_znorm, paa_points)
            PAA_array = np.asarray(PAA_array)
            PAA_array = np.float32(PAA_array)
            # Get breakpoints to convert segments into SAX string
            breakPointsArray = pm.getBreakPointsArray(PAA_array, alphabet_size)
            sax_output = ts_to_string(PAA_array, breakPointsArray)

            ## Convert to numeric representation
            numericSaxConversionArray = pm.getNumericSaxArray(breakPointsArray)
            numeric_SAX_flux = []

            for i in range(len(sax_output)):
                letter_represented_as_int = pm.getAlfabetToNumericConverter(sax_output[i], numericSaxConversionArray)
                numeric_SAX_flux.append(letter_represented_as_int)

            numeric_SAX_flux= np.asarray(numeric_SAX_flux)
            numeric_SAX_flux = np.float32(numeric_SAX_flux)

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


            if(alphabet_size==11 and segment_size== 11):
                print("halfway done with computation!")

            ## PERIODOGRAM FOR NUMERIC SAX REPRESENTATION
            BLS = BoxLeastSquares(numeric_SAX_time, numeric_SAX_flux)
            periodogram = BLS.autopower(actual_duration_days)
            #Find period with highest power in periodogram
            best_period = np.argmax(periodogram.power)
            period = periodogram.period[best_period]
            #period_array.append(period)
            period_array[cell_counter]=period
            if(cell_counter % 5 ==0):
               loadData(period_array)
            cell_counter +=1
    loadData(period_array)
    print("done building period matrix")

    """
     x, y = np.meshgrid(range(1, 15),range(3, 20))
    # Subtract by actual period and round down each value in array
    periods = np.asarray(period_array) - actual_period
    periods = np.around(periods, decimals = 6)

    # Convert this grid to columnar data expected by Altair
    source = pd.DataFrame({'PAA_division_Int': x.ravel(),
                        'SAX_Alfabet_size': y.ravel(),
                        'z': periods})

    #eel.send_data_to_frontend(source)
    """

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
    data = np.asarray(data_arr) #np.load("./data/mean_period_arr_100%.npy")
    #print(data)
    #print(data.tolist())
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
        page = 'index.html'

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
