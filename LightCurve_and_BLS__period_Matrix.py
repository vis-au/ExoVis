#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: sebastianloeschcke
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import altair as alt
import pandas as pd
from scipy import signal
from sklearn.cluster import KMeans
import math


# Astropy, Mast, LightKurve libraries
from astropy.io import fits
from astropy.table import Table 
import matplotlib.pyplot as plt
import lightkurve as lk
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


def getBreakPointsArray(dat_Paa_array, number_of_breakpoints):
    #2d array simulting 1d - Kmeans takes 2d array as input
    x1 = np.ones(len(dat_Paa_array))
    x2 = np.asarray(dat_Paa_array)
    # create new plot and data
    X = np.array(list(zip(x1, x2))).reshape(len(x1), 2)
    # KMeans algorithm 
    K = number_of_breakpoints -1 # one less than alphabet size
    kmeans_model = KMeans(n_clusters=K).fit(X)
    centers = np.array(kmeans_model.cluster_centers_)
    centers_y_values = np.sort(centers[:,1])
    # insert inf and -inf as in saxpy breakpoint arrays
    a = np.array([-np.inf])
    b = np.array([np.inf])
    centers_y_values= np.concatenate((a,centers_y_values, b))
    return centers_y_values


def getNumericSaxArray(centroids = []):
    numeric_SAX_array = []
    for j in range(1,len(centroids)-2):
        tempValue = (centroids[j]+ centroids[j+1])/2
        numeric_SAX_array.append(tempValue)
    firstItem = numeric_SAX_array[0] - (centroids[2]-centroids[1])
    lastItem = numeric_SAX_array[len(numeric_SAX_array)-1] + (centroids[len(centroids)-2]-centroids[len(centroids)-3])
    numeric_SAX_array.insert(0,firstItem)
    numeric_SAX_array.append(lastItem)
    return numeric_SAX_array 
     
    
## fucntion for chooseing alfabet converter
def getAlfabetToNumericConverter(saxString):
    if saxString == "a" :
        return  numericSaxConversionArray[0] 
    if saxString == "b":
        return  numericSaxConversionArray[1] 
    if saxString == "c":
        return numericSaxConversionArray[2] 
    if saxString == "d":
        return  numericSaxConversionArray[3] 
    if saxString == "e":
        return numericSaxConversionArray[4] 
    if saxString == "f":
        return numericSaxConversionArray[5]   
    if saxString == "g":
        return numericSaxConversionArray[6] 
    if saxString == "h":
        return numericSaxConversionArray[7] 
    if saxString == "i":
        return numericSaxConversionArray[8] 
    if saxString == "j":
        return numericSaxConversionArray[9] 
    if saxString == "k":
        return numericSaxConversionArray[10] 
    if saxString == "l":
        return numericSaxConversionArray[11] 
    if saxString == "m":
        return numericSaxConversionArray[12] 
    if saxString == "n":
        return numericSaxConversionArray[13] 
    if saxString == "o":
        return numericSaxConversionArray[14] 
    if saxString == "p":
        return numericSaxConversionArray[15] 
    if saxString == "q":
        return numericSaxConversionArray[16]
    if saxString == "r":
        return numericSaxConversionArray[17] 
    if saxString == "s":
        return numericSaxConversionArray[18] 
    if saxString == "t":
        return numericSaxConversionArray[19] 
    if saxString == "u":
        return numericSaxConversionArray[20] 
    if saxString == "v":
        return numericSaxConversionArray[21] 

def divide_array_in_chunks(list_, elements_in_chunks): 
        # looping till length l 
        for i in range(0, len(list_), elements_in_chunks):  
            yield list_[i:i + elements_in_chunks] 



#array with periods for each parameter combination for PAA/SAX 
period_array = []

import lightkurve as lk
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

for alphabet_size in range(3, 20):
    for paa_division_integer in range(1, 15):
        
        #PAA transformation - 
        paa_points = int(dat_size/paa_division_integer)
        
        ## PAA transformation of data
        PAA_array = paa(dat_znorm, paa_points)
        PAA_array = np.asarray(PAA_array)
        PAA_array = np.float32(PAA_array)
        breakPointsArray = getBreakPointsArray(PAA_array, alphabet_size)
        sax_output = ts_to_string(PAA_array, breakPointsArray)

        ## Convert to numeric representation 
        numericSaxConversionArray = getNumericSaxArray(breakPointsArray)
        numeric_SAX_flux = []
        for i in range(len(sax_output)):
            letter_represented_as_int = getAlfabetToNumericConverter(sax_output[i])
            numeric_SAX_flux.append(letter_represented_as_int)
        numeric_SAX_flux= np.asarray(numeric_SAX_flux)
        numeric_SAX_flux = np.float32(numeric_SAX_flux)

        numeric_SAX_time = time
        # Repeat each element in array x times, where x is the number of PAA points
        repeated_x_array= np.repeat(numeric_SAX_time,paa_points)
        # How many elements each list should have 
        n = int(len(repeated_x_array)/paa_points)
        final_x_array=[]
        lists = list(divide_array_in_chunks(repeated_x_array, n)) 
        #take mean of all chunks
        for l in lists:
            final_x_array.append(np.mean(l))                                 
        numeric_SAX_time= final_x_array

     
        if(alphabet_size==11 and paa_division_integer== 11):
            print("halfway done with computation!")
       
        ## PERIODOGRAM FOR NUMERIC SAX REPRESENTATION
        BLS = BoxLeastSquares(numeric_SAX_time, numeric_SAX_flux)
        periodogram = BLS.autopower(actual_duration_days)
        #Find period with highest power in periodogram
        best_period = np.argmax(periodogram.power)  
        period = periodogram.period[best_period]  
        period_array.append(period)
        

### MATRIX VISUALIZATION
## x and y equal to dimensions of the PAA and Alphabets parameters that have been tested
x, y = np.meshgrid(range(1, 15),range(3, 20))

# Subtract by actual period and round down each value in array 
periods = np.asarray(period_array) - actual_period
periods = np.around(periods, decimals = 6)

# Convert this grid to columnar data expected by Altair
source = pd.DataFrame({'PAA_division_Int': x.ravel(),
                     'SAX_Alfabet_size': y.ravel(),
                     'z': periods})
colorr = alt.Chart(source).mark_rect().encode(
    x='PAA_division_Int:O',
    y='SAX_Alfabet_size:O',
    color=alt.Color('difference:Q',scale=alt.Scale(domain=[0.00112, 1],type='log', base=10))
).transform_calculate(
    difference ="abs(datum.z )"
).properties(
    width=1000,
    height=300
)
text = alt.Chart(source).mark_text().encode(
    x='PAA_division_Int:O',
    y='SAX_Alfabet_size:O',
    text='difference:Q'
).transform_calculate(
    difference ="datum.z"
)
colorr + text
