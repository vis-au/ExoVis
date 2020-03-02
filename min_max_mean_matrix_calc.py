#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 09:54:42 2020

@author: sebastianloeschcke
"""

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
     
    
## function for choosing alfabet converter
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

from statsmodels import api as sm
import lightkurve as lk
## Search for lightcurve file of Exoplanet with LightKurve library - choose the corrected PDCSAP_FLUX and remove NaNs

#get 10 lighcurve files from 10 known exoplanets

ids_of_exoplanets_to_download= [2,3,4,5,6,7,8,9,10,11]
time_flux_tuple_arr =[]

#array with periods for each parameter combination for PAA/SAX 
periods = []
#array holding periods for every lc file
periods_array=[]

for i in ids_of_exoplanets_to_download:
    kepler_name = "kepler-"+str(i)
    lc = lk.search_lightcurvefile(kepler_name, quarter=2).download().PDCSAP_FLUX.remove_nans()
    time = lc.time    
    fluxes = lc.flux
    norm_fluxes = stats.zscore(fluxes)
    time_flux_tuple = (time, norm_fluxes)
    time_flux_tuple_arr.append(time_flux_tuple)    

        
#get ground truth values for all lc's with autocorrelation 
ground_truth_arr=[]
for time_flux_tuple in time_flux_tuple_arr:
    
    time = time_flux_tuple[0]
    norm_fluxes = time_flux_tuple[1]
    # get the autocorrelation coefficient
    acf = sm.tsa.acf(norm_fluxes, nlags=len(norm_fluxes))
    
    period = np.where(acf[20:] == np.amax(acf[20:]))
    if(period[0].size > 1):
        #take average of period if several values in period array
        resPeriod= sum(map(float, filter(None, period[0][1:])))/(len(period[0])-1) +20
    else:
        resPeriod = str(period[0] +20)
        resPeriod = int(resPeriod[1:-1])
    avg_distance_in_time_arr = (time[(time.size)-1]-time[0])/(time.size-1) 
    resPeriod= resPeriod*avg_distance_in_time_arr
    actual_period = resPeriod
    ground_truth_arr.append(actual_period)


#calculate matrix values for all lighcurves
for i in range(len (time_flux_tuple_arr)):
    # get flux, time, and ground thruth for i'th tuple
    ground_truth_period = ground_truth_arr[i]
    time_flux_tuple = time_flux_tuple_arr[i]
    time = time_flux_tuple[0]
    norm_fluxes = time_flux_tuple[1]
    dat_size = norm_fluxes.size
    
    for alphabet_size in range(3, 15): 
        for paa_division_integer in range(1, 15):
            
            ###PAA transformation procedure
            #Determine number of PAA points from the datasize devided by the paa_division_integer(number of points per segment) 
            paa_points = int(dat_size/paa_division_integer)
            
            ## PAA transformation of data
            PAA_array = paa(norm_fluxes, paa_points)
            PAA_array = np.asarray(PAA_array)
            PAA_array = np.float32(PAA_array)
            # Get breakpoints to convert segments into SAX string
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
    
          
            # declaring magnitude of repetition 
            K = paa_division_integer
            # using list comprehension 
            # repeat elements K times 
            res =  [ele for ele in numeric_SAX_flux for i in range(K)] 
            # get the autocorrelation coefficient
            tempNorm = res
            acf = sm.tsa.acf(tempNorm, nlags=len(tempNorm))
    
            #lag = np.arange(0,len(tempNorm), 0.5)
            lag = np.arange(len(tempNorm)) 
    
            period = np.where(acf[20:] == np.amax(acf[20:]))
            if(period[0].size > 1):
                resPeriod= sum(map(float, filter(None, period[0][1:])))/(len(period[0])-1) +20
            else:
                #print(" paa:"+str(q)+ " alfa"+str(alfa)+ str(period[0]))
                resPeriod = str(period[0] +20)
                resPeriod = int(resPeriod[1:-1])
            avg_distance_in_time_arr = (time[(time.size)-1]-time[0])/(time.size-1) 
            period= resPeriod*avg_distance_in_time_arr
            periods.append(period)
    periods_array.append(periods)
    periods=[]
        
     


ten_ints_to_compare = []
min_period_arr=[]
max_period_arr=[]
mean_period_arr=[]
for i in range(len(periods_array[0])):
    ten_ints_to_compare = []
    for array in periods_array:
        ten_ints_to_compare.append(array[i])
    min_period_arr.append(min(ten_ints_to_compare))
    max_period_arr.append(max(ten_ints_to_compare))
    mean_period_arr.append(np.mean(ten_ints_to_compare))

   

#generate matrix vizualization
x, y = np.meshgrid(range(1, 15),range(3, 15))

# Subtract by actual period and round down each value in array 
periods = np.asarray(period_array) - actual_period
periods = np.around(periods, decimals = 5)
min_value = min(np.abs(periods))
if min_value == 0:
    min_value +=0.001
max_value = max(np.abs(periods))

# Convert this grid to columnar data expected by Altair
source = pd.DataFrame({'Segment size ω': x.ravel(),
                     'Alphabet size α': y.ravel(),
                     'z': periods})
color = alt.Chart(source).mark_rect().encode(
    x='Segment size ω:O',
    y='Alphabet size α:O',
    color=alt.Color('difference:Q',scale=alt.Scale(domain=[min_value, max_value],type='log', base=10))
).transform_calculate(
    difference ="abs(datum.z )"
).properties(
    width=1000,
    height=300
)
text = alt.Chart(source).mark_text(

    baseline='middle',
    
    ).encode(
        x='Segment size ω:O',
        y='Alphabet size α:O',
        text='difference:Q',
    color=alt.condition(
            alt.datum.z/max_value < 0.10,
            alt.value('black'),
            alt.value('white')
        )
        
    ).transform_calculate(
difference ="datum.z"
   )

color + text