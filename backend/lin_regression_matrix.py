#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: sebastianloeschcke
"""

import numpy as np
from sklearn.cluster import KMeans
import string
import lightkurve as lk
from scipy import stats
from astropy.timeseries import BoxLeastSquares


#Returns array with tuples of all loaded weather data files
def get_weather_data_from_files():
    time_flux_tuple_arr = [] 
    time_flux_tuple_arr.append( np.load('./weather_data/dk_year_weather_time_avgTemp.npy'))
    time_flux_tuple_arr.append(np.load('./weather_data/afghan_year_weather_time_avgTemp.npy'))
    time_flux_tuple_arr.append(np.load('./weather_data/swe_year_weather_time_avgTemp.npy')) 
    time_flux_tuple_arr.append(np.load('./weather_data/ned_year_weather_time_avgTemp.npy'))
    time_flux_tuple_arr.append( np.load('./weather_data/ger_year_weather_time_avgTemp.npy')) 
    time_flux_tuple_arr.append(np.load('./weather_data/ch_year_weather_time_avgTemp.npy'))
    time_flux_tuple_arr.append(np.load('./weather_data/fr_year_weather_time_avgTemp.npy'))
    time_flux_tuple_arr.append( np.load('./weather_data/lux_year_weather_time_avgTemp.npy'))
    time_flux_tuple_arr.append(np.load('./weather_data/spain_year_weather_time_avgTemp.npy')) 
    time_flux_tuple_arr.append(np.load('./weather_data/por_year_weather_time_avgTemp.npy')) 
    return time_flux_tuple_arr

#Calculates ground truth value for slope of lin. regression model for each weather data file
def get_ground_truth_values(_time_avg_temp_arr):
    time_avg_temp_arr = _time_avg_temp_arr 
    ground_truth_arr = []
    for weather_tuple in time_avg_temp_arr:
        slope, intercept, r_value, p_value, std_err = stats.linregress(weather_tuple[0],weather_tuple[1])
        ground_truth_arr.append(slope)
    return ground_truth_arr  

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
alphabet_string = string.ascii_lowercase
alphabet_list = list(alphabet_string)

## function for choosing alfabet converter
def getAlfabetToNumericConverter(saxString, numeric_sax_conversion_array):
    numericSaxConversionArray = numeric_sax_conversion_array
    res = [let for let in alphabet_list if let == saxString]
    return numericSaxConversionArray[alphabet_list.index(res[0])]

def divide_array_in_chunks(list_, elements_in_chunks):
        # looping till length l
        for i in range(0, len(list_), elements_in_chunks):
            yield list_[i:i + elements_in_chunks]