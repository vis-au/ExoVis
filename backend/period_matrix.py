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

"""
def getAlfabetToNumericConverter(saxString, numericSaxConversionArray):
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
"""
def divide_array_in_chunks(list_, elements_in_chunks):
        # looping till length l
        for i in range(0, len(list_), elements_in_chunks):
            yield list_[i:i + elements_in_chunks]

ids_of_exoplanets_to_download= [2,3,4,5,6,7]
time_flux_tuple_arr =[]

def get_lightcurve_data():
    try: 
        time_flux_tuple_arr = np.load('./time_flux_tuple_arr.npy')
    except: 
        print("Lightcurve data is not downloaded")
        print("Downloading Lightcurve data")
        #Download Lighcurve data
        for i in ids_of_exoplanets_to_download:
            kepler_name = "kepler-"+str(i)
            lc = lk.search_lightcurvefile(kepler_name, quarter=0).download().PDCSAP_FLUX.remove_nans()
            time = lc.time    
            fluxes = lc.flux
            norm_fluxes = stats.zscore(fluxes)
            time_flux_tuple = (time, norm_fluxes)
            time_flux_tuple_arr.append(time_flux_tuple)    
            print("finished downloading " + kepler_name)
        np.save('./time_flux_tuple_arr', time_flux_tuple_arr)
    return time_flux_tuple_arr




#transform durations from exoplanet archieve from hours to days
actual_duration_arr=[3.88216/24, 2.36386/24, 3.98235/24 , 4.56904/24 ,3.60111/24, 5.16165/24, 3.19843/24 ] ##kepler-2,3,4,5,6,7,8 https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=cumulative

def get_ground_truth_values(time_flux_tuple_values):
    time_flux_tuple_arr = time_flux_tuple_values
    counter=0
    ground_truth_arr=[]
    #get ground truth values for all lc's with BLS 
    for time_flux_tuple in time_flux_tuple_arr:
        counter +=1
        time = time_flux_tuple[0]
        norm_fluxes = time_flux_tuple[1]
        ## PERIODOGRAM FOR NUMERIC SAX REPRESENTATION
        BLS = BoxLeastSquares(time, norm_fluxes)
        periodogram = BLS.autopower(actual_duration_arr[counter])
        #Find period with highest power in periodogram
        best_period = np.argmax(periodogram.power)  
        period = periodogram.period[best_period]  
        ground_truth_arr.append(period)
    return ground_truth_arr