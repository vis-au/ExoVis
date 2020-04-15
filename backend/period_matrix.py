#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: sebastianloeschcke
"""

import numpy as np
from sklearn.cluster import KMeans

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

def divide_array_in_chunks(list_, elements_in_chunks):
        # looping till length l
        for i in range(0, len(list_), elements_in_chunks):
            yield list_[i:i + elements_in_chunks]