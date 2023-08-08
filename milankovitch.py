############################################# IMPORTING NEEDED PACKAGES #############################################
import numpy as np
from scipy.stats import zscore

############################################# LOADING MILANKOVITCH DATA ##############################################

#Data file has columns are time, latitude, insolation, eccentricity, precession, longitude, obliquity
dir_milan = '' #directory where data is stored june 21
milan_data = np.loadtxt(dir_milan) #loading data

yrs = np.sort(milan_data[:,0])/1000 #data converted to units of kyr
insolation = zscore(milan_data[:,2]) #data zero meaned


def milankovitch(t, paillard_truncation = False):
    '''
    Input: time (t) in units of kyrs
    Output: insolation at that time
    This function takes the Berger & Loutre 1991 insolation values (variable name insolation)
    for the last 10 million years and interpolates for any given time (t) in kyrs
    '''

    milankovitch_interpolated = np.interp(t, yrs , insolation) 
    
    if paillard_truncation != False:
    #pallaird truncation of insolation timeseries
        a = 1
        milankovitch_pt = 0.5*(milankovitch_interpolated + np.sqrt(4*a**2 + milankovitch_interpolated**2))
        average = 1.1167871593192462 #1.1097818163097832
        rms = 0.5292873319224887
        Ihat_trunc = (milankovitch_pt - average)/rms
        min_Ihat = 1.45875129556558
        max_Ihat = 5.546414620499402
        I_trunc = (Ihat_trunc + min_Ihat)/(max_Ihat) #nondimensionalizes the truncated time series
        return I_trunc #nondimensionalized and paillard truncation applied for time t
    else:
        max_ins = max(insolation)
        min_ins = min(insolation)
        I = (milankovitch_interpolated - min_ins)/(max_ins - min_ins) #nondimensionalizing time series
        return I #nondimensionalized insolation for time t
