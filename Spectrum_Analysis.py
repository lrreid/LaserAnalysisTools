"""
Functions to analyse laser spectrum data

None of these have been tested.
"""

import numpy as np
from os.path import splitext
from os import listdir
from LaserAnalysisTools.core.rounding import find_nearest
from LaserAnalysisTools.core.statistics import w_std
from LaserAnalysisTools.utils.background_removal import remove_spectrometer_background

def select_wavelength_range(wavelength_data, intensity_data, wavelength_min, wavelength_max):
    # Need to add some error handling. What if selected range is outwith the data in the file?
       
    index_min, lw_min = find_nearest(wavelength_data, wavelength_min)
    index_max, lw_max = find_nearest(wavelength_data, wavelength_max)
    
    if abs(lw_min-wavelength_min) > 5:
        print("Warning! Minimum wavelength >5 nm from user specified value")
        
    if abs(lw_max-wavelength_max) > 5:
        print("Warning! Maximum wavelength >5 nm from user specified value")
    
    wavelength_select = wavelength_data[index_min:index_max]
    intensity_select  = intensity_data[index_min:index_max] 
    
    return wavelength_select, intensity_select, [index_min, index_max]


def central_wavelength_thres(int_threshold, wavelength_data, intenisty_data):
    
    index_peak = np.where(intenisty_data == np.max(intenisty_data))[0][0]
    
    blue_max = intenisty_data[0:index_peak]
    red_max  = intenisty_data[index_peak:-1]
    
    thes_ind_blue = find_nearest(blue_max, int_threshold*np.max(intenisty_data))[0]
    thes_ind_red  = find_nearest(red_max,  int_threshold*np.max(intenisty_data))[0] + index_peak

    central_wavelength = wavelength_data[thes_ind_blue] + (( wavelength_data[thes_ind_red] - wavelength_data[thes_ind_blue] )/2)
    
    return central_wavelength, wavelength_data[thes_ind_blue], wavelength_data[thes_ind_red]


def analyse_spectrum(wavelength_data, intensity_data, int_threshold = 1/np.exp(1)**2):

    mean_wavelength, RMS_Spread = w_std( wavelength_data, intensity_data )
    stat_bandwidth              = 4*RMS_Spread
    # print("\nWeighted average wavelength: %0.2f nm" % (np.round(mean_wavelength,2)))
    # print("Bandwidth: %0.2f nm\n" % (np.round(4*RMS_Spread,2)) )

    central_wavelength_e2, e2_blue, e2_red = central_wavelength_thres(int_threshold=int_threshold, wavelength_data=wavelength_data, intenisty_data=intensity_data)
    theshold_bandwidth                     = e2_red-e2_blue
    # print("\nWavelength at centre of 1/e^2 boundaries: %0.2f nm" % (np.round(central_wavelength_e2,2)))
    # print("Width to threshold intensity: %0.2f nm\n" % (e2_red-e2_blue))

    index_peak      = np.where(intensity_data == np.max(intensity_data))[0][0]
    wavelength_peak = wavelength_data[index_peak]
    # print("\nWavelength at peak intensity: %0.2f nm\n" % (np.round(wavelength_peak,2)))

    return [mean_wavelength, stat_bandwidth], [central_wavelength_e2, theshold_bandwidth, e2_blue, e2_red], wavelength_peak


# This is partially pseudo-code, it needs fixed.
def do_spectrum_analysis(txt_dir, file_number, remove_background: bool = False, bg_limits = [925, 950], select_data: bool = False, select_limits = [650, 950], int_threshold = 1/np.exp(1)**2):
    # I don't think I like this method for looping through lots of data files. This function will re-find a list of all the text files in a folder numerous times. 
    # I would also have to find text_files to see the number of text files in a given folder to see how many times to loop. Think about stacking functions in a different way?

    text_files            = [splitext(f)[0] for f in listdir(txt_dir) if f.endswith('.txt')]        # List all files is folder that end in .txt
    wavelength, intensity = np.loadtxt((open("%s\\%s.txt" % (txt_dir, text_files[file_number]),'rt').readlines()[:-1]), skiprows=0, unpack=True)

    if remove_background is True:
        _, intensity = remove_spectrometer_background(wavelength, intensity, bg_limits[0], bg_limits[1])

    if select_data is True:
        wavelen_select, intensity_select, _ = select_wavelength_range(wavelength, intensity, select_limits[0], select_limits[1])

    [mean_wavelength, stat_bandwidth], [central_wavelength_e2, theshold_bandwidth, e2_blue, e2_red], wavelength_peak = analyse_spectrum(wavelength_data = wavelen_select, intensity_data = intensity_select, int_threshold = int_threshold)

    return [wavelen_select, intensity_select], [mean_wavelength, stat_bandwidth], [central_wavelength_e2, theshold_bandwidth, e2_blue, e2_red], wavelength_peak