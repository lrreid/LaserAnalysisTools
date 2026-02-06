"""
Read and analyse text files containing laser pointing data from a far field measurement. 
Files are assumed to be a .txt file with a matrix of pointing values in x, y and total.
"""

import numpy as np
from os.path import splitext
from os import listdir
from LaserAnalysisTools.core.statistics import w_std
from LaserAnalysisTools.core.plot_limits import plt_limits_absolute
from LaserAnalysisTools.utils.find_skiprows import find_skiprows

def analyse_pointing_data_1D(pointing_data, centre_data: bool = True):
    mean, rms = w_std( pointing_data, None )

    if centre_data == True:
        # Centre data to mean position.
        centered_data = pointing_data - mean
        mean, rms     = w_std( centered_data, None ) # update mean position of pointing after centering, should be almost exactly zero.
    else:
        centered_data = pointing_data
    
    return mean, rms, centered_data

def analyse_pointing_data(data_x, data_y, centre_data: bool = True, rounding = 5):

    mean_x, rms_x, centered_data_x = analyse_pointing_data_1D(pointing_data = data_x, centre_data = centre_data)
    mean_y, rms_y, centered_data_y = analyse_pointing_data_1D(pointing_data = data_y, centre_data = centre_data)

    rms_total = np.sqrt( (rms_x**2) + (rms_y**2) )
    plt_max   = plt_limits_absolute(np.concatenate((centered_data_x, centered_data_y)), rounding = rounding)

    return [mean_x, rms_x, centered_data_x], [mean_y, rms_y, centered_data_y], rms_total, plt_max

def do_pointing_analysis(txt_dir, file_number, magnification = 1, centre_data: bool = True, rounding = 5):

    text_files = [splitext(f)[0] for f in listdir(txt_dir) if f.endswith('.txt')]        # List all files is folder that end in .txt
    print("Reading file: %s\n" % (text_files[file_number]))

    skiprows_value      = find_skiprows("%s\\%s.txt" % (txt_dir, text_files[file_number]))
    point_x, point_y, _ = np.loadtxt((open("%s\\%s.txt" % (txt_dir, text_files[file_number]),'rt').readlines()), usecols=(1, 2, 3), skiprows=skiprows_value+1, unpack=True)
    
    data_x, data_y, rms_total, plt_max = analyse_pointing_data(point_x/magnification, point_y/magnification, centre_data = centre_data, rounding = rounding)

    return data_x, data_y, rms_total, plt_max

