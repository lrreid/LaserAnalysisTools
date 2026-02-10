"""
Read and analyse text files containing laser pointing data from a far field measurement. 
Files are assumed to be a .txt file with a matrix of pointing values in x, y and total.
"""

import numpy as np
from os.path import splitext
from os import listdir
import matplotlib.pyplot as plt

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

def do_pointing_analysis(data_x, data_y, centre_data: bool = True, rounding = 5):
    # This function needs testing. 

    mean_x, rms_x, centered_data_x = analyse_pointing_data_1D(pointing_data = data_x, centre_data = centre_data)
    mean_y, rms_y, centered_data_y = analyse_pointing_data_1D(pointing_data = data_y, centre_data = centre_data)

    rms_total = np.sqrt( (rms_x**2) + (rms_y**2) )
    plt_lims  = plt_limits_absolute(np.concatenate((centered_data_x, centered_data_y)), rounding = rounding)

    return [mean_x, rms_x, centered_data_x], [mean_y, rms_y, centered_data_y], rms_total, plt_lims

def plot_pointing(txt_dir: str, magnification: float = 1, centre_data: bool = True, rounding = 5, Loop_files: bool = False, file_no=0, Save_Plots: bool = False, DPI = 150, fsize = 12):

    if Save_Plots is True:
        from matplotlib import use
        use("Agg")                  # "Agg" makes it possible to save plots without a display attached. Useful for analysis on remote computing cluster or saving plots in a loop without them opening.

    text_files = [splitext(f)[0] for f in listdir(txt_dir) if f.endswith('.txt')]        # List all files is folder that end in .txt

    if Loop_files is True:
        Q = len(text_files)
    else:
        Q = len([0])

    for i in range(Q):

        if Loop_files is False:
            i = file_no 

        skiprows_value      = find_skiprows("%s\\%s.txt" % (txt_dir, text_files[i]))
        point_x, point_y, _ = np.loadtxt((open("%s\\%s.txt" % (txt_dir, text_files[i]),'rt').readlines()), usecols=(1, 2, 3), skiprows=skiprows_value+1, unpack=True)
        
        data_x, data_y, rms_total, plt_lims = do_pointing_analysis(point_x/magnification, point_y/magnification, centre_data = centre_data, rounding = rounding)

        # Add plot
        fig, ax = plt.subplots(figsize=[6, 6])
        ax.set_aspect('equal')
        plt.plot(data_x[2], data_y[2], '.',color='dodgerblue', label="Pointing data")
        plt.axis([plt_lims[0], plt_lims[1], plt_lims[0], plt_lims[1]])         # set axes [xmin, xmax, ymin, ymax]
        ax.set_title('%s\nRMS stability: %0.2f $\\mu$rad' % (text_files[i], np.round(rms_total, 2)), fontsize=fsize)
        ax.set_xlabel('x ($\\mu$rad)', fontsize=fsize)
        ax.set_ylabel('y ($\\mu$rad)', fontsize=fsize)
        plt.grid(True)
        plt.tight_layout()

        if Save_Plots == True:
            plt.savefig(("%s\\%s_plot.png" % (txt_dir, text_files[i]) ),dpi=DPI,format="png")
            plt.close()
        else:
            plt.show()

    return None

