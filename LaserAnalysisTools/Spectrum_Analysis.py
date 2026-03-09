"""
Functions to analyse laser spectrum data

None of these have been tested.
"""

import sys
import numpy as np
from os.path import splitext
from os import listdir
import matplotlib.pyplot as plt

from LaserAnalysisTools.core.rounding import find_nearest
from LaserAnalysisTools.core.statistics import w_std
from LaserAnalysisTools.utils.background_removal import remove_spectrometer_background
from LaserAnalysisTools.core.plot_limits import plt_limits_log, plt_limits

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

def do_spectrum_analysis(wavelength, intensity, remove_background: bool = False, bg_limits = [925, 950], select_data: bool = False, select_limits = [650, 950], int_threshold = 1/np.exp(1)**2):

    if remove_background is True:
        _, intensity = remove_spectrometer_background(wavelength, intensity, bg_limits[0], bg_limits[1])

    if select_data is True:
        wavelen_select, intensity_select, _ = select_wavelength_range(wavelength, intensity, select_limits[0], select_limits[1])

    [mean_wavelength, stat_bandwidth], [central_wavelength_e2, theshold_bandwidth, e2_blue, e2_red], wavelength_peak = analyse_spectrum(wavelength_data = wavelen_select, intensity_data = intensity_select, int_threshold = int_threshold)

    return [wavelen_select, intensity_select], [mean_wavelength, stat_bandwidth], [central_wavelength_e2, theshold_bandwidth, e2_blue, e2_red], wavelength_peak


def plot_spectrum(txt_dir: str, remove_background: bool = False, bg_limits = [925, 950], select_data: bool = False, select_limits = [650, 950], int_threshold = 1/np.exp(1)**2, Loop_files: bool = False, file_no=0, Save_Plots: bool = True, DPI = 150, fsize = 12):

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

        wavelength, intensity = np.loadtxt((open("%s\\%s.txt" % (txt_dir, text_files[i]),'rt').readlines()[:-1]), skiprows=0, unpack=True)

        data_select, stat, bandwidth, wavelength_peak = do_spectrum_analysis(wavelength, intensity, remove_background=remove_background, bg_limits=bg_limits, select_data=select_data, select_limits=select_limits, int_threshold=int_threshold)

        _, Int_max = plt_limits_log(data_select[1])      # Find max intenisty for plot limit
        Int_min = 0
        plt_min, plt_max = plt_limits(data_select[0], 10)

        fig, ax = plt.subplots()
        plt.plot([bandwidth[2], bandwidth[2]], [0, Int_max], '--', linewidth = 2, color="seagreen")
        plt.plot([bandwidth[3], bandwidth[3]], [0, Int_max], '--', linewidth = 2, color="seagreen")
        plt.plot([bandwidth[0], bandwidth[0]],   [0, Int_max], '-', linewidth = 2, color="seagreen", label="centre of 1/e$^{2}$\nboundaries:\n%0.2f nm" % (np.round(bandwidth[0],2)))
        plt.plot([wavelength_peak, wavelength_peak], [0, Int_max], '-', linewidth = 2, color="dodgerblue", label="$\\lambda$ at peak intensity:\n%0.2f nm" % (np.round(wavelength_peak,2)))
        plt.plot([stat[0], stat[0]], [0, Int_max], '-', linewidth = 2, color="magenta", label="weighted average:\n%0.2f nm" % (np.round(stat[0],2)))
        plt.plot(data_select[0], data_select[1], color='firebrick')
        plt.text((0.86*(plt_max-plt_min))+plt_min, (0.62*(Int_max-Int_min))+Int_min, "Width to 1/e$^{2}$\n%0.2f nm" % (bandwidth[1]) ,  color='seagreen', fontweight='normal', fontsize=12, horizontalalignment ="center" )
        ax.set_title('%s' % (text_files[i]), fontsize=fsize)
        
        plt.axis([plt_min, plt_max, Int_min, Int_max])         # set axes [xmin, xmax, ymin, ymax]
        ax.set_xlabel('Wavelength (nm)', fontsize=fsize)
        ax.set_ylabel('Intensity (a.u.)', fontsize=fsize)
        ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0)) # put yticks in scientific notation
        plt.grid(True)
        plt.legend(loc='upper left', prop={'size': 9})
        plt.tight_layout()
        
        if Save_Plots is True:
            plt.savefig(("%s\\%s_Analysed.png" % (txt_dir, text_files[i]) ), dpi=DPI, format="png")
            plt.close()

            # Save text file with summary of analysis
            fp =  open("%s/%s_Analysis_Summary.txt" % (txt_dir, text_files[i]),"w")     # create text file, 'w' for write, 'a' for append.
            fp.write("Script used for analysis:\t%s" % sys.argv[0])
            fp.write("\nFrom statistical analysis of data\n")
            fp.write("Central wavelength:\t%0.2f nm\n" % np.round(stat[0], 2)  )
            fp.write("Bandwidth:\t%0.2f nm\n" % np.round(stat[1], 2)  )
            fp.write("\nFrom intensity threshold analysis method\n")
            fp.write("Intensity threshold at %0.3f of peak intensity\n" % np.round(int_threshold, 3)  )
            fp.write("Central wavelength:\t%0.2f nm\n" % np.round(bandwidth[0], 2)  )
            fp.write("Bandwidth:\t%0.2f nm\n" % np.round(bandwidth[1], 2)  )
            fp.write("\nWavelength at peak intensity:\t%0.2f nm" % np.round(wavelength_peak, 2) )
            fp.close()
        else:
            plt.show()

    return None