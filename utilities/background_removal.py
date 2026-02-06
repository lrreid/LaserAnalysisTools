"""
Functions to remove backgrounds from 1D and 2D datasets
"""

import numpy as np
from LaserAnalysisTools.core.rounding import find_nearest

def remove_spectrometer_background(wavelength_data, intenisty_data, wavelength_min, wavelength_max):
    
    index_min = find_nearest(wavelength_data, wavelength_min)[0]
    index_max = find_nearest(wavelength_data, wavelength_max)[0]
    
    background = np.average(intenisty_data[index_min:index_max] )
    int_bg = intenisty_data - background
    
    return background, np.where(int_bg<0, 0, int_bg)

def remove_image_background(image, bg_len, location):
    y_len, x_len = image.shape

    if bg_len > min(image.shape):
        bg_len = min(image.shape) - 1
    
    if location == "top left":
        bg = np.mean(image[0:bg_len, 0:bg_len] ) 
    elif location == "top right":
        bg = np.mean(image[0:bg_len, x_len-bg_len:x_len] )
    elif location == "bottom left":
        bg = np.mean(image[y_len-bg_len:y_len, 0:bg_len] ) 
    elif location == "bottom right":
        bg = np.mean(image[y_len-bg_len:y_len, x_len-bg_len:x_len] ) 
    else:
        print("Lcoation incorrectly defined. Top left used as default")
        bg = np.mean(image[0:bg_len, 0:bg_len] )

    # print("\n\ncalcululated background: %0.2f\n\n" % np.round(bg,2) )
    
    image_bg             = np.round(image - bg,0)
    image_bg_neg_removed =  np.where(image_bg<0, 0, image_bg)
    
    return image_bg_neg_removed