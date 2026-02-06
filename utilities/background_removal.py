"""
Functions to remove backgrounds from 1D and 2D datasets
"""

import numpy as np


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