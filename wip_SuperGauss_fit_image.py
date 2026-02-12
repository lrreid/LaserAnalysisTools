

def SuperGauss_fit_image(image, smooth_data: bool = False, smooth_sigma: int = 5, Lineout_Bandwidth: bool = False, Lineout_depth: int = 2, init_vals_x = None, init_vals_y = None, calibration: float = 1):
    """
    I want to update how the lineout bandwidths work so you define the total width not either side of the centroid.

    I shouldn't do everything for x and y separetly. This should be moved to a function and called twice.
    """
    
    if smooth_data is True:
        image = gaussian_filter( image, sigma = smooth_sigma) # will this cause problems? Do I want to use .copy() to copy image but not reference it

    index_max = np.where(image == np.max(image))
    centroid  = D4S_centroid(image)


    if Lineout_Bandwidth is True:
        x_lineout = np.sum(image[int(round(centroid[1],0)) - Lineout_depth : int(round(centroid[1],0)) + Lineout_depth+1 ,:], axis=0) / ((2*Lineout_depth)+1) # Divide normalises data back to single pixel intensity values
        y_lineout = np.sum(image[:, int(round(centroid[0],0)) - Lineout_depth : int(round(centroid[0],0)) + Lineout_depth+1], axis=1) / ((2*Lineout_depth)+1)
    else:    
        x_lineout = image[int(round(centroid[1],0)),:]
        y_lineout = image[:,int(round(centroid[0],0))]

    x_pixels = np.arange(0, len(x_lineout), 1)
    y_pixels = np.arange(0, len(y_lineout), 1)

    if init_vals_x is None:                                                                     # If guess values for fit not provided, use automatic guesses.
        init_vals_x = [0.75*max(x_lineout), int(round(centroid[0],0)), len(x_lineout)//2, 4]    # initial (guess) values for [amplitude, centroid, width, order]

    if init_vals_y is None:
        init_vals_y = [0.75*max(y_lineout), int(round(centroid[1],0)), len(y_lineout)//2, 4]    # initial (guess) values for [amplitude, centroid, width, order]

    x_fit, Amp_x, order_x, centroid_x, FWHM_x = SuperGauss_lst_sq(x_data=x_pixels, y_data=x_lineout, guess_vals=init_vals_x)
    y_fit, Amp_y, order_y, centroid_y, FWHM_y = SuperGauss_lst_sq(x_data=y_pixels, y_data=y_lineout, guess_vals=init_vals_y)

    sizes_x = superGauss_beam_conversions_FWHM(order = order_x, FWHM = FWHM_x/calibration)
    sizes_y = superGauss_beam_conversions_FWHM(order = order_y, FWHM = FWHM_y/calibration)

    x_pixels_norm = np.arange(0, len(x_lineout), 1) - centroid_x
    y_pixels_norm = np.arange(0, len(y_lineout), 1) - centroid_y

    x_mm = x_pixels_norm/calibration
    y_mm = y_pixels_norm/calibration

    return [x_fit, Amp_x, order_x, centroid_x, FWHM_x], [y_fit, Amp_y, order_y, centroid_y, FWHM_y], [x_mm, y_mm], sizes_x, sizes_y





