
import numpy as np
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter

from LaserAnalysisTools.core.statistics import D4S_centroid
from LaserAnalysisTools.core.distributions import superGauss_beam_conversions_FWHM

def SuperGauss(x_data, Amp, x0, FWHM, order):		# Define function for Guassian distribution
		return Amp * np.exp(-np.log(2)*( ((4*((x_data-x0)**2)) / (FWHM**2))  )**order )
    
def SuperGauss_lst_sq(x_data, y_data, guess_vals):
    Coeff, _ = curve_fit(SuperGauss, x_data, y_data, p0=guess_vals)
    fit_line = SuperGauss(x_data=x_data, Amp=Coeff[0], x0=Coeff[1], FWHM=Coeff[2], order=Coeff[3])
    
    return fit_line, Coeff[0], Coeff[3], Coeff[1], Coeff[2]


def get_size_from_px_intensity(int_data, x_data, target_counts, counts_rounding):
    
    n_points = np.size(int_data)
    
    int_1 = np.where(np.round(int_data[0:n_points//2], counts_rounding) == target_counts)[0][-1]
    int_2 = np.where(np.round(int_data[n_points//2: n_points], counts_rounding) == target_counts)[0][0] + n_points//2
    
    width = x_data[int_2]-x_data[int_1]
    
    return width, x_data[int_1], x_data[int_2]


def SuperGauss_fit_line(lineout, init_vals = None, calibration: float = 1):

    pixels = np.arange(0, len(lineout), 1)

    if init_vals is None:                                                                       # If guess values for fit not provided, use automatic guesses.
        init_vals = [0.75*max(lineout), len(lineout)//2, len(lineout)//3, 4]                    # initial (guess) values for [amplitude, centroid, width, order]

    fit, Amp, order, centroid, FWHM = SuperGauss_lst_sq(x_data=pixels, y_data=lineout, guess_vals=init_vals)

    sizes = superGauss_beam_conversions_FWHM(order = order, FWHM = FWHM/calibration)

    pixels_norm = np.arange(0, len(lineout), 1) - centroid
    line_cal    = pixels_norm/calibration

    return [fit, Amp, order, centroid, FWHM], line_cal, sizes


def SuperGauss_fit_image(image, smooth_data: bool = False, smooth_sigma: int = 5, Lineout_Bandwidth: bool = False, Lineout_depth: int = 2, init_vals_x = None, init_vals_y = None, calibration: float = 1, target_counts: int = 2, counts_rounding: int = 0):
    """
    I want to update how the lineout bandwidths work so you define the total width not either side of the centroid.
    """
    
    if smooth_data is True:
        image = gaussian_filter( image, sigma = smooth_sigma) # will this cause problems? Do I want to use .copy() to copy image for smoothing but not reference it

    centroid = D4S_centroid(image)

    if Lineout_Bandwidth is True:
        x_lineout = np.sum(image[int(round(centroid[1],0)) - Lineout_depth : int(round(centroid[1],0)) + Lineout_depth+1 ,:], axis=0) / ((2*Lineout_depth)+1) # Divide normalises data back to single pixel intensity values
        y_lineout = np.sum(image[:, int(round(centroid[0],0)) - Lineout_depth : int(round(centroid[0],0)) + Lineout_depth+1], axis=1) / ((2*Lineout_depth)+1)
    else:    
        x_lineout = image[int(round(centroid[1],0)),:]
        y_lineout = image[:,int(round(centroid[0],0))]

    fit_x, x_mm, sizes_x = SuperGauss_fit_line(x_lineout, init_vals = init_vals_x, calibration = calibration)
    fit_y, y_mm, sizes_y = SuperGauss_fit_line(y_lineout, init_vals = init_vals_y, calibration = calibration)

    thres_width_x, x_low_bound, x_high_bound = get_size_from_px_intensity(int_data = x_lineout, x_data = x_mm, target_counts = target_counts, counts_rounding = counts_rounding)
    thres_width_y, y_low_bound, y_high_bound = get_size_from_px_intensity(int_data = y_lineout, x_data = y_mm, target_counts = target_counts, counts_rounding = counts_rounding)

    return fit_x, fit_y, [x_mm, y_mm], sizes_x, sizes_y, [thres_width_x, x_low_bound, x_high_bound], [thres_width_y, y_low_bound, y_high_bound]

