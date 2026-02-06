
"""
Functions to define bins for creating histograms of data. 
"""

import numpy as np
from LaserAnalysisTools.core.rounding import roundup, rounddn

def get_bins(array, bin_res, bin_round):
    """
    Get numpy array of containing bins to be used for creating histograms of data.

    Parameters
    ----------
    array : array_like
        List/numpy array containing data.
    bin_res : float
        Resolution/bin width of data for histogram
    bin_round : float
        Value to round up and down the max and min values of the array to find max and min bins.

    Returns
    -------
    bins_min: float
        Minimum bin value
    bins_max: float
        Maximum bin value
    hist_bins: numpy array
        array of histogram bins
    """

    bins_max  = roundup(np.max(array), bin_round)
    bins_min  = rounddn(np.min(array), bin_round)
    hist_bins = np.arange(bins_min, bins_max+bin_res, bin_res)

    return bins_min, bins_max, hist_bins

def get_bins_absolute(array, bin_res, bin_round):
    """
    Get numpy array of containing bins to be used for creating histograms of data where the min = -1*max.

    Parameters
    ----------
    array : array_like
        List/numpy array containing data.
    bin_res : float
        Resolution/bin width of data for histogram
    bin_round : float
        Value to round up and down the max and min values of the array to find max and min bins.

    Returns
    -------
    bins_min: float
        Minimum bin value
    bins_max: float
        Maximum bin value
    hist_bins: numpy array
        array of histogram bins
    """

    bins_max  = roundup(np.max(abs(array)), bin_round)
    hist_bins = np.arange(-bins_max, bins_max+bin_res, bin_res)

    return -bins_max, bins_max, hist_bins