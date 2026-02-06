"""
Functions for creating histograms.

TO DO:
    - Add option for no weights. 

"""

# import numpy as np
from numpy import histogram, histogram2d
from LaserAnalysisTools.core import binning


def get_histogram_1d(data, bin_res, data_round, weights, equal_bins = False):
    """
    Compute the histogram of a dataset.

    Parameters
    ----------
    data : array_like
        Input data
    bin_res : float
        Width of single histogram bin
    data_round : float
        Rounding applied to max and min data value to form pleasing values for the largest and smallest bin value
    weights : array_like
        An array of weights for the data.
    equal_bins : bool, optional
        If true then Bin_Min = -1*Bin_max, by default False

    Returns
    -------
    hist_1d: array_like
        The values of the histogram
    Bin_Min: float
        Lower bin range
    Bin_Max: float
        Upper bin range
    """

    if equal_bins == False:
        Bin_Min, Bin_Max, Bins = binning.get_bins(data, bin_res, data_round)
    elif equal_bins == True: # makes Bin_Min = -Bin_Max
         Bin_Min, Bin_Max, Bins = binning.get_bins_absolute(data, bin_res, data_round)

    hist_1d = histogram(data, bins=Bins, weights=weights)

    return hist_1d, Bin_Min, Bin_Max, Bins


def get_histogram_2d(data_x, bin_res_x, data_round_x, data_y, bin_res_y, data_round_y, weights, equal_bins_x = False, equal_bins_y = False):
    """
    Compute the bi-dimensional histogram of two data samples.

    Parameters
    ----------
    data_x : array_like
        An array containing the x coordinates of the points to be histogrammed.
    bin_res_x : float
        Width of single histogram bin in the x coordinate
    data_round_x : float
        Rounding applied to max and min x-axis data value to form pleasing values for the largest and smallest bin value
    data_y : array_like
        An array containing the y coordinates of the points to be histogrammed.
    bin_res_y : float
        Width of single histogram bin in the y coordinate
    data_round_y : float
        Rounding applied to max and min y-axis data value to form pleasing values for the largest and smallest bin value
    weights : array_like
        An array of weights for the data.
    equal_bins_x : bool, optional
        If true then Bin_Min = -1*Bin_max, by default False, by default False
    equal_bins_y : bool, optional
        If true then Bin_Min = -1*Bin_max, by default False, by default False

    Returns
    -------
    hist_2d: ndarray
        The bi-dimensional histogram of samples x and y. Values in x are histogrammed along the first dimension and values in y are histogrammed along the second dimension.
    Bin_Min_x: float
        Lower bin range in x coordinate
    Bin_Max_x: float
        Upper bin range in x coordinate
    Bins_x: array_like
        Array of bin values for x coordinate used to make the histogram
    Bin_Min_y: float
        Lower bin range in y coordinate
    Bin_Max_y: float
        Upper bin range in y coordinate
    Bins_y: array_like
        Array of bin values for y coordinate used to make the histogram
    """

    if equal_bins_x == False:
        Bin_Min_x, Bin_Max_x, Bins_x = binning.get_bins(data_x, bin_res_x, data_round_x)
    elif equal_bins_x == True: # makes Bin_Min = -Bin_Max
         Bin_Min_x, Bin_Max_x, Bins_x = binning.get_bins_absolute(data_x, bin_res_x, data_round_x)

    if equal_bins_y == False:
        Bin_Min_y, Bin_Max_y, Bins_y = binning.get_bins(data_y, bin_res_y, data_round_y)
    elif equal_bins_y == True: # makes Bin_Min = -Bin_Max
         Bin_Min_y, Bin_Max_y, Bins_y = binning.get_bins_absolute(data_y, bin_res_y, data_round_y)

    hist_2d, _, _ = histogram2d(data_x, data_y, bins=(Bins_x, Bins_y), weights=weights )

    return hist_2d, Bin_Min_x, Bin_Max_x, Bins_x, Bin_Min_y, Bin_Max_y, Bins_y



