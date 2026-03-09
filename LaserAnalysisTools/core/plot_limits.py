"""
Define functions for calculating plot axes limits
"""

import numpy as np
from LaserAnalysisTools.core.rounding import roundup, rounddn

def get_data_order_of_magnitude(value):
    """
    Get the order of magnitude of a number

    Parameters
    ----------
    value : float
        Number that you want to get the order of magnitude

    Returns
    -------
    order: int
        Order of magnitude of value
    """
    order = np.floor(np.log10(value)) # order of magnitude of max

    return order

def plt_limits(array, rounding):
    """
    Calculate upper and lower plot limit axes for a given dataset.

    Parameters
    ----------
    array : array_like
        List/numpy array containing data.
    rounding : float
        Number that you want to round the data to.

    Returns
    -------
    plt_min : float
        Value to be used as lower plot limit.
    plt_max : float
        Value to be used as upper plot limit.

    """
    
    plt_min = rounddn(np.min(array), rounding)
    plt_max = roundup(np.max(array), rounding)
    
    return plt_min, plt_max

def plt_limits_absolute(array, rounding):
    """
    Calculate upper and lower plot limit axes for a given dataset where the min = -1*max.

    Parameters
    ----------
    array : array_like
        List/numpy array containing data.
    rounding : float
        Number that you want to round the data to.

    Returns
    -------
    plt_min : float
        Value to be used as lower plot limit.
    plt_max : float
        Value to be used as upper plot limit.
    """
    
    plt_max = roundup( np.max(abs(array)) , rounding)

    return -plt_max, plt_max


def plt_limits_log(array, min_offset: int = 0, max_offset: int = 0):
    """
    Calculate upper and lower plot limit axes for a given dataset rounded to the nearest order of magnitude.
    Name is a bit confusing. Called "Log" as it was initially written for finding good looking plot limits on log axis plots.

    Parameters
    ----------
    array : array_like
        List/numpy array containing data.
    min_offset : int, optional
        Offset on order of magnitude rounding for lower plot limit, by default 0
    max_offset : int, optional
        Offset on order of magnitude rounding for upper plot limit, by default 0

    Returns
    -------
    plt_min : float
        Value to be used as lower plot limit.
    plt_max : float
        Value to be used as upper plot limit.
    """
    
    if np.min(array) == 0:
        plt_min = 1                                 # order of magnitude of min
    elif np.min(array) < 0:
        data_inv_min = np.max(-1*array)
        order_min = get_data_order_of_magnitude(data_inv_min)
        plt_min   = -1*roundup(data_inv_min, 10**(order_min+min_offset))
    else:
        data_min  = np.min(array)
        order_min = get_data_order_of_magnitude(data_min) # order of magnitude of max
        plt_min   = rounddn(data_min, 10**(order_min+min_offset))
        
    data_max  = np.max(array)
    order_max = get_data_order_of_magnitude(data_max) # order of magnitude of max
    plt_max   = roundup(data_max, 10**(order_max+max_offset))

    return plt_min, plt_max


def plt_limits_log_absolute(array, offset: int = 0):
    """
    Calculate upper and lower plot limit axes for absolute value of given dataset rounded to the nearest order of magnitude.
    Name is a bit confusing. Called "Log" as it was initially written for finding good looking plot limits on log axis plots.

    Parameters
    ----------
    array : _type_
        List/numpy array containing data.
    min_offset : int, optional
        Offset on order of magnitude rounding for plot limit, by default 0

    Returns
    -------
    plt_max : float
        Value to be used as upper plot limit.
    """

    if np.max(array) == 0:
        plt_max = 1 
    else:
        data_max = np.max(abs(array))
        order    = get_data_order_of_magnitude(data_max)
        plt_max  = roundup(data_max, 10**(order+offset))

    return plt_max