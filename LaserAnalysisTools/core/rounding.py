"""
Rounding functions
"""

import numpy as np

def roundup(x,n):
    """
    Round up a number up to the closest value of n

    Parameters
    ----------
    x : float
        Number that you want to appy rounding to
    n : float
        Value that you want to round to

    Returns
    -------
    float
        x rounded up to the closest value of n
    """
    return int(np.ceil(x / n)) * n

def rounddn(x,n):
    """
    Round up a number down to the closest value of n

    Parameters
    ----------
    x : float
        Number that you want to appy rounding to
    n : float
        Value that you want to round to

    Returns
    -------
    float
        x rounded down to the closest value of n
    """
    return int(np.floor(x / n)) * n


def round_nearest(x, n):
    """
    Round a number to the closest value of n

    Parameters
    ----------
    x : float
        Number that you want to appy rounding to
    n : float
        Value that you want to round to

    Returns
    -------
    float
        x rounded to the closest value of n
    """
    return n * np.round(x/n)

def normalise(array):
    """
    Normalise a dataset

    Parameters
    ----------
    array : array_like
        Input dataset

     Returns
    -------
    array_like
        Normalised dataset  
    """
    peak = np.max(abs(array))
    norm = array/peak
    return(norm)


def find_nearest(array, target):
    """
    Find the value in an array closest to a target value

    TO DO:
        What does this function do if it finds two equally good candidates?

    Parameters
    ----------
    array : array_like
        Input dataset_
    target : float
        Target value that you want to find in the dataset

    Returns
    -------
    idx: int
        Index in array corresponding to value closest to the target
    value: float
        Value in array closest to the target.
    """
    array = np.asarray(array)
    idx = (np.abs(array - target)).argmin()
    return idx, array[idx]