'''
Statistical calculations

TO DO:
    - Create 1D version of D4S centroid calculation.

'''


import numpy as np

def w_std( a, weights ):
    """
    Calculate the weighted standard deviation of a dataset.
    Adapted from w_std from openpmd_viewer.addons.pic.lpa_diagnostics

    Parameters
    ----------
    a : array_like
        Calculate the weighted standard deviation for these a.

    weights : array_like
        An array of weights for the values in a.

     Returns
     -------
     average: float
         Average value of dataset
    std_dev: float
        Weighted standard deviation of dataset
    """

    average = np.average(a, weights=weights)
    std_dev = np.sqrt(np.average((a - average) ** 2, weights=weights))
    return average, std_dev


def D4S_centroid(image, rtn_int: bool = False):
    """
    Find the coordinates of the centroid of a 2D array according to a D4-sigma method

    Parameters
    ----------
    image : ndarray
        image data
    rtn_int : bool, optional
        Round coordinates to closest value and return and integer, by default False

    Returns
    -------
    xc: float
        centroid of data in x coordinate
    yc: float
        centroid of data in y coordinate
    """
    v, h = image.shape

    p = np.sum(image)                       # float avoids integer overflow
    # find the centroid
    hh = np.arange(h, dtype=float)          # float avoids integer overflow
    vv = np.arange(v, dtype=float)          # ditto
    xc = sum(np.dot(image, hh))/p
    yc = sum(np.dot(image.T, vv))/p
    
    if rtn_int is True:
        return int(np.round(xc,0)), int(np.round(yc,0))
    else:
        return xc, yc