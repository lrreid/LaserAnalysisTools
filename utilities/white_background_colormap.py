"""
Functions for creating custom matplotlib colormaps
"""

from matplotlib.colors import ListedColormap
from numpy import asarray, linspace
from matplotlib.pyplot import get_cmap

def alpha_varying_cmap(cmap, thresh: int, thresh_position: str = 'lower'):
    """
    Create a colormap with continuously varying transparency (alpha)

    Parameters
    ----------
    cmap : matplotlib.colors.ListedColormap
        matplotlib colormap. For example: matplotlib.pyplot.cm.inferno
    thresh : int
        Threshold pixel value where colormap transparency begins to change. Should be set between 1 and 254.
    thresh_position : str, optional
        Option for varying transparency of colormap from upper (high values) or lower (low values), by default 'lower'

    Returns
    -------
    matplotlib.colors.ListedColormap
        Colormap with varying transparency (alpha)
    """

    color_map = get_cmap(cmap, 255)
    cmap_array = asarray(color_map(linspace(0, 1, 255)))

    if thresh < 1:
        thresh = 1
    elif thresh > 254:
        thresh = 254
    
    if thresh_position == 'upper':
        for i in range(thresh):
            cmap_array[254 - i , 3] = float(i) / thresh
    else:
        for i in range(thresh):
            cmap_array[i , 3] = float(i) / thresh

    return ListedColormap(cmap_array)


def cmap_white(cmap):
    """
    Create a colormap where the minimum color value is white

    Parameters
    ----------
    cmap : matplotlib.colors.ListedColormap
        matplotlib colormap. For example: matplotlib.pyplot.cm.inferno

    Returns
    -------
    matplotlib.colors.ListedColormap
        Colormap with minimum value coloured white
    """
    color_map  = get_cmap(cmap , 255)
    cmap_array = asarray(color_map(linspace(0, 1, 255)))

    cmap_array[0 , :] = 0

    return ListedColormap(cmap_array)
