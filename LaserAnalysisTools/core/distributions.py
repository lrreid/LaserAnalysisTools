"""
Define standard distributions useful for laser and plasma simulations and calculations.
Super-Gauss derivations for conversions between FWHM, 1/e^2 etc. in bound notebook, 12/04/2023.

To Do:
    - Add functions for 2D distributions
"""

import numpy as np

def Gauss_1D_sigma(Amp: float, centroid: float, sigma: float, x_data):
    """
    Calculate a Gaussian distribution for a given sigma (rms)

    Parameters
    ----------
    Amp : float
        Amplitude of the distribution
    centroid : float
        Centroid of the distribution
    sigma : float
         sigma / RMS/ standard deviation of the distribution
    x_data : array_like
        Array containing range of x-values of the distributuion 

    Returns
    -------
    Gauss: array_like
        Array containing Gaussian distribution. Intenstiy values for each value of x_data
    FWHM: float
        Full width half maximum corresponding to given value of sigma of distribution
    """
    Gauss  = Amp * np.exp( -1 * ( ((x_data-centroid)**2)/(2*(sigma**2))  )    )
    FWHM   = (2*np.sqrt(2*np.log(2))*sigma)
    return Gauss, FWHM

def Gauss_1D_FWHM(Amp: float, centroid: float, FWHM: float, x_data):
    """
    Calculate a Gaussian distribution for a given FWHM

    Parameters
    ----------
    Amp : float
        Amplitude of the distribution
    centroid : float
        Centroid of the distribution
    FWHM : float
        Full width half maximum of the distribution
    x_data : array_like
        Array containing range of x-values of the distributuion 

    Returns
    -------
    Gauss: array_like
        Array containing Gaussian distribution. Intenstiy values for each value of x_data
    sigma : float
         sigma / RMS/ standard deviation of the distribution calculated from FWHM
    """
    sigma  = FWHM/(2*np.sqrt(2*np.log(2)))
    Gauss  = Amp * np.exp( (-4*np.log(2)) * ( ((x_data-centroid)**2)/(FWHM**2) )    )
    return Gauss, sigma

def superGauss_1D(Amp: float, centroid: float, sigma: float, order: float, x_data):
    """
    Calculate a Super-Gaussian distribution from a given sigma

    Parameters
    ----------
    Amp : float
        Amplitude of the distribution
    centroid : float
        Centroid of the distribution
    sigma : float
        Sigma of the distribution. Note this is not equivalent to the RMS value.
    order : float
        Order of the super-Gaussian
    x_data : array_like
        Array containing range of x-values of the distributuion 

    Returns
    -------
    Super_Gauss_1D : array_like
        Array containing super-Gaussian distribution. Intenstiy values for each value of x_data
    rms_radius : float
        RMS radius corresponding to sigma of distribution (1/sqrt(e) intensity)
    FWHM : float
        Full width half maximum to sigma of distribution
    e2_radius : float
        1/e^2 radius corresponding to sigma of distribution (1/e^2 intensity)
    """
    Super_Gauss_1D = Amp * np.exp(-1*( ((x_data-centroid)**2) / (2*sigma**2) )**order )

    rms_radius = np.sqrt(2)*sigma*(np.log(np.sqrt(np.exp(1))))**(1/(2*order))
    FWHM       = 2*np.sqrt(2)*sigma*(np.log(1/0.5))**(1/(2*order))
    e2_radius  = np.sqrt(2)*sigma*(np.log(np.exp(1)**2))**(1/(2*order))
    
    return Super_Gauss_1D, rms_radius, FWHM, e2_radius

def superGauss_1D_FWHM(Amp: float, centroid: float, FWHM: float, order: float, x_data):
    """
    Calculate a Super-Gaussian distribution from a given FWHM

    Parameters
    ----------
    Amp : float
        Amplitude of the distribution
    centroid : float
        Centroid of the distribution
    FWHM : float
        Full width half maximum of the distribution
    order : float
        Order of the super-Gaussian
    x_data : array_like
        Array containing range of x-values of the distributuion 

    Returns
    -------
    Super_Gauss_1D : array_like
        Array containing super-Gaussian distribution. Intenstiy values for each value of x_data
    rms_radius : float
        RMS radius corresponding to sigma of distribution (1/sqrt(e) intensity)
    e2_radius : float
        1/e^2 radius corresponding to sigma of distribution (1/e^2 intensity)
    """
    Super_Gauss_1D = Amp * np.exp(-np.log(2)*( ((4*((x_data-centroid)**2)) / (FWHM**2))  )**order )
    rms_radius     = (FWHM/2)*( (np.log( np.sqrt(np.exp(1)))/np.log(2))**(1/(2*order)) )
    e2_radius      = (FWHM/2)*( (np.log( np.exp(1)**2)/np.log(2))**(1/(2*order)) )

    return Super_Gauss_1D, rms_radius, e2_radius

def superGauss_width_conversion_sigma(order: float, sigma: float, norm_int: float):
    """
    Find the radius of a Super-Gaussian distribution from a given sigma at a given intenisty value

    Parameters
    ----------
    order : float
        Order of the super-Gaussian
    sigma : float
        Sigma of the distribution. Note this is not equivalent to the RMS value.
    norm_int : float
        Intensity value (normalised to peak intensity) that you want to find the radius of

    Returns
    -------
    radius : float
        Radius of distribution at normalised intenstiy norm_int
    """

    radius = np.sqrt(2)*sigma*(np.log(1/norm_int))**(1/(2*order))

    return radius

def superGauss_width_conversion_FWHM(order: float, FWHM: float, norm_int: float):
    """
    Find the radius of a Super-Gaussian distribution from a given FWHM at a given intenisty value

    Parameters
    ----------
    order : float
        Order of the super-Gaussian
    FWHM : float
        Full width half maximum of the distribution
    norm_int : float
        Intensity value (normalised to peak intensity) that you want to find the radius of

    Returns
    -------
    radius : float
        Radius of distribution at normalised intenstiy norm_int
    """

    radius = (FWHM/2)*( (np.log( 1/norm_int)/np.log(2))**(1/(2*order)) ) 

    return radius

def superGauss_beam_conversions_FWHM(order: float, FWHM: float, show_prints: bool = False):
    """
    Calculate common distribution width definitions of a Super-Gaussian from a given Full width half max

    Parameters
    ----------
    order : float
        Order of the super-Gaussian
    FWHM : float
        Full width half maximum of the distribution
    show_prints : bool, optional
        Print results to terminal, by default False

    Returns
    -------
    rms_radius : float
        Beam radius at RMS (1/sqrt(e)) intensity
    e2_radius : float
        Beam radius at 1/e^2 of peak intensity contour
    one_pc_contour : flaot
        Beam radius at 1% of peak intensity contour
    clear_aperture : flaot
        Beam radius at 0.1% of peak intensity contour. The "clear aperture" of the distribution
    """

    rms_radius     = (FWHM/2)*( (np.log( np.sqrt(np.exp(1)))/np.log(2))**(1/(2*order)) ) 
    e2_radius      = (FWHM/2)*( (np.log( np.exp(1)**2)/np.log(2))**(1/(2*order)) ) 
    one_pc_contour = (FWHM/2)*( (np.log( 1/0.01  )/np.log(2))**(1/(2*order)) ) 
    clear_aperture = (FWHM/2)*( (np.log( 1/0.001 )/np.log(2))**(1/(2*order)) ) 

    if show_prints is True:
        print("For order %d and FWHM: %d" % (order, FWHM) )
        print("RMS: %0.2f" % np.round(rms_radius, 2) )
        print("1/e2 dia: %0.2f" % np.round(e2_radius*2, 2) )
        print("1 %% intensity contour dia: %0.2f" % np.round(one_pc_contour*2, 2) )
        print("Clear aperture dia: %0.2f" % np.round(clear_aperture*2, 2) )

    return rms_radius, e2_radius, one_pc_contour, clear_aperture


def superGauss_beam_conversions_e2(order: float, e2_radius: float, show_prints: bool = False):
    """
    Calculate common distribution width definitions of a Super-Gaussian from a given 1/e^2 radius

    Parameters
    ----------
    order : float
        Order of the super-Gaussian
    e2_radius : float
        Beam radius at 1/e^2 of peak intensity contour
    show_prints : bool, optional
        Print results to terminal, by default False

    Returns
    -------
    rms_radius : float
        Beam radius at RMS (1/sqrt(e)) intensity
    FWHM : float
        Full width half maximum of the distribution
    one_pc_contour : flaot
        Beam radius at 1% of peak intensity contour
    clear_aperture : flaot
        Beam radius at 0.1% of peak intensity contour. The "clear aperture" of the distribution
    """

    rms_radius     = (e2_radius)*( (np.log(np.sqrt(np.exp(1)))/np.log(np.exp(1)**2))**(1/(2*order)) )
    FWHM           = (2*e2_radius)*( (np.log(2)/np.log(np.exp(1)**2))**(1/(2*order)) )
    one_pc_contour = (e2_radius)*( (np.log(1/0.01 )/np.log(np.exp(1)**2))**(1/(2*order)) )
    clear_aperture = (e2_radius)*( (np.log(1/0.001)/np.log(np.exp(1)**2))**(1/(2*order)) )

    if show_prints is True:
        print("For order %d and 1/e^2 radius: %d" % (order, e2_radius) )
        print("RMS: %0.2f" % np.round(rms_radius, 2) )
        print("FWHM: %0.2f" % np.round(FWHM, 2) )
        print("1 %% intensity contour dia: %0.2f" % np.round(one_pc_contour*2, 2) )
        print("Clear aperture dia: %0.2f" % np.round(clear_aperture*2, 2) )

    return rms_radius, FWHM, one_pc_contour, clear_aperture


def superGauss_beam_conversions_CA(order: float, clear_aperture_rad: float, show_prints: bool = False):
    """
    Calculate common distribution width definitions of a Super-Gaussian from a given clear aperture radius. Defined at the 0.1% intenstiy contour.

    Parameters
    ----------
    order : float
        Order of the super-Gaussian
    clear_aperture : flaot
        Beam radius at 0.1% of peak intensity contour. The "clear aperture" of the distribution
    show_prints : bool, optional
        Print results to terminal, by default False

    Returns
    -------
    rms_radius : float
        Beam radius at RMS (1/sqrt(e)) intensity
    FWHM : float
        Full width half maximum of the distribution
    e2_radius : float
        Beam radius at 1/e^2 of peak intensity contour
    one_pc_contour : flaot
        Beam radius at 1% of peak intensity contour
    """

    rms_radius     = (clear_aperture_rad)*( (np.log(np.sqrt(np.exp(1)))/np.log(1/0.001))**(1/(2*order)) )
    FWHM           = (2*clear_aperture_rad)*( (np.log(2)/np.log(1/0.001))**(1/(2*order)) )
    e2_radius      = (clear_aperture_rad)*( (np.log(np.exp(1)**2)/np.log(1/0.001))**(1/(2*order)) )
    one_pc_contour = (clear_aperture_rad)*( (np.log(1/0.01)/np.log(1/0.001))**(1/(2*order)) )

    if show_prints is True:
        print("For order %d and clear aperture radius: %d" % (order, clear_aperture_rad) )
        print("RMS: %0.2f" % np.round(rms_radius, 2) )
        print("FWHM: %0.2f" % np.round(FWHM, 2) )
        print("1/e2 dia: %0.2f" % np.round(e2_radius*2, 2) )
        print("1 %% intensity contour dia: %0.2f" % np.round(one_pc_contour*2, 2) )

    return rms_radius, FWHM, e2_radius, one_pc_contour