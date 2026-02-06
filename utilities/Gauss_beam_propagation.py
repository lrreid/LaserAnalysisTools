"""
Functions related to propagation of Gaussian laser beams
"""


import numpy as np
from scipy.constants import pi
from LaserAnalysisTools.core.unit_conversions import magnitude_conversion


def rayleigh_length(w0, lambda0, spot_unit: str = "micro", wavelength_unit: str = "nano", rayleigh_unit: str = "milli"):
    """
    Calculate the Rayleigh length associated with a known laser focus.

    Parameters
    ----------
    w0 : float
        Laser spot size (radius to 1/e^2 of peak intensity) at focus. Default unit: microns
    lambda0 : float
        Laser central wavelength. Default unit: nanometers
    spot_unit : str, optional
        Order of magnitude of laser spot size, by default "micro"
    wavelength_unit : str, optional
        Order of magnitude of wavelength unit, by default "nano"
    rayleigh_unit : str, optional
        Order of magnitude of Rayleigh length unit, by default "milli"

    Returns
    -------
    Zr : float
        Rayleigh length of laser focus. Default unit: mm
    """

    w0_SI       = magnitude_conversion(w0, spot_unit, "")
    lambda0_SI  = magnitude_conversion(lambda0, wavelength_unit, "")

    zr = (pi*(w0_SI**2))/lambda0_SI    

    return magnitude_conversion(zr, "", rayleigh_unit)


def Gaussian_propagation(w0, lambda0, z, spot_unit: str = "micro", wavelength_unit: str = "nano", z_unit: str = "milli"):
    """
    Calculate the transverse size of a Gaussian laser beam as it propagates from its focus.

    Parameters
    ----------
    w0 : float
        Laser spot size (radius to 1/e^2 of peak intensity) at focus. Default unit: microns
    lambda0 : float
        Laser central wavelength. Default unit: nanometers
    z : float or (numpy) array of floats
        Distance relative to the focal plane of the laser. Default unit: mm
    spot_unit : str, optional
        Order of magnitude of laser spot size, by default "micro"
    wavelength_unit : str, optional
        Order of magnitude of wavelength unit, by default "nano"
    z_unit : str, optional
        Order of magnitude of propagation distance unit, by default "milli"

    Returns
    -------
    w: float
        Transverse size of laser at distance z from the focus. Default unit: microns
    """

    w0_SI       = magnitude_conversion(w0, spot_unit, "")
    lambda0_SI  = magnitude_conversion(lambda0, wavelength_unit, "")
    z_SI        = magnitude_conversion(z, z_unit, "")

    w = w0_SI*np.sqrt(1+((lambda0_SI*z_SI)/(pi*w0_SI**2))**2)

    return magnitude_conversion(w, "", spot_unit)


def prop_for_beam_size(w0, lambda0, beam_radius, spot_unit: str = "micro", wavelength_unit: str = "nano", beam_unit: str = "micro", z_unit: str = "milli"):
    """
    Calculate the distance from the waist of a laser to reach a given transverse beam size.

    Parameters
    ----------
    w0 : float
        Laser spot size (radius to 1/e^2 of peak intensity) at focus. Default unit from spot_unit: microns
    lambda0 : float
        Laser central wavelength. Default unit: nanometers
    beam_radius : float
        Size of laser spot away from focus. Default unit from beam_unit: microns
    spot_unit : str, optional
        Order of magnitude of laser spot size, by default "micro"
    wavelength_unit : str, optional
        Order of magnitude of wavelength unit, by default "nano"
    beam_unit : str, optional
        Order of magnitude of beam size, by default "micro"
    z_unit : str, optional
        Order of magnitude of propagation distance unit, by default from z_unit "milli"

    Returns
    -------
    Z: float.
        Distance from waist to achieve chosen beam size. Default unit: microns
    """

    w0_SI          = magnitude_conversion(w0, spot_unit, "")
    lambda0_SI     = magnitude_conversion(lambda0, wavelength_unit, "")
    beam_radius_SI = magnitude_conversion(beam_radius, beam_unit, "")

    if w0 == beam_radius:
        Z = 0
    else:
        Z = ((pi*w0_SI**2)/lambda0_SI) * np.sqrt((beam_radius_SI/w0_SI)**2-1)
    
    return magnitude_conversion(Z, "", z_unit)


def waist_from_focal_length(focal_length, lambda0, beam_radius, f_unit: str = "milli", wavelength_unit: str = "nano", beam_unit: str = "micro", spot_unit: str = "micro", rayleigh_unit: str = "milli"):
    """
    Caclulate the size of the laser focal spot of a Gaussian beam for a given focal length optic and beam size.

    Parameters
    ----------
    focal_length : float
        Focal length of focusing optic (mirror/OAP). Default unit from f_unit: mm
    lambda0 : float
        Laser central wavelength. Default unit from wavelength unit: nanometers
    beam_radius : float
        Radius of laser pulse on focusing optic. Default unit from beam_unit: microns
    f_unit : str, optional
        Order of magnitude of the focal length unit, by default "milli"
    wavelength_unit : str, optional
        Order of magnitude of wavelength unit, by default "nano"
    beam_unit : str, optional
        Order of magnitude of the beam radius unit, by default "micro"
    spot_unit : str, optional
        Order of magnitude of laser spot size, by default "micro"
    rayleigh_unit : str, optional
        Order of magnitude of Rayleigh length unit, by default "milli"

    Returns
    -------
    w0 : float
        Laser spot size (radius to 1/e^2 of peak intensity) at focus. Default unit from spot_unit: microns
    zr : float
        Rayleigh length of laser focus. Default unit from rayleigh_unit: mm
    """

    focal_length_SI = magnitude_conversion(focal_length, f_unit, "")
    lambda0_SI      = magnitude_conversion(lambda0, wavelength_unit, "")
    beam_radius_SI  = magnitude_conversion(beam_radius, beam_unit, "")

    w0 = (focal_length_SI*lambda0_SI)/(pi*beam_radius_SI)
    zr = rayleigh_length(w0, lambda0_SI, spot_unit = "", wavelength_unit = "", rayleigh_unit = rayleigh_unit)

    return magnitude_conversion(w0, "", spot_unit), zr


def focal_length_from_waist(beam_radius, w0, lambda0, beam_unit: str = "micro", spot_unit: str = "micro", wavelength_unit: str = "nano", f_unit: str = "milli"):
    """
    Calculate the focal length of optic required to give spot size w0 from known beam size on the focusing optic for a Gaussian beam.

    Parameters
    ----------
    beam_radius : float
        Radius of laser pulse on focusing optic. Default unit from beam_unit: microns
    w0 : float
        Laser spot size (radius to 1/e^2 of peak intensity) at focus. Default unit from spot_unit: microns
    lambda0 : float
        Laser central wavelength. Default unit from wavelength unit: nanometers
    beam_unit : str, optional
        Order of magnitude of the beam radius unit, by default "micro"
    spot_unit : str, optional
        Order of magnitude of laser spot size, by default "micro"
    wavelength_unit : str, optional
        Order of magnitude of wavelength unit, by default "nano"
    f_unit : str, optional
        Order of magnitude of the focal length unit, by default "milli"

    Returns
    -------
    focal_length : float
        Focal length of focusing optic (mirror/OAP). Default unit from f_unit: mm
    """
    
    beam_radius_SI = magnitude_conversion(beam_radius, beam_unit, "")
    w0_SI          = magnitude_conversion(w0, spot_unit, "")
    lambda0_SI     = magnitude_conversion(lambda0, wavelength_unit, "")

    focal_length = (pi*beam_radius_SI*w0_SI)/lambda0_SI

    return magnitude_conversion(focal_length, "", f_unit)