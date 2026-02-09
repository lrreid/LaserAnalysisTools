"""
Simple plasma calculations and conversions
"""

import numpy as np
from scipy.constants import c, e, pi, epsilon_0, m_e
from LaserAnalysisTools.utils.unit_conversions import magnitude_conversion, magnitude_conversion_vol

class PlasmaDen_Conversions():

    def __init__(self, n_e, den_unit: str = "centi", wavelength_unit: str = "micro", freq_unit: str = "Tera", skin_unit: str = "micro", time_unit: str = "femto"):
        """
        Convert plasma density to other useful related quantities

        Parameters
        ----------
        n_e : float
            Plasma electron density. Default unit: cm^-3
        den_unit : str, optional
            Order of magnitude of plasma density unit, by default "centi"
        wavelength_unit : str, optional
            Order of magnitude of plasma wavelength unit, by default "micro"
        freq_unit : str, optional
            Order of magnitude of plasma frequency unit, by default "Tera"
        skin_unit : str, optional
            Order of magnitude of plasma skin depth unit, by default "micro"
        time_unit : str, optional
            Order of magnitude of plasma period unit, by default "femto"
        """

        self.n_e             = n_e
        self.den_unit        = den_unit
        self.wavelength_unit = wavelength_unit
        self.freq_unit       = freq_unit
        self.skin_unit       = skin_unit
        self.time_unit       = time_unit
        self.n_e_SI          = magnitude_conversion_vol(self.n_e, self.den_unit, "", reciprocal_units = True)        # convert plasma density to si units (m^-3)
        self.w_p, self.f_p   = self.plasma_frequency_from_density()                         # Plasma (angular) frequency
        self.lambda_p        = self.plasma_wavelength_from_density()                        # Plasma wavelength
        self.k_p, self.skin_depth = self.plasma_wavevector_from_density()                   # plasma wavevector and skin depth
        self.T_p             = self.plasma_period_from_density()                            # Plasma period

    def plasma_wavelength_from_density(self):
        """
        Calculate plasma wavelength from density

        Returns
        -------
        lambda_p: float
            Plasma wavelength. Default unit: microns
        """

        lambda_p  = (2*pi*c)*np.sqrt( (epsilon_0*m_e)/(e**2*self.n_e_SI) )
        
        return magnitude_conversion(lambda_p, "", self.wavelength_unit)


    def plasma_frequency_from_density(self):
        """
        Calculate plasma linear and angular frequency from density

        Returns
        -------
        w_p: float
            Plasma angular frequency. Unit rad/s
        f_p: float
            Plasma linear frequency. Default unit: THz
        """

        w_p = np.sqrt( (self.n_e_SI*e**2) / (epsilon_0*m_e) )                # Plasma (angular) frequency (rad/s)
        f_p = (1/(2*pi)) * np.sqrt( (self.n_e_SI*e**2) / (epsilon_0*m_e) )   # Plasma (linear) frequency (Hz)

        return w_p, magnitude_conversion(f_p, "", self.freq_unit)


    def plasma_period_from_density(self):
        """
        Calculate plasma frequency from density

        Returns
        -------
        T_p: float
            Plasma period. Default unit: fs
        """

        T_p = ((2*pi)*np.sqrt((epsilon_0*m_e)/(self.n_e_SI*e**2)) )          # Plasma period (s)

        return magnitude_conversion(T_p, "", self.time_unit)


    def plasma_wavevector_from_density(self):
        """
        Calculate plasma wavevector and skin depth from density

        Returns
        -------
        k_p: float
            Plasma wavevector. Unit: 1/m
        skin_depth: float
            Plasma skin depth. Default unit: microns
        """

        k_p = (e/c) * np.sqrt(self.n_e_SI/(epsilon_0*m_e))
        skin_depth = 1/k_p

        return k_p, magnitude_conversion(skin_depth, "", self.skin_unit)

#### end of class ####

def plasma_density_from_frequency(w_p, den_unit: str = "centi"):
    """
    Calculate plasma density from the angular frequency

    Parameters
    ----------
    w_p: float
        Plasma angular frequency. Unit rad/s
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"

    Returns
    -------
    n_e: float
        Plasma electron density. Default unit: cm^-3
    """

    n_e = (w_p**2*epsilon_0*m_e)/(e**2)

    return magnitude_conversion_vol(n_e, "", den_unit, reciprocal_units = True)


def plasma_density_from_wavelength(lambda_p, wavelength_unit: str = "micro", den_unit: str = "centi"):
    """
    Calculate the plasma density from the wavelength

    Parameters
    ----------
    lambda_p : type
        Plasma wavelength. Default unit: microns
    wavelength_unit : str, optional
        Order of magnitude of plasma wavelength unit, by default "micro"
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"

    Returns
    -------
    n_e: float
        Plasma electron density. Default unit: cm^-3
    """

    lambda_p_SI = magnitude_conversion(lambda_p, wavelength_unit, "" ) # convert plasma period to time in SI units (s).
    n_e         = epsilon_0*m_e*((2*pi*c)/(lambda_p_SI*e))**2

    return magnitude_conversion_vol(n_e, "", den_unit, reciprocal_units = True)


def plasma_density_from_period(T_p, time_unit: str = "femto", den_unit: str = "centi"):
    """
    Calculate the plasma density from the period

    Parameters
    ----------
    T_p : _type_
        Plasma period. Default unit: fs
    time_unit : str, optional
        Order of magnitude of plasma period unit, by default "femto"
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"

    Returns
    -------
    n_e: float
        Plasma electron density. Default unit: cm^-3
    """

    T_p_SI = magnitude_conversion(T_p, time_unit, "" ) # convert plasma period to time in SI units (s).
    n_e    = epsilon_0*m_e* ((2*pi)/(e*T_p_SI))**2

    return magnitude_conversion_vol(n_e, "", den_unit, reciprocal_units = True)

def plasma_density_from_wavevector(k_p, den_unit: str = "centi"):
    """
    Calculate the plasma density from the wavevector

    Parameters
    ----------
    k_p : float
        Plasma wavevector. Unit: 1/m
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"

    Returns
    -------
    n_e: float
        Plasma electron density. Default unit: cm^-3
    """

    n_e = ((c*k_p)/e)**2 * epsilon_0*m_e

    return magnitude_conversion_vol(n_e, "", den_unit, reciprocal_units = True)

def plasma_density_from_skin_depth(skin_depth, skin_unit: str = "micro", den_unit: str = "centi"):
    """
    Calculate the plasma density from the skin depth

    Parameters
    ----------
    skin_depth : float
        Plasma skin depth. Default unit: microns
    skin_unit : str, optional
        Order of magnitude of skin depth units, by default "micro"
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"

    Returns
    -------
    n_e: float
        Plasma electron density. Default unit: cm^-3
    """

    skin_depth_SI = magnitude_conversion(skin_depth, skin_unit, "")         # Convert skin depth to SI units (m)
    n_e           = (c/(e*skin_depth_SI))**2 * epsilon_0*m_e

    return magnitude_conversion_vol(n_e, "", den_unit, reciprocal_units = True)


def critical_power(n_e, lambda0: float = 800, den_unit: str = "centi", wavelength_unit: str = "nano", power_unit: str = "Tera"):
    """
    Calculate the critical power for laser self-focusing in a plasma.

    Parameters
    ----------
    n_e: float
        Plasma electron density. Default unit: cm^-3
    lambda0 : float, optional
        Laser central wavelength, by default 800 nm
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"
    wavelength_unit : str, optional
        Order of magnitude of laser wavelength unit, by default "nano"
    power_unit : str, optional
        Orer of magnitude of laser power unit, by default "Tera"

    Returns
    -------
    Pcrit: float
        Critical power for self-focusing. Default unit: TW
    """

    n_e_SI     = magnitude_conversion_vol(n_e, den_unit, "", reciprocal_units = True)
    lambda0_SI = magnitude_conversion(lambda0, wavelength_unit, "")                     # Get laser wavelength in SI unit (m)
    Pcrit      = (32*pi**3*epsilon_0**2*m_e**3*c**7)/(n_e_SI*e**4*lambda0_SI**2)        # critical power in W

    return magnitude_conversion(Pcrit, "", power_unit)


def plasma_density_from_critical_power(Pcrit, lambda0, power_unit: str = "Tera", wavelength_unit: str = "nano", den_unit: str = "centi" ):
    """
    Calculate the threshold density for self-focusing in a plasma from the laser power.

    Parameters
    ----------
    Pcrit: float
        Laser peak power. Default unit: TW
    lambda0 : float, optional
        Laser central wavelength, by default 800 nm
    power_unit : str, optional
        Orer of magnitude of laser power unit, by default "Tera"
    wavelength_unit : str, optional
        Order of magnitude of laser wavelength unit, by default "nano"
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"

    Returns
    -------
    n_e: float
        Plasma electron density where critical power begins to occur. Default unit: cm^-3
    """

    Pcrit_SI   = magnitude_conversion(Pcrit, power_unit, "")                            # Get critical power in SI unit (Watt)
    lambda0_SI = magnitude_conversion(lambda0, wavelength_unit, "")                     # Get laser wavelength in SI unit (m)
    n_e        = (32*pi**3*epsilon_0**2*m_e**3*c**7)/(Pcrit_SI*e**4*lambda0_SI**2)      # Plasma density at critical power (m^-3)

    return magnitude_conversion_vol(n_e, "", den_unit, reciprocal_units = True)


def critical_density(lambda0, wavelength_unit: str = "nano", den_unit: str = "centi" ):
    """
    Calculate the critical density where the plasma transitions from underdense to overdense from the laser wavelength

    Parameters
    ----------
    lambda0 : float, optional
        Laser central wavelength, by default 800 nm
    wavelength_unit : str, optional
        Order of magnitude of laser wavelength unit, by default "nano"
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"

    Returns
    -------
    n_c: float
        Plasma critical density. Default unit: cm^-3
    """

    lambda0_SI = magnitude_conversion(lambda0, wavelength_unit, "")         # Get laser wavelength in SI unit (m)
    n_c        = m_e*epsilon_0*((2*pi*c)/(lambda0_SI*e))**2                 # Critical density (m^-3)

    return magnitude_conversion_vol(n_c, "", den_unit, reciprocal_units = True)


def critical_density_from_frequency(omega_l, den_unit: str = "centi"):
    """
    Calculate the critical density where the plasma transitions from underdense to overdense from the laser angular frequency

    Parameters
    ----------
    omega_l : float
        Laser angular frequency (rad/s)
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"

    Returns
    -------
    n_c: float
        Plasma critical density. Default unit: cm^-3
    """

    n_c = (m_e*epsilon_0*omega_l**2)/(e**2)                                 # Critical density (m^-3)

    return magnitude_conversion_vol(n_c, "", den_unit, reciprocal_units = True)


def plasma_refractive_index(lambda0, n_e, wavelength_unit: str = "nano", den_unit: str = "centi"):
    """
    Calculate the plasma refractive index from the laser wavelength and plasma density

    Parameters
    ----------
    lambda0 : float, optional
        Laser central wavelength, by default 800 nm
    n_e: float
        Plasma electron density. Default unit: cm^-3
    wavelength_unit : str, optional
        Order of magnitude of laser wavelength unit, by default "nano"
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"

    Returns
    -------
    eta: float
        Plasma refractive index
    """

    lambda0_SI = magnitude_conversion(lambda0, wavelength_unit, "")                         # Get laser wavelength in SI unit (m)
    n_e_SI     = magnitude_conversion_vol(n_e, den_unit, "", reciprocal_units = True)       # Plasma density in SI unit (m^-3)
    eta        = np.sqrt(1 - ( (n_e_SI/(4*epsilon_0*m_e)) * ((e*lambda0_SI)/(pi*c))**2 ) )

    return eta

def laser_group_velocity(lambda0, n_e, wavelength_unit: str = "nano", den_unit: str = "centi"):
    """
    Calculate the laser group velocity in plasma from the laser wavelength and plasma density

    Parameters
    ----------
    lambda0 : float, optional
        Laser central wavelength, by default 800 nm
    n_e: float
        Plasma electron density. Default unit: cm^-3
    wavelength_unit : str, optional
        Order of magnitude of laser wavelength unit, by default "nano"
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"

    Returns
    -------
    v_group: float
        Laser group velocity in plasma. Unit: m/s
    """
    
    v_group = c * plasma_refractive_index(lambda0, n_e, wavelength_unit, den_unit)

    return v_group


def wavebreaking_field(n_e, den_unit: str = "centi", field_unit: str = "Giga"):
    """
    Calculate the plasma wavebreaking field from the plasma density

    Parameters
    ----------
    n_e: float
        Plasma electron density. Default unit: cm^-3
    den_unit : str, optional
        Order of magnitude of plasma density unit, by default "centi"
    field_unit : str, optional
        Order of magnitude of electric field unit, by default "Giga"

    Returns
    -------
    E_wb = float
        Plasma wavebreaking field. Default unit: GV/m
    """

    n_e_SI = magnitude_conversion_vol(n_e, den_unit, "", reciprocal_units = True)
    E_wb   = c * np.sqrt((m_e*n_e_SI)/epsilon_0)

    return magnitude_conversion(E_wb, "", field_unit) 
 