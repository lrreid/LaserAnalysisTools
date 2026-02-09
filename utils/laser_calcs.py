"""

Simple laser calculations and conversions

To Do:
    - Add function for focusing of top hat beams
    - Tidy up all the functions
    - Add doc strings


"""

import numpy as np
from scipy.constants import c, e, pi, epsilon_0, m_e, h
from LaserAnalysisTools.core.unit_conversions import magnitude_conversion, magnitude_conversion_area # type: ignore
from LaserAnalysisTools.core.basic_calcs import area_circle, area_circle_projection

class laser_wavelength_conversions():

    def __init__(self, lambda0: float, wavelength_unit: str = "nano", freq_unit: str = "Tera", wavenumber_unit: str = "centi", time_unit: str = "femto", eV_unit: str = "", energy_unit: str = "" ):
        """
        Convert laser wavelength to other useful related quantities

        Parameters
        ----------
        lambda0 : float
            Central wavelength of laser radiation. Default unit: nanometers
        wavelength_unit : str, optional
            Order of magnitude of wavelength unit, by default "nano"
        freq_unit : str, optional
            Order of magnitude of frequency unit, by default "Tera"
        wavenumber_unit : str, optional
            Order of magnitude of wavenumber unit, by default "centi"
        time_unit : str, optional
            Order of magnitude of time units, by default "femto"
        eV_unit : str, optional
            Order of magnitude of energy in electronvolts units, by default ""
        energy_unit : str, optional
            Order of magnitude of energy in SI units, by default ""
        """

        self.lambda0          = lambda0
        self.wavelength_unit  = wavelength_unit
        self.freq_unit        = freq_unit
        self.wavenumber_unit  = wavenumber_unit
        self.time_unit        = time_unit
        self.eV_unit          = eV_unit
        self.energy_unit      = energy_unit
        self.lambda0_SI       = magnitude_conversion(self.lambda0, self.wavelength_unit, "", reciprocal_units = False)
        self.freq, self.omega = self.get_frequency()
        self.wavenumber       = self.get_wavenumber()
        self.period           = self.get_period()
        self.Eph              = self.get_photon_energy()
    

    def get_frequency(self):
        """
        Get the frequency and angular frequency of the laser from the wavelength

        Returns
        -------
        freq: float
            Frequency of the laser radiation. Default unit: THz
        omega: float
            Angular frequency of the laser radiation. Default unit: rad/s 
        """

        freq  = c/self.lambda0_SI
        omega = (2*pi*c)/self.lambda0_SI

        return magnitude_conversion(freq, "", self.freq_unit), omega
    
    def get_wavenumber(self):
        """
        Get wavenumber of the laser from the central wavelength

        Returns
        -------
        k float
            wavenumber of laser. Default unit: 1/m
        """

        k = 1/self.lambda0_SI

        return magnitude_conversion(k, "", self.wavenumber_unit, reciprocal_units = True)
    

    def get_period(self):
        """
        Get the period of the laser radiation from the wavelength

        Returns
        -------
        T: float
            Laser period. Default unit: fs
        """

        T = self.lambda0_SI/c

        return magnitude_conversion(T, "", self.time_unit)
    
    def get_photon_energy(self):
        """
        Get the energy of a photon from the laser wavelength

        Returns
        -------
        Eph_eV: float
            Photon energy in electronvolts. Default unit: eV
        Eph_J: float
            Photon energy in Joules: Default unit: J
        """

        Eph_eV = (h*c)/(self.lambda0*e)
        Eph_J  = (h*c)/(self.lambda0)

        return magnitude_conversion(Eph_eV, "", self.eV_unit), magnitude_conversion(Eph_J, "", self.energy_unit)

#%% Focusing of top hat beams

def a0_from_intensity(Int, lambda0: float = 800, int_unit: str = "centi", wavelength_unit: str = "nano"):
    """
    Get the laser normalised vector potential from the peak intensity

    Parameters
    ----------
    Int : float
        Laser peak intensity: Default unit: Wcm^-2
    lambda0 : int, optional
        Laser central wavelength, by default 800 nm
    int_unit : str, optional
        Order of magnitude of intensity units, by default "centi"
    wavelength_unit : str, optional
        Order of magnitude of wavelength units, by default "nano"

    Returns
    -------
    a0: float
        Laser normalised vector potential.
    """

    lambda0_SI = magnitude_conversion(lambda0, wavelength_unit, "")
    Int_SI     = magnitude_conversion_area(Int, int_unit, "", reciprocal_units = True)
    a0         = ((e*lambda0_SI)/(pi*m_e)) * np.sqrt((Int_SI)/(2*epsilon_0*c**5))           # laser a0
    
    return a0


def intensity_from_a0(a0, lambda0: float = 800, wavelength_unit: str = "nano", int_unit: str = "centi"):
    """
    Get laser peak intensity from the normalised vector potential

    Parameters
    ----------
    a0 : float
        normalised vector potential
    lambda0 : int, optional
        Laser central wavelength, by default 800 nm
    wavelength_unit : str, optional
        Order of magnitude of wavelength units, by default "nano"
    int_unit : str, optional
        Order of magnitude of intensity units, by default "centi"

    Returns
    -------
    Int: float
        Laser peak intensity. Default unit: Wcm^-2
    """

    lambda0_SI = magnitude_conversion(lambda0, wavelength_unit, "")

    Int = (2*a0**2*pi**2*epsilon_0*m_e**2*c**5)/(e**2*lambda0_SI**2)
    
    return magnitude_conversion_area(Int, "", int_unit, reciprocal_units = True)


def Gaussian_laser_intensitiy(Energy, tau_FWHM, w0, lambda0: float = 800, energy_unit: str = "milli", time_unit: str = "femto", spot_unit: str = "micro", wavelength_unit: str = "nano", int_unit: str = "centi"):
    """
    Calculate peak intensity of a Gaussian beam from experimentally measurable values.

    Parameters
    ----------
    Energy : float
        Laser energy within focus. Default unit: mJ
    tau_FWHM : float
        Pulse duration to FWHM. Default unit: fs
    w0 : float
        Laser spot size to 1/e^2 radius. Default unit: microns
    lambda0 : int, optional
        Laser central wavelength, by default 800 nm
    energy_unit : str, optional
        Order of magnitude of energy unit, by default "milli"
    time_unit : str, optional
        Order of magnitude of time unit, by default "femto"
    spot_unit : str, optional
        Order of magnitude of spot size unit, by default "micro"
    wavelength_unit : str, optional
        Order of magnitude of wavelength units, by default "nano"
    int_unit : str, optional
        Order of magnitude of intensity units, by default "centi"

    Returns
    -------
    Int: float
        Laser peak intensity. Default unit: Wcm^-2
    a0: float
        Normalised vector potential
    """

    Energy_SI   = magnitude_conversion(Energy, energy_unit, "")
    tau_FWHM_SI = magnitude_conversion(tau_FWHM, time_unit, "")
    w0_SI       = magnitude_conversion(w0, spot_unit, "")
    lambda0_SI  = magnitude_conversion(lambda0, wavelength_unit, "")

    Int = (2*Energy_SI)/(tau_FWHM_SI*pi*w0_SI**2)                                                # Intensity (W/m^2). Gaussian profile assumed.
    a0  = ((e*lambda0_SI)/(pi*m_e*c**2*w0_SI)) * np.sqrt( Energy_SI/(pi*epsilon_0*c*tau_FWHM_SI) )  # Laser amplitude. Gaussian profile assumed.

    return magnitude_conversion_area(Int, "", int_unit, reciprocal_units = True), a0


def Gaussian_spot_size_from_intensity(Energy, tau_FWHM, Int, energy_unit: str = "milli", time_unit: str = "femto", int_unit: str = "centi", spot_unit: str = "micro" ):
    
    Energy_SI   = magnitude_conversion(Energy, energy_unit, "")
    tau_FWHM_SI = magnitude_conversion(tau_FWHM, time_unit, "") 
    Int_SI      = magnitude_conversion_area(Int, int_unit, "", reciprocal_units = True)
    
    spot_size = (2*Energy_SI)/(tau_FWHM_SI*pi*Int_SI)
    
    return magnitude_conversion(spot_size, "", spot_unit)


def Gaussian_Energy_from_intenstiy(tau_FWHM, w0, Int, time_unit: str = "femto", spot_unit: str = "micro", int_unit: str = "centi", energy_unit: str = "milli"):

    tau_FWHM_SI = magnitude_conversion(tau_FWHM, time_unit, "")
    w0_SI       = magnitude_conversion(w0, spot_unit, "")
    Int_SI      = magnitude_conversion_area(Int, int_unit, "", reciprocal_units = True)
    
    Energy = 0.5*Int_SI*tau_FWHM_SI*pi*(w0_SI**2)
    
    return magnitude_conversion(Energy, "", energy_unit)


def get_laser_power(Energy, tau_FWHM, energy_unit: str = "milli", time_unit: str = "femto", power_unit: str = "Tera"):

    Energy_SI   = magnitude_conversion(Energy, energy_unit, "")
    tau_FWHM_SI = magnitude_conversion(tau_FWHM, time_unit, "") 
    
    Power = Energy_SI/tau_FWHM_SI
    
    return magnitude_conversion(Power, "", power_unit)

#%% Focusing of top hat beams

def top_hat_laser_intensitiy(Energy, tau_FWHM, beam_rad, lambda0, focal_length, energy_unit: str = "milli", time_unit: str = "femto", beam_unit: str = "milli", f_unit: str = "milli", wavelength_unit: str = "nano", int_unit: str = "centi", spot_unit: str = "micro"):
    """
    Calculate peak intensity of a Top hat beam from experimentally measurable values using the beam size on the focusing optic.

    Parameters
    ----------
    Energy : float
        Laser energy within focus. Default unit: mJ
    tau_FWHM : float
        Pulse duration to FWHM. Default unit: fs
    beam_rad : float
        Radius of laser beam on focusing optic. Default unit: mm
    lambda0 : int, optional
        Laser central wavelength, by default 800 nm
    focal_length : float
        Focal length of focusing optic. Default unit: mm
    energy_unit : str, optional
        Order of magnitude of energy unit, by default "milli"
    time_unit : str, optional
        Order of magnitude of time unit, by default "femto"
    beam_unit : str, optional
         Order of magnitude of laser transverse size before focusing, by default "milli"
    f_unit : str, optional
        Order of magnitude of focal length unit, by default "milli"
    wavelength_unit : str, optional
        Order of magnitude of wavelength units, by default "nano"
    int_unit : str, optional
        Order of magnitude of intensity units, by default "centi"
    spot_unit : str, optional
        Order of magnitude of spot size unit, by default "micro"

    Returns
    -------
    Int: float
        Laser peak intensity. Default unit: Wcm^-2
    a0: float
        Normalised vector potential
    w0: float
        Laser spot size at focus. Default unit: microns
    """

    Energy_SI       = magnitude_conversion(Energy, energy_unit, "")
    tau_FWHM_SI     = magnitude_conversion(tau_FWHM, time_unit, "")
    beam_rad_SI     = magnitude_conversion(beam_rad, beam_unit, "")
    lambda0_SI      = magnitude_conversion(lambda0, wavelength_unit, "")
    focal_length_SI = magnitude_conversion(focal_length, f_unit, "")
    
    omega_0 = 0.823*lambda0_SI*(focal_length_SI/(2*beam_rad_SI))
    Int     = (Energy_SI*pi*beam_rad_SI**2)/(tau_FWHM_SI*(lambda0_SI**2)*(focal_length_SI**2))    
    a0      = a0_from_intensity(Int, lambda0=lambda0_SI, int_unit="", wavelength_unit="")
    
    return magnitude_conversion_area(Int, "", int_unit, reciprocal_units = True), a0, magnitude_conversion(omega_0, "", spot_unit)

def top_hat_laser_intensity_from_spot(Energy, tau_FWHM, lambda0, omega_0, energy_unit: str = "milli", time_unit: str = "femto", wavelength_unit: str = "nano", spot_unit: str = "micro", int_unit: str = "centi"):
    """
    Calculate peak intensity of a Top hat beam from experimentally measurable values using the spot size at focus.

    Parameters
    ----------
    Energy : float
        Laser energy within focus. Default unit: mJ
    tau_FWHM : float
        Pulse duration to FWHM. Default unit: fs
    lambda0 : int, optional
        Laser central wavelength, by default 800 nm
    omega_0 : float
        Laser spot size to 1/e^2 radius. Default unit: microns
    energy_unit : str, optional
        Order of magnitude of energy unit, by default "milli"
    time_unit : str, optional
        Order of magnitude of time unit, by default "femto"
    wavelength_unit : str, optional
        Order of magnitude of wavelength units, by default "nano"
    spot_unit : str, optional
        Order of magnitude of spot size unit, by default "micro"
    int_unit : str, optional
        Order of magnitude of intensity units, by default "centi"

    Returns
    -------
    Int: float
        Laser peak intensity. Default unit: Wcm^-2
    a0: float
        Normalised vector potential
    """

    Energy_SI   = magnitude_conversion(Energy, energy_unit, "")
    tau_FWHM_SI = magnitude_conversion(tau_FWHM, time_unit, "")
    lambda0_SI  = magnitude_conversion(lambda0, wavelength_unit, "")
    omega_0_SI  = magnitude_conversion(omega_0, spot_unit, "")
    
    Int = ( (0.823**2)*Energy_SI*pi ) / ( 4*tau_FWHM_SI*(omega_0_SI)**2  )
    a0  = a0_from_intensity(Int, lambda0=lambda0_SI, int_unit="", wavelength_unit="")
    
    return magnitude_conversion_area(Int, "", int_unit, reciprocal_units = True), a0

#%% Other calculations

def photon_energy_from_wavelength(lambda0: float = 800, wavelength_unit: str = "nano", eV_unit: str = "", energy_unit: str = ""):
    """
    Get the energy of a photon from the laser wavelength.
    Should I use the laser_wavelength_conversions() class for this? Seems silly to do lots of pointless calculations when I only need one of them.

    Parameters
    ----------
    lambda0 : float, optional
        Laser central wavelength, by default 800 nm
    wavelength_unit : str, optional
        Order of magnitude of wavelength unit, by default "nano"
    eV_unit : str, optional
        Order of magnitude of energy in electronvolts units, by default ""
    energy_unit : str, optional
        Order of magnitude of energy in SI units, by default ""

    Returns
    -------
    Eph_eV: float
        Photon energy in electronvolts. Default unit: eV
    Eph_J: float
        Photon energy in Joules: Default unit: J
    """

    lambda0_SI = magnitude_conversion(lambda0, wavelength_unit, "")
    Eph_eV     = (h*c)/(lambda0_SI*e)
    Eph_J      = (h*c)/(lambda0_SI)

    return magnitude_conversion(Eph_eV, "", eV_unit), magnitude_conversion(Eph_J, "", energy_unit)


def wavelength_from_photon_energy(Eph: float, energy_unit_type: str = "eV", energy_unit: str = "", wavelength_unit: str = "nano"):
    """
    Calculate the laser wavelength from the photon energy.

    Parameters
    ----------
    Eph : float
        Photon energy
    energy_unit_type : str, optional
        Type of units of photon energy. Choose either "Joule" or "eV", by default "eV"
    energy_unit : str, optional
        Order of magnitude of energy unit, by default ""
    wavelength_unit : str, optional
        Order of magnitude of wavelength unit, by default "nano"

    Returns
    -------
    lambda0: float
        Laser central wavelength. Default units: nm
    """

    if energy_unit_type == "Joule":
        lambda0 = (h*c) / magnitude_conversion(Eph, energy_unit, "")
    elif energy_unit_type == "eV":
        lambda0 = (h*c) / (e*magnitude_conversion(Eph, energy_unit, ""))

    return magnitude_conversion(lambda0, "", wavelength_unit)

def ponderomotive_energy(a0, energy_unit: str = "kilo"):
    """
    Calculate the energy of an electron (in electronvolts) after being given a ponderomotive kick by a laser field.
    Taken from: https://arxiv.org/pdf/2207.05773 and https://arxiv.org/pdf/2408.09052v1

    Parameters
    ----------
    a0 : float
        Normalised laser vector potential
    energy_unit : str, optional
        Order of magnitude of energy unit, by default "kilo"

    Returns
    -------
    E_pond_max: float
        Maximum kinetic energy of electrons accelerated in the ponderomotive potential. Default unit: keV
    E_pond_avg: float
        Cycle-averaged ponderomotive electron energy. Default unit: keV
    """

    E_pond_max = ((m_e*c**2)/e)*(np.sqrt(1+(a0**2))-1)
    E_pond_avg = ((m_e*c**2)/e)*(np.sqrt(1+((a0**2)/2))-1)

    return magnitude_conversion(E_pond_max, "", energy_unit), magnitude_conversion(E_pond_avg, "", energy_unit)

def photons_per_pulse(Energy, lambda0: float = 800, energy_unit: str = "milli", wavelength_unit: str = "nano"):
    """
    Calculate the number of photons per laser pulse

    Parameters
    ----------
    Energy : float
        Energy within a laser pulse. Default unit: mJ
    lambda0 : float, optional
        Laser central wavelength, by default 800 nm
    energy_unit : str, optional
        Order of magnitude of energy unit, by default "milli"
    wavelength_unit : str, optional
        Order of magnitude of wavelength unit, by default "nano"

    Returns
    -------
    N_photons: float
        Number of photons per pulse
    """

    E_photon  = photon_energy_from_wavelength(lambda0, wavelength_unit, eV_unit = "", energy_unit = "")[1]
    Energy_SI = magnitude_conversion(Energy, energy_unit, "")
    N_photons = Energy_SI/E_photon

    return N_photons


def ND_filter_transmission(OD):
    """
    Calculate the transmission through a neutral density filter

    Parameters
    ----------
    OD : float
        Optical density of filter

    Returns
    -------
    Transmission: float
        Transmission through filter
    """

    Transmission = 10**(-OD)

    return Transmission

def ND_filter_Energy_transmission(OD, Energy_in, energy_in_unit: str = "milli", energy_out_unit: str = "micro"):
    """
    Calculate the laser energy passing through or reflected by a neutral density filter.

    Parameters
    ----------
    OD : float
        Optical density of filter
    Energy_in : float
        Energy incident on neutral density filter
    energy_in_unit : str, optional
        Order of magnitude of energy units of incident beam, by default "milli"
    energy_out_unit : str, optional
        Order of magnitude of energy units of beam after filter, by default "micro"

    Returns
    -------
    Energy_out: float
        Energy of laser after neutral density filter
    """

    Energy_in_SI = magnitude_conversion(Energy_in, energy_in_unit, "")
    Energy_out   = Energy_in_SI*ND_filter_transmission(OD)

    return magnitude_conversion(Energy_out, "", energy_out_unit)



#%% Calculate laser beam energy and power denisty

def energy_density(Laser_Energy, Beam_diameter, Beam_profile: str = "Super-Gauss", energy_unit: str = "milli", dia_unit: str = "milli", Eden_energy_unit: str = "milli", Eden_area_unit: str = "centi" ):
    """
    For a Super-Gaussian beam, use the FWHM for the diameter.
    For a Gaussian beam, use the 1/e^2 for the diameter.

    """

    Laser_Energy_SI = magnitude_conversion(Laser_Energy, energy_unit, "", reciprocal_units = False)
    
    Energy_Density = Laser_Energy_SI/area_circle(Beam_diameter, dia_unit = dia_unit, area_unit=Eden_area_unit)

    if Beam_profile == "Super-Gauss":
        Energy_Density = Laser_Energy_SI/area_circle(Beam_diameter, dia_unit = dia_unit, area_unit=Eden_area_unit)
    elif Beam_profile == "Gaussian":
        Energy_Density = (2*Laser_Energy_SI)/area_circle(Beam_diameter, dia_unit = dia_unit, area_unit=Eden_area_unit)
    else:
        print("Beam profile incorrectly defined. Choose Super-Gauss or Gaussian.\nDefaulted to Super-Gauss.")
        Energy_Density = Laser_Energy_SI/area_circle(Beam_diameter, dia_unit = dia_unit, area_unit=Eden_area_unit)
        

    return magnitude_conversion(Energy_Density, "", Eden_energy_unit, reciprocal_units = False)


def peak_power_density(Laser_Energy, Beam_diameter, tau_FWHM, Beam_profile: str = "Super-Gauss", energy_unit: str = "milli", time_unit: str = "femto", dia_unit: str = "milli", Pden_power_unit: str = "Mega", Pden_area_unit: str = "centi"):

    Energy_density = energy_density(Laser_Energy, Beam_diameter, Beam_profile = Beam_profile, energy_unit = energy_unit, dia_unit = dia_unit, Eden_energy_unit = "", Eden_area_unit = Pden_area_unit )

    tau_FWHM_SI   = magnitude_conversion(tau_FWHM, time_unit, "", reciprocal_units = False)
    Power_density = Energy_density/tau_FWHM_SI

    return magnitude_conversion(Power_density, "", Pden_power_unit, reciprocal_units = False) 

def average_power_denstiy(Laser_Energy, Beam_diameter, rep_rate, Beam_profile: str = "Super-Gauss", energy_unit: str = "milli", rep_unit: str = "", dia_unit: str = "milli", Pden_power_unit: str = "", Pden_area_unit: str = "centi"):

    rep_rate_Hz   = magnitude_conversion(rep_rate, rep_unit, "", reciprocal_units = False)

    Energy_density = energy_density(Laser_Energy, Beam_diameter, Beam_profile = Beam_profile, energy_unit = energy_unit, dia_unit = dia_unit, Eden_energy_unit = "", Eden_area_unit = Pden_area_unit )
    Power_density  = Energy_density * rep_rate_Hz

    return magnitude_conversion(Power_density, "", Pden_power_unit, reciprocal_units = False)



#%% TO DO: Add similar laser beam energy and power denisty calculations for a beam projected at an angle