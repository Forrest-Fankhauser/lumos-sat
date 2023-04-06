""" Conversions which are used throughout Lumos """

import numpy as np
import lumos.constants
import lumos.functions

def intensity_to_ab_mag(intensity : float | np.ndarray, clip : bool = True) -> float | np.ndarray:
    """
    Converts intensity to AB Magnitude

    Parameters:
        intensity (np.ndarray or float) : Intensity measured in watts / meter^2
    
    Returns:
        ab_mag (np.ndarray or float) : Brightness in AB Magnitude
    """
    log_val = intensity * lumos.constants.WAVELENGTH / (lumos.constants.SPEED_OF_LIGHT * 3631e-26)
    if clip:
        log_val = np.clip(log_val, 10e-6, None)
    ab_mag = -2.5 * np.log10( log_val )
    return ab_mag

def altaz_to_unit(altitude : float | np.ndarray, azimuth : float | np.ndarray) \
    -> tuple[ float | np.ndarray, ...]:
    """
    Converts altitude and azimuth to a unit vector.

    Parameters:
        altitude (float or np.ndarray) : Altitude in HCS frame (degrees)
        azimuth (float or np.ndarray) : Azimuth in HCS frame (degrees)
    
    Returns:
        x, y, z (tuple) : Unit vector components
    """
    phi = np.deg2rad(90 - altitude)
    theta = np.deg2rad(azimuth)
    x = np.sin(phi) * np.cos(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(phi)
    return x, y, z

def unit_to_spherical(x : float | np.ndarray, y : float | np.ndarray, z : float | np.ndarray) \
    -> tuple[float | np.ndarray]:
    """
    Converts a unit vector to spherical coordinates.

    Parameters:
        x (float or np.ndarray)
        y (float or np.ndarray)
        z (float or np.ndarray)
    
    Returns:
        phi, theta (tuple) : Spherical coordinate representation
    """
    phi = np.arccos(z)
    theta = np.arctan2(y, x)
    phi = np.rad2deg(phi)
    theta = np.rad2deg(theta)
    return phi, theta

def spherical_to_unit(phi : float | np.ndarray, theta : float | np.ndarray) \
    -> tuple[float | np.ndarray]:
    """
    Converts a unit vector to spherical coordinates.

    Parameters:
        phi (float or np.ndarray)
        theta (float or np.ndarray)
    
    Returns:
        x, y, z (tuple) : Unit vector representation
    """
    x = np.sin(phi) * np.cos(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(phi)
    return x, y, z