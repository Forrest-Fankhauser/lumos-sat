""" Conversions which are used throughout Lumos """

import numpy as np
import lumos.constants
import lumos.functions

def intensity_to_ab_mag(intensity, clip = True):
    """
    Converts from intensity to AB Magnitude.
    If clip is set to True, outputs below 12 AB Magnitude will be clipped.

    :param intensity: Intensity :math:`\\frac{W}{m^2}`
    :type intensity: :class:`np.ndarray` or float
    :param clip: Whether or not to clip very small intensities
    :type clip: bool, optional
    :return: AB Magnitude
    :rtype: :class:`np.ndarray` or float
    """
    log_val = intensity * lumos.constants.WAVELENGTH / (lumos.constants.SPEED_OF_LIGHT * 3631e-26)
    if clip:
        log_val = np.clip(log_val, 10e-6, None)
    ab_mag = -2.5 * np.log10( log_val )
    return ab_mag

def altaz_to_unit(altitude, azimuth):
    """
    Converts altitude and azimuth to a unit vector.

    :param altitude: Altitude in HCS frame (degrees)
    :type altitude: :class:`np.ndarray` or float
    :param azimuth: Azimuth in HCS frame (degrees)
    :type azimuth: :class:`np.ndarray` or float
    :return: Unit vector components :math:`(x, y, z)`
    :rtype: tuple
    """
    phi = np.deg2rad(90 - altitude)
    theta = np.deg2rad(azimuth)
    x = np.sin(phi) * np.cos(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(phi)
    return x, y, z

def unit_to_spherical(x, y, z):
    """
    Converts a unit vector :math:`(x, y, z)` to spherical coordinates
    :math:`(\\phi, \\theta)`

    :param x:
    :type x: :class:`np.ndarray` or float
    :param y:
    :type y: :class:`np.ndarray` or float
    :param z:
    :type z: :class:`np.ndarray` or float
    :return: :math:`(\\phi, \\theta)`
    :rtype: tuple
    """
    phi = np.arccos(z)
    theta = np.arctan2(y, x)
    phi = np.rad2deg(phi)
    theta = np.rad2deg(theta)
    return phi, theta

def spherical_to_unit(phi, theta):
    """
    Converts from spherical coordinates :math:`(\\phi, \\theta)`
    to a unit vector :math:`(x, y, z)`

    :param phi:
    :type phi: :class:`np.ndarray` or float
    :param theta:
    :type theta: :class:`np.ndarray` or float
    :return: :math:`(x, y, z)`
    :rtype: tuple
    """
    x = np.sin(phi) * np.cos(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(phi)
    return x, y, z