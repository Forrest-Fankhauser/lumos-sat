"""
Intensity calculations in an observer's reference frame.
"""

import numpy as np
import lumos.calculators.brightness_frame as brightness_frame
import lumos.constants
import lumos.functions
import lumos.conversions
import lumos.geometry

def get_brightness_coords(
    sat_alt, 
    sat_az, 
    sat_height, 
    sun_alt, 
    sun_az):
    """
    Converts from HCS observer coordinates to brightness coordinates.

    :param sat_alt: Altitude of satellite (degrees)
    :type sat_alt: float or :class:`np.ndarray`
    :param sat_az: Altitude of satellite (degrees)
    :type sat_az: float or :class:`np.ndarray`
    :param sat_height: Geodetic height of satellite (meters)
    :type sat_height: float or :class:`np.ndarray`
    :param sun_alt: Altitude of sun (degrees)
    :type sun_alt: float
    :param sun_az: Azimuth angle of sun (degrees)
    :type sun_az: float
    :returns:
        - (obs_x, obs_y, obs_z) - Vector from satellite to observer in brightness frame (meters)
        - angle_past_terminator - Angle of satellite past terminator (radians)
    """

    sat_x, sat_y, sat_z = lumos.conversions.altaz_to_unit(sat_alt, sat_az)
    sun_x, sun_y, sun_z = lumos.conversions.altaz_to_unit(sun_alt, sun_az)

    phi = np.arccos(sat_z) - np.arcsin(lumos.constants.EARTH_RADIUS / (lumos.constants.EARTH_RADIUS + sat_height) * np.sqrt(sat_x**2 + sat_y**2) )
    theta = np.arctan2(sat_y, sat_x)

    Zx = np.sin(phi) * np.cos(theta)
    Zy = np.sin(phi) * np.sin(theta)
    Zz = np.cos(phi)

    dot = Zx * sun_x + Zy * sun_y + Zz * sun_z
    beta = 1 / np.sqrt(1 - dot**2)
    alpha = - beta * dot

    Yx = alpha * Zx + beta * sun_x
    Yy = alpha * Zy + beta * sun_y
    Yz = alpha * Zz + beta * sun_z

    flip = np.sign(Yx * sun_x + Yy * sun_y + Yz * sun_z)
    Yx, Yy, Yz = flip * Yx, flip * Yy, flip * Yz

    Xx = Yy * Zz - Zy * Yz
    Xy = Zx * Yz - Yx * Zz
    Xz = Yx * Zy - Zx * Yy

    T11, T12, T13, \
    T21, T22, T23, \
    T31, T32, T33  = lumos.functions.inv_3(Xx, Yx, Zx, 
                                           Xy, Yy, Zy, 
                                           Xz, Yz, Zz)
    
    dist_to_sat = np.sqrt(lumos.constants.EARTH_RADIUS**2 
    + (lumos.constants.EARTH_RADIUS + sat_height)**2 
    - 2 * lumos.constants.EARTH_RADIUS * (lumos.constants.EARTH_RADIUS + sat_height) * np.cos(phi))

    obs_x = - dist_to_sat * (T11 * sat_x + T12 * sat_y + T13 * sat_z)
    obs_y = - dist_to_sat * (T21 * sat_x + T22 * sat_y + T23 * sat_z)
    obs_z = (lumos.constants.EARTH_RADIUS + sat_height) - dist_to_sat * (T31 * sat_x + T32 * sat_y + T33 * sat_z)

    angle_past_terminator = -np.arcsin(T31 * sun_x + T32 * sun_y + T33 * sun_z)

    return obs_x, obs_y, obs_z, angle_past_terminator

def calculate_intensity(
    sat_surfaces, 
    sat_heights, 
    sat_altitudes, 
    sat_azimuths,
    sun_altitude, 
    sun_azimuth,
    include_sun = True,
    include_earthshine = True,
    earth_panel_density = 150,
    earth_brdf = None):

    """
    Calculates the flux of a satellite seen by an observer.

    :param sat_surfaces: List of surfaces on satellite
    :type sat_surfaces: list[lumos.geometry.Surface]
    :param sat_heights: Heights of satellite above geodetic nadir (meters)
    :type sat_heights: float or :class:`np.ndarray`
    :param sat_altitudes: Satellite altitudes in HCS frame (degrees)
    :type sat_altitudes: float or :class:`np.ndarray`
    :param sat_azimuths: Satellite azimuths in HCS frame (degrees)
    :type sat_azimuths: float or :class:`np.ndarray`
    :param sun_alitutde: Altitude of sun in HCS frame (degrees)
    :type sun_altitude: float
    :param sun_azimuth: Azimuth of sun in HCS frame (degrees)
    :type sun_azimuth: float
    :param include_sun: Whether to include flux scattered by the satellite from the sun
    :type include_sun: bool, optional
    :param include_earthshine: Whether to include flux scattered by the satellite from earthshine
    :type include_earthshine: bool, optional
    :param earth_panel_density: There will be earth_panel_density x earth_panel_density panels in the earthshine mesh
    :type earth_panel_density: int, optional
    :param earth_brdf: A function representing the BRDF of Earth's surface
    :type earth_brdf: function
    :return: Flux of light incident on the observer (W / m^2)
    :rtype: float or :class:`np.ndarray`
    """
    
    if sun_altitude > 0:
        raise ValueError("Observatory is in daylight!")
    
    sat_alts, sat_azs = np.asarray(sat_altitudes), np.asarray(sat_azimuths)
    output_shape = sat_alts.shape

    if not isinstance(sat_heights, np.ndarray):
        sat_heights = sat_heights * np.ones(output_shape)

    sat_heights, sat_alts, sat_azs = sat_heights.flatten(), sat_alts.flatten(), sat_azs.flatten()
    intensity = np.zeros_like(sat_alts)

    obs_x, obs_y, obs_z, angle_past_terminator = get_brightness_coords(sat_alts, sat_azs, sat_heights, 
        sun_altitude, sun_azimuth)

    for (i, (x, y, z, alpha, h)) in enumerate(zip(obs_x, obs_y, obs_z, angle_past_terminator, sat_heights)):
        intensity[i] = brightness_frame.calculate_intensity(
            sat_surfaces, h, alpha, (x, y, z),
            include_sun, include_earthshine, earth_panel_density, earth_brdf)

    intensity = intensity.reshape(output_shape)
    return intensity