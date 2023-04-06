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
    sat_alt : float | np.ndarray, 
    sat_az : float | np.ndarray, 
    sat_height : float | np.ndarray, 
    sun_alt : float, 
    sun_az : float) \
    -> tuple[float | np.ndarray, ...]:
    '''
    Converts from HCS observer coordinates to brightness coordinates.
        Parameters:
            sat_alt (float or np.ndarray) : Altitude of satellite (degrees)
            sat_az (float or np.ndarray) : Azimuth of satellite (degrees)
            sat_height (float or np.ndarray) : Geodetic height of satellite (meters)
            sun_alt (float) : Altitude of sun (degrees)
            sun_az (float) : Azimuth of sun (degrees)
        Returns:
            obs_x (np.ndarray) : x-component of vector from satellite to observer in brightness frame (meters)
            obs_y (np.ndarray) : y-component of vector from satellite to observer in brightness frame (meters)
            obs_z (np.ndarray) : z-component of vector from satellite to observer in brightness frame (meters)
            angle_past_terminator (np.ndarray) : Angle of satellite past the terminator (radians)
    '''

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
    sat_surfaces : list[lumos.geometry.Surface], 
    sat_heights : np.ndarray | float, 
    sat_altitudes : np.ndarray | float, 
    sat_azimuths : np.ndarray | float,
    sun_altitude : float, 
    sun_azimuth : float, 
    **kwargs) -> np.ndarray | float:
    '''
    Calculates the intensity of a satellite as seen in a ground observer's reference frame.
        Parameters:
            sat_surfaces (list[Surface]) : A list of lumos.geometry.general.Surface objects
            sat_heights (np.ndarray or float) : Heights of satellites above geodetic nadir (meters)
            sat_altitudes (np.ndarray or float) : Altitude of satellites (degrees)
            sat_azimuths (np.ndarray or float) : Azimuths of satellites (degrees)
            sun_altitude (float) : Altitude of sun
            sun_azimuth (float) : Azimuth of sun
        Keyword Arguments:
            include_sun (bool) : Whether to include brightness directly reflected from sun
            include_earthshine (bool) : Whether to include brightness caused by earthshine
            earth_panel_density (int) : There will be earth_panel_density^2 panels in the earthshine mesh
            earth_brdf (callable) : A function representing the BRDF of earth
        Returns:
            intensity (np.ndarray) : The intensity of light incident on the observer (W / m^2)
    '''

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
        intensity[i] = brightness_frame.calculate_intensity(sat_surfaces, h, alpha, (x, y, z), **kwargs)
    
    intensity = intensity.reshape(output_shape)
    return intensity