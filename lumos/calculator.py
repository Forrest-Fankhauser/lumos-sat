"""
Main brightness calculator for Lumos
"""

import numpy as np
import lumos.constants
import lumos.geometry
import astropy.coordinates

def get_earthshine_panels(sat_z, angle_past_terminator, density):
    """
    Creates a mesh of pixels on Earth's surface which are visible to the satellite and illuminated
    by the sun.

    :param sat_z: The height of the satellite above the center of Earth (meters)
    :type sat_z: float
    :param angle_past_terminator: The angle of the satellite past the terminator (radians)
    :type angle_past_terminator: float
    :param density: The density of the pixels. Grid will have size density x density.
    :type density: int
    :returns: 
        - (x, y, z) - Positions of pixels (meters)
        - (nx, ny, nz) - Normal vectors of pixels
        - areas - Areas of pixels (:math:`m^2`)
    """

    R = lumos.constants.EARTH_RADIUS

    max_angle = np.arccos(R / sat_z)

    angles_off_plane = np.linspace(-max_angle, max_angle, density)
    angles_on_plane = np.linspace(angle_past_terminator, max_angle, density)

    d_phi = abs( angles_off_plane[1] - angles_off_plane[0] )
    d_theta = abs( angles_on_plane[1] - angles_on_plane[0] )

    angles_on_plane, angles_off_plane = np.meshgrid(angles_on_plane, angles_off_plane)
    angles_on_plane, angles_off_plane = angles_on_plane.flatten(), angles_off_plane.flatten()

    # Set up panel positions
    nz = 1 / np.sqrt( 1 + np.tan(angles_on_plane)**2 + np.tan(angles_off_plane)**2 )
    nx = np.tan(angles_off_plane) * nz
    ny = np.tan(angles_on_plane) * nz

    # Clip the panels which aren't visible to the satellite
    visible_to_sat = np.arccos(nz) < max_angle
    angles_off_plane, angles_on_plane = angles_off_plane[visible_to_sat], angles_on_plane[visible_to_sat]
    nx, ny, nz = nx[visible_to_sat], ny[visible_to_sat], nz[visible_to_sat]

    # Calculate Jacobian determinant to get panel areas
    x, y, z = nx * R, ny * R, nz * R

    phi = angles_off_plane
    theta = angles_on_plane
        
    dx_dr = nx / nz * z / R
    dx_dphi = z**3 / (R**2 * np.cos(phi)**2 * np.cos(theta)**2)
    dx_dtheta = - (ny / nz * nx / nz* z**3 ) / (R**2 * np.cos(theta)**2 )
    
    dy_dr = np.tan(theta) * z / R
    dy_dphi = - ( ny / nz * nx / nz * z**3 ) / (R**2 * np.cos(phi)**2 )
    dy_dtheta = dx_dphi
    
    dz_dr = z / R
    dz_dphi = - (nx / nz * z**3 ) / (R**2 * np.cos(phi)**2)
    dz_dtheta = - (ny / nz * z**3 ) / (R**2 * np.cos(theta)**2)

    determinant = (
        dx_dr * (dy_dphi * dz_dtheta - dy_dtheta * dz_dphi) -
        dy_dr * (dx_dphi * dz_dtheta - dx_dtheta * dz_dphi) +
        dz_dr * (dx_dphi * dy_dtheta - dx_dtheta * dy_dphi) )
    
    areas = determinant * d_phi * d_theta
    
    return x, y, z, nx, ny, nz, areas

def get_intensity_satellite_frame(
    sat_surfaces, 
    sat_height, 
    angle_past_terminator,
    observer_position,
    include_sun = True, 
    include_earthshine = True,
    earth_panel_density = 150, 
    earth_brdf = None):
    '''
    Calculates flux scattered by a satellite and seen by an observer.

    :param sat_surfaces: List of surfaces on satellite
    :type sat_surfaces: list[lumos.geometry.Surface]
    :param sat_height: Height of satellite above geodetic nadir (meters)
    :type sat_height: float
    :param angle_past_terminator: Angle of satellite past terminator (radians)
    :type angle_past_terminator: float
    :param observer_position: Position of an observer, measured in brightness frame (meters)
    :type observer_position: :class:`np.ndarray`
    :param include_sun: Whether to include flux scattered by the satellite from the sun
    :type include_sun: bool, optional
    :param include_earthshine: Whether to include flux scattered by the satellite from earthshine
    :type include_earthshine: bool, optional
    :param earth_panel_density: There will be earth_panel_density x earth_panel_density panels in the earthshine mesh
    :type earth_panel_density: int, optional
    :param earth_brdf: A function representing the BRDF of Earth's surface
    :type earth_brdf: function
    :return: Flux of light incident on the observer (W / m^2)
    :rtype: float
    '''
    
    horizon_angle = np.arccos(lumos.constants.EARTH_RADIUS / (lumos.constants.EARTH_RADIUS + sat_height))
    if angle_past_terminator > horizon_angle:
        # Inside earth's shadow
        return 0
    
    observer_x, observer_y, observer_z = observer_position[0], observer_position[1], observer_position[2]

    if np.arccos(observer_z / lumos.constants.EARTH_RADIUS) > horizon_angle + np.deg2rad(1):
        # Not visible to observer
        # Fudging this by 1 degree makes plots much nicer on edges
        return 0

    vector_2_sun = (0, np.cos(angle_past_terminator), - np.sin(angle_past_terminator))

    sat_z = sat_height + lumos.constants.EARTH_RADIUS

    # Distances from observers to satellite
    dist_sat_2_obs = np.sqrt( observer_x**2 + observer_y**2 
                            + (observer_z - sat_z)**2 )
    
    # Unit vectors from observers to satellite
    sat_obs_x = observer_x / dist_sat_2_obs
    sat_obs_y = observer_y / dist_sat_2_obs
    sat_obs_z = (observer_z - sat_z) / dist_sat_2_obs

    if include_earthshine:
        panel_x, panel_y, panel_z, panel_nx, panel_ny, panel_nz, panel_areas = \
            get_earthshine_panels(sat_z, angle_past_terminator, earth_panel_density)
        
        # Distances from earthshine panels to satellite
        dist_panels_2_sat = np.sqrt(panel_x**2 + panel_y**2 + (sat_z - panel_z)**2 )

        panel_sat_x = - panel_x / dist_panels_2_sat
        panel_sat_y = - panel_y / dist_panels_2_sat
        panel_sat_z = (sat_z - panel_z) / dist_panels_2_sat

        # Panel Normalization
        panel_normalizations = (panel_ny * vector_2_sun[1] + panel_nz * vector_2_sun[2])

        panel_brdfs = earth_brdf(vector_2_sun,
                                (panel_nx, panel_ny, panel_nz),
                                (panel_sat_x, panel_sat_y, panel_sat_z))

    intensity = 0

    for surface in sat_surfaces:
        
        if callable(surface.normal):
            surface_normal = surface.normal(angle_past_terminator)
        else:
            surface_normal = surface.normal
        
        sun_contribution = 0
        earth_contribution = 0
        
        surface_normalization = surface_normal[1] * vector_2_sun[1] + surface_normal[2] * vector_2_sun[2]
        surface_normalization = np.clip(surface_normalization, 0, None)

        observer_normalization = surface_normal[0] * sat_obs_x + surface_normal[1] * sat_obs_y + surface_normal[2] * sat_obs_z
        observer_normalization = np.clip(observer_normalization, 0, None)

        if include_sun:
            surface_brdf = surface.brdf(vector_2_sun,
                                        surface_normal,
                                        (sat_obs_x, sat_obs_y, sat_obs_z))

            sun_contribution = sun_contribution + surface_brdf * surface_normalization * observer_normalization

        if include_earthshine:

            surface_normalizations = np.clip(
                - surface_normal[0] * panel_sat_x \
                - surface_normal[1] * panel_sat_y \
                - surface_normal[2] * panel_sat_z,
                0,
                None
            )

            panel_observing_normalization = np.clip(
                panel_nx * panel_sat_x
                + panel_ny * panel_sat_y
                + panel_nz * panel_sat_z,
                0,
                None
            )

            surface_brdf = surface.brdf( (-panel_sat_x, -panel_sat_y, -panel_sat_z),
                                          surface_normal,
                                         (sat_obs_x, sat_obs_y, sat_obs_z) )
            
            earth_contribution = earth_contribution \
                  + np.sum(panel_brdfs * surface_brdf 
                         * panel_normalizations * observer_normalization 
                         * surface_normalizations * panel_observing_normalization
                         * panel_areas / dist_panels_2_sat**2 )
        
        intensity = intensity + lumos.constants.SUN_INTENSITY * surface.area * (sun_contribution + earth_contribution) / dist_sat_2_obs ** 2

    return intensity

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

def get_intensity_observer_frame(
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
        intensity[i] = get_intensity_observer_frame(
            sat_surfaces, h, alpha, (x, y, z),
            include_sun, include_earthshine, earth_panel_density, earth_brdf)

    intensity = intensity.reshape(output_shape)
    return intensity

def get_sun_alt_az(time, observer_location):
    
    """
    Convenience function for finding the altitude and azimuth of the sun

    :param time: Time of observation
    :type time: :class:`astropy.time.Time`
    :param observer_location: Location of observation
    :type observer_location: :class:`astropy.coordinates.EarthLocation`
    """

    aa_frame = astropy.coordinates.AltAz(obstime = time, location = observer_location)
    sun_altaz = astropy.coordinates.get_sun(time).transform_to(aa_frame)
    return sun_altaz.alt.degree, sun_altaz.az.degree

