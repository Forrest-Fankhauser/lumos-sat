"""
Intensity calculations in the brightness reference frame.
"""

import numpy as np
import lumos.constants
import lumos.geometry

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

def calculate_intensity(
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
    :param earth_panel_density: There will be earth_panel_density x earth_panel_density 
    panels in the earthshine mesh
    :type earth_panel_density: int
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
            surface_normalizations = np.clip( surface_normal[0] * panel_sat_x 
                             + surface_normal[1] * panel_sat_y
                             + surface_normal[2] * panel_sat_z, 0, None)**2

            surface_brdf = surface.brdf( (-panel_sat_x, -panel_sat_y, -panel_sat_z),
                                          surface_normal,
                                         (sat_obs_x, sat_obs_y, sat_obs_z) )
            
            earth_contribution = earth_contribution \
                  + np.sum(panel_brdfs * surface_brdf 
                         * panel_normalizations * surface_normalizations
                         * panel_areas / dist_panels_2_sat**2 )
        
        intensity = intensity + lumos.constants.SUN_INTENSITY * surface.area * (sun_contribution + earth_contribution) / dist_sat_2_obs ** 2

    return intensity