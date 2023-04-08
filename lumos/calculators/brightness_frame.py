"""
Intensity calculations in the brightness reference frame.
"""
import numpy as np
import lumos.constants
import lumos.geometry

def get_earthshine_panels(sat_z : float, angle_past_terminator : float, density : int) \
      -> tuple[np.ndarray, ...]:
    """
    Creates mesh of panels on Earth's surface.

    Parameters:
        sat_z (float): The height of the satellite above geodetic nadir (meters)
        angle_past_terminator (float): The angle of the satellite past terminator (radians)
        density (int): The density of the panels. Grid will have size density x 
    
    Returns:
        x (np.ndarray): x positions of panels (meters)
        y (np.ndarray): y positions of panels (meters)
        z (np.ndarray): z positions of panels (meters)
        nx (np.ndarray): x component of panel normals
        ny (np.ndarray): y component of panel normals
        nz (np.ndarray): z component of panel normals
        areas (np.ndarray): areas of panels (meters ^ 2)
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

def calculate_intensity(sat_surfaces : list[lumos.geometry.Surface], 
                        sat_altitude : float, 
                        angle_past_terminator : float,
                        observer_position : float,
                        include_sun : bool = True, 
                        include_earthshine : bool = True,
                        earth_panel_density : int = 150, 
                        earth_brdf : callable = None) -> float:
    '''
        Parameters:
            sat_surfaces (list[lumos.geometry.Surface]) : A list of lumos.geometry.general.Surface objects
            sat_altitude (float) : Height of satellite above geodetic nadir (meters)
            angle_past_terminator (float) : Angle of satellite past terminator (radians)
            include_sun (bool) : Whether to include brightness directly reflected from sun
            include_earthshine (bool) : Whether to include brightness caused by earthshine
            earth_panel_density (int) : There will be earth_panel_density^2 panels in the earthshine mesh
            earth_brdf (callable) : A function representing the BRDF of earth
        
        Returns:
            intensity (float) : The intensity of light incident on the observer (W / m^2)
    '''
    
    horizon_angle = np.arccos(lumos.constants.EARTH_RADIUS / (lumos.constants.EARTH_RADIUS + sat_altitude))
    if angle_past_terminator > horizon_angle:
        # Inside earth's shadow
        return 0
    
    observer_x, observer_y, observer_z = observer_position[0], observer_position[1], observer_position[2]

    if np.arccos(observer_z / lumos.constants.EARTH_RADIUS) > horizon_angle + np.deg2rad(1):
        # Not visible to observer
        # Fudging this by 1 degree makes plots much nicer on edges
        return 0

    vector_2_sun = (0, np.cos(angle_past_terminator), - np.sin(angle_past_terminator))

    sat_z = sat_altitude + lumos.constants.EARTH_RADIUS

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

        if include_sun:
            surface_brdf = surface.brdf(vector_2_sun,
                                                surface_normal,
                                                (sat_obs_x, sat_obs_y, sat_obs_z))

            sun_contribution = sun_contribution + surface_brdf * surface_normalization

        if include_earthshine:
            surface_normalizations = -( surface_normal[0] * panel_sat_x 
                             + surface_normal[1] * panel_sat_y
                             + surface_normal[2] * panel_sat_z)

            surface_brdf = surface.brdf( (-panel_sat_x, -panel_sat_y, -panel_sat_z),
                                                   surface_normal,
                                                   (sat_obs_x, sat_obs_y, sat_obs_z) )
            
            earth_contribution = earth_contribution \
                  + np.sum(panel_brdfs * surface_brdf 
                         * panel_normalizations * surface_normalizations
                         * panel_areas / dist_panels_2_sat**2 )
        
        intensity = intensity + lumos.constants.SUN_INTENSITY * surface.area * (sun_contribution + earth_contribution) / dist_sat_2_obs ** 2

    return intensity