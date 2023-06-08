"""
Geometry objects which are used throughout Lumos
"""

import numpy as np
import lumos.constants
import lumos.calculator
class Surface:
    """
    Container to hold area, normal vector, and BRDF of a surface 

    :param area: Area of surface :math:`m^2`
    :type area: float
    :param normal: Normal vector of surface. Measured in brightness frame. A function of
        the angle past terminator may be passed, which must return the surface normal as 
        a :class:`np.ndarray` for given angles past the terminator.
    :type normal: :class:`np.ndarray` or function
    :param brdf: Bidirectional Reflectance Distribution Function (BRDF) of surface.
    :type brdf: function
    """
    def __init__(self, area, normal, brdf):
        """
        Constructor Method
        """
        self.area = area
        self.normal = normal
        self.brdf = brdf
    
    def __str__(self):
        """
        Printing Method
        """
        return f'| Surface \n' \
               f'|-- Area: {self.area:.2f} m^2 \n' \
               f'|-- Normal Vector: ' \
               f'<{self.normal[0]:.2f}, {self.normal[1]:.2f}, {self.normal[2]:.2f}> \n'

class EarthMesh:
    """
    A mesh of points on Earth's surface

    :param angles_off_plane: 1D array of the angles-off-plane of the mesh.
    :type angles_off_plane: :class:`np.ndarray`
    :param angles_on_plane: 1D array of the angles-on-plane of the mesh.
    :type angles_on_plane: :class:`np.ndarray`
    """
    def __init__(self, angles_off_plane, angles_on_plane):
        """
        Constructor method
        """
        self.d_phi = np.abs(angles_off_plane[1] - angles_off_plane[0])
        self.d_theta = np.abs(angles_on_plane[1] - angles_on_plane[0])
        
        self.angles_off_plane, self.angles_on_plane = np.meshgrid(angles_off_plane,
                                                                  angles_on_plane)
        
        self.dists_off_plane = lumos.constants.EARTH_RADIUS * self.angles_off_plane
        self.dists_on_plane = lumos.constants.EARTH_RADIUS * self.angles_on_plane
        
        # Set up panel positions
        self.z = ( lumos.constants.EARTH_RADIUS 
                  / np.sqrt(1 + np.tan(self.angles_on_plane)**2 
                            + np.tan(self.angles_off_plane)**2 ) )
                  
        self.x = np.tan(self.angles_off_plane) * self.z
        self.y = np.tan(self.angles_on_plane) * self.z
        
        # Get normal vectors of panels
        self.nx = self.x / lumos.constants.EARTH_RADIUS
        self.ny = self.y / lumos.constants.EARTH_RADIUS
        self.nz = self.z / lumos.constants.EARTH_RADIUS

        self.shape = self.z.shape
        
class GroundObservers(EarthMesh):
    """
    A mesh of observers visible to the satellite and on the night side of earth

    :param sat_height: Geodetic height of satellite (meters)
    :type sat_height: float
    :param angle_past_terminator: Angle of satellite past the terminator (radians)
    :type angle_past_terminator: float
    :param density: Mesh will have size density x density
    :type density: int
    """

    def __init__(self, sat_height, angle_past_terminator, density):
        """
        Constructor Method
        """
        self.max_angle = np.arccos(lumos.constants.EARTH_RADIUS 
                                   / (lumos.constants.EARTH_RADIUS + sat_height))
        self.min_angle = angle_past_terminator
        
        angles_off_plane = np.linspace(-self.max_angle, self.max_angle, density)
        angles_on_plane = np.linspace(-self.max_angle, self.min_angle, density)
        
        super().__init__(angles_off_plane, angles_on_plane)
        self.sat_height = sat_height
        self.angle_past_terminator = angle_past_terminator
    
    def calculate_intensity(self, surfaces,
                            include_sun = True, include_earthshine = False, 
                            earth_panel_density = 151, earth_brdf = None):
        """
        Calculates intensity for observers on ground
        
        :param surfaces: List of satellite surfaces
        :type surfaces: List[:py:func:`lumos.geometry.Surface`]
        :param include_sun: Whether to include contribution of brightness due to direct sunlight
        :type include_sun: bool, optional
        :param include_earthshine: Whether to include contribution of brightness due to earthshine
        :type include_earthshine: bool, optional
        :param earth_panel_density: Earthshine discretization has earth_panel_density squared pixels
        :type earth_panel_density: int, optional
        :param earth_brdf: The BRDF of Earth's surface
        :type earth_brdf: callable
        """

        # Sets up an array to hold the intensity seen by each ground observer
        self.intensities = np.zeros(self.shape)

        for i in range(self.shape[0]):
            for j in range(self.shape[1]):
                
                self.intensities[i, j] = \
                    lumos.calculator.get_intensity_satellite_frame(
                        surfaces,
                        self.sat_height,
                        self.angle_past_terminator,
                        (self.x[i, j], self.y[i, j], self.z[i, j]),
                        include_sun,
                        include_earthshine,
                        earth_panel_density,
                        earth_brdf)