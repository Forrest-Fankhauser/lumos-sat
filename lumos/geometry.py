"""
Geometry objects which are used throughout Lumos
"""
import numpy as np
import lumos.constants
class Surface:
    """ Object to hold area, normal vector, and BRDF of a surface """
    def __init__(self, area : float, normal : np.ndarray, brdf : callable):
        """
        Parameters:
            area (float) : Area of surface (meters ^ 2)
            normal (np.ndarray) : Normal vector of surface. Measured in brightness frame.
            brdf (callable) : BRDF of surface
        """
        self.area = area
        self.normal = normal
        self.brdf = brdf
    
    def __str__(self):
        return f'| Surface \n' \
               f'|-- Area: {self.area:.2f} m^2 \n' \
               f'|-- Normal Vector: ' \
               f'<{self.normal[0]:.2f}, {self.normal[1]:.2f}, {self.normal[2]:.2f}> \n'

class EarthMesh:
    """
    A mesh of points on Earth's surface
    """
    def __init__(self, angles_off_plane : np.ndarray, angles_on_plane : np.ndarray):
        """
        Parameters:
            angles_off_plane (np.ndarray) : 1D array of the angles-off-plane of the mesh
            angles_on_plane (np.ndarray) : 1D array of the angles-on-plane of the mesh
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
    """

    def __init__(self, sat_height : float, angle_past_terminator : float, density : int):
        """
        Parameters:
            sat_height (float) : Geodetic height of satellite (meters)
            angle_past_terminator (float) : Angle of satellite past the terminator (radians)
            density (int) : Mesh will have size density x density
        """
        self.max_angle = np.arccos(lumos.constants.EARTH_RADIUS 
                                   / (lumos.constants.EARTH_RADIUS + sat_height))
        self.min_angle = angle_past_terminator
        
        angles_off_plane = np.linspace(-self.max_angle, self.max_angle, density)
        angles_on_plane = np.linspace(-self.max_angle, self.min_angle, density)
        
        super().__init__(angles_off_plane, angles_on_plane)
    
    def __iter__(self):
        for i in range(self.shape[0]):
            for j in range(self.shape[1]):
                yield i, j, (self.x[i, j], self.y[i, j], self.z[i, j])