""" Library of Satellite Constellations """

import numpy as np
import sgp4.api
import astropy.time
import astropy.units
import astropy.coordinates
import matplotlib.pyplot as plt

class Constellation:
    """
    Base class for all constellations
    """
    
    def get_teme_position(self, time):
        """
        Gets positions of all satellites in TEME coordinate frame

        :param time: Time at which to get positions
        :type time: :class:`astropy.time.Time`
        :return: (x, y, z) position measured in meters
        :rtype: tuple[:class:`np.ndarray`]
        """

        jd1 = np.array([time.jd1])
        jd2 = np.array([time.jd2])

        _, teme_pos, _ = self.constellation.sgp4(jd1, jd2)
        teme_pos = teme_pos[:,0,:]
        x, y, z = teme_pos[:,0], teme_pos[:,1], teme_pos[:,2]

        return x * 1000, y * 1000, z * 1000
    
    def get_hcs_position(self, time, earth_location):
        """
        Gets position of all satellites in constellation in HCS frame
        
        :param time: Time at which to get positions
        :type time: :class:`astropy.time.Time`
        :param earth_location: Location at which to get positions in HCS frame
        :type earth_location: :class:`astropy.coordinates.EarthLocation`
        :return: altitudes (degrees), azimuths (degrees), heights (meters)
        :rtype: tuple[:class:`np.ndarray`]
        """

        x, y, z = self.get_teme_position(time)
        teme = astropy.coordinates.TEME(
            x = x * astropy.units.meter,
            y = y * astropy.units.meter,
            z = z * astropy.units.meter,
            representation_type = 'cartesian',
            obstime = time)
        
        itrs_geo = teme.transform_to(astropy.coordinates.ITRS(obstime = time))
        topo_itrs_repr = itrs_geo.cartesian.without_differentials() \
                        - earth_location.get_itrs(time).cartesian
        itrs_topo = astropy.coordinates.ITRS(topo_itrs_repr, obstime = time, location = earth_location)
        aa = itrs_topo.transform_to(
            astropy.coordinates.AltAz(obstime = time, 
                                      location = earth_location,
                                      pressure = 0)
                                      )
        
        heights = itrs_geo.earth_location.geodetic.height.value

        return aa.alt.degree, aa.az.degree, heights
    
    def plot_teme(self, ax, time):
        """
        Plots satellite constellation in TEME coordinate frame

        :param ax: Axis for plotting onto (must be 3D projection)
        :type ax: :class:`matplotlib.pyplot.Axes`
        :param time: Time at which to plot constellation
        :type time: :class:`astropy.time.Time`
        """

        x, y, z = self.get_teme_position(time)
        x, y, z = x / 1000, y / 1000, z / 1000

        ax.scatter(x, y, z, s = 1, alpha = 0.5)
        ax.set_box_aspect((np.ptp(x), np.ptp(y), np.ptp(z)))
        ax.xaxis.set_pane_color((1, 1, 1, 0))
        ax.yaxis.set_pane_color((1, 1, 1, 0))
        ax.zaxis.set_pane_color((1, 1, 1, 0))
        ax.set_xlabel("x (km)")
        ax.set_ylabel("y (km)")
        ax.set_zlabel("z (km)")
    
    def plot_hcs(self, ax, time, location):
        """
        Plots satellite constellation in HCS coordinate frame

        :param ax: Axis for plotting onto (must be polar projection)
        :type ax: :class:`matplotlib.pyplot.Axes`
        :param time: Time at which to plot constellation
        :type time: :class:`astropy.time.Time`
        :param location: Location on Earth at which to plot constellation
        :type location: :class:`astropy.coordinates.EarthLocation`
        """

        altitudes, azimuths, _ = self.get_hcs_position(time, location)
        ax.plot( np.deg2rad(azimuths), 90 - altitudes, '.')
        ax.set_rmax(90)
        ax.set_yticklabels([])
        ax.set_theta_zero_location('N')
        ax.set_rticks([10, 20, 30, 40, 50, 60, 70, 80, 90])
        ax.set_xticks(np.deg2rad([0, 90, 180, 270]))
        ax.set_xticklabels(['N', 'E', 'S', 'W'])

class ConstellationFromTLE(Constellation):
    """
    This class builds a constellation from a TLE file

    :param tle_file_path: Path to a TLE text file
    :type tle_file_path: str
    """

    def __init__(self, tle_file_path):
        """
        Constructor Method
        """
        
        with open(tle_file_path, 'r') as file:
            lines = file.read().splitlines()

        TLES = [ (l0, l1, l2) for l0, l1, l2 in zip(lines[::3], lines[1::3], lines[2::3]) ]
    
        satellite_ids = []
        constellation = []

        for tle in TLES:
            satellite_ids.append( tle[0].strip() )
            constellation.append(
                sgp4.api.Satrec.twoline2rv(tle[1], tle[2])
                )
        
        self.constellation = sgp4.api.SatrecArray(constellation)
        self.satellite_ids = tuple(satellite_ids)
        self.constellation_size = len( self.satellite_ids )