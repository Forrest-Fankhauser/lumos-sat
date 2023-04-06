import lumos.geometry
import lumos.brdf.library as brdf_library
import lumos.calculators.ground_frame as calculator
import numpy as np
import lumos.plot.ground_frame as lumos_plot

nadir1 = lumos.geometry.Surface(
    area = 10.0,
    normal = np.array([0, 0, -1]), 
    brdf = brdf_library.LAMBERTIAN(1)
    )

nadir2 = lumos.geometry.Surface(
    area = 10.0,
    normal = np.array([0, 0, -1]),
    brdf = brdf_library.ABG(0.1, 0.1, 2)
)

satellite_surfaces = [nadir1, nadir2]
angle_past_terminator = np.deg2rad(10)
satellite_height = 550 * 1000

altitudes = np.linspace(0, 90, 45)
azimuths = np.linspace(0, 360, 90)
altitudes, azimuths = np.meshgrid(altitudes, azimuths)

sun_altitude = -20
sun_azimuth = 270

intensities = calculator.calculate_intensity(
    satellite_surfaces,
    satellite_height, 
    sat_altitudes = altitudes,
    sat_azimuths = azimuths, 
    sun_altitude = sun_altitude,
    sun_azimuth = sun_azimuth,
    include_sun = True,
    include_earthshine = False,
    earth_panel_density = 101, 
    earth_brdf = None)

lumos_plot.plot_AB_Mag_contour(
    altitudes,
    azimuths,
    intensities,
    sun_altitude,
    sun_azimuth,
    levels = (0, 8),
    filter = True,
    title = "Simple Satellite in Ground Frame",
    save_file = "ground_frame.png"
    )