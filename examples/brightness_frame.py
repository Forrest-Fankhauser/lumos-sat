import lumos.geometry
import lumos.brdf.library as brdf_library
import lumos.calculators.brightness_frame as calculator
import numpy as np
import lumos.plot.brightness_frame as lumos_plot

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

# Define observer mesh
observers = lumos.geometry.GroundObservers(satellite_height, angle_past_terminator, 51)

# Calculations
intensities = np.zeros(observers.shape)

for i, j, position in observers:
    intensities[i, j] = calculator.calculate_intensity(
        satellite_surfaces,
        satellite_height,
        angle_past_terminator,
        position,
        include_sun = True,
        include_earthshine = False,
        earth_panel_density = 101,
        earth_brdf = None)

lumos_plot.plot_AB_Mag(
    observers,
    intensities,
    title = "Simple Satellite in Brightness Frame",
    save_file = "brightness_frame.png",
    levels = (0, 8))