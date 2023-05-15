import numpy as np
import matplotlib.pyplot as plt
import satellite
import lumos.calculator
from lumos.geometry import GroundObservers
import lumos.plot
import lumos.conversions
import lumos.brdf.library

sat_height = 550 * 1000 # Geodetic height of satellite (meters)

def plot_brightness():
    fig, axs = plt.subplots(1, 4)
    angles_past_terminator = (-10, 0, 10, 20)

    for ax, angle_past_terminator in zip(axs, angles_past_terminator):
        
        # Creates a mesh of observers on Earth's surface
        observers = GroundObservers(
            sat_height, 
            np.deg2rad(angle_past_terminator), 
            density = 30)
        
        # Sets up an array to hold the intensity seen by each ground observer
        intensities = np.zeros(observers.shape)

        # Loops through each observer and calculates intensity seen
        # by that observer
        for i, j in observers:
            intensities[i, j] = \
                lumos.calculator.get_intensity_satellite_frame(
                    satellite.SURFACES,
                    sat_height,
                    angle_past_terminator,
                    (observers.x[i, j], observers.y[i, j], observers.z[i, j]),
                    include_sun = True,
                    include_earthshine = True,
                    earth_panel_density = 150,
                    earth_brdf = lumos.brdf.library.LAMBERTIAN(1)
                    )
        
        # Convert intensity to AB Magnitude
        ab_magnitudes = lumos.conversions.intensity_to_ab_mag(intensities)

        # Plots intensity
        lumos.plot.contour_satellite_frame(
            ax,
            observers,
            ab_magnitudes,
            levels = (2, 10),
            cmap = 'plasma_r'
            )
        
        ax.set_title(f"{angle_past_terminator}")

    plt.show()