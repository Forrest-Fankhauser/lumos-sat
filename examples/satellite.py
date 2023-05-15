# Satellite Brightness Model

from lumos.geometry import Surface
import lumos.brdf.library
import numpy as np

x_hat = np.array([1, 0, 0])
y_hat = np.array([0, 1, 0])
z_hat = np.array([0, 0, 1])

chassis = Surface(
    area = 1.0, # meters^2
    normal = -z_hat, # To geodetic nadir
    brdf = lumos.brdf.library.PHONG(1, 1, 2) # BRDF of chassis
    )

solar_array = Surface(
    area = 1.0, # meters^2
    normal = y_hat, # Perpendicular to the chassis
    brdf = lumos.brdf.library.PHONG(1, 1, 1) # BRDF of solar array
    )

SURFACES = [chassis, solar_array]