from lumos.geometry import Surface
from lumos.brdf.library import ABG
import numpy as np

area = 5.0     # meters^2
normal = np.array([0, 0, -1])
brdf = ABG(A = 3.362e-03, B = 4.309e-06, g = 2.068)

aluminum_plate = Surface(area, normal, brdf)

SURFACES = [aluminum_plate]