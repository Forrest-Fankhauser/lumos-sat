import simple_sat
import lumos.plot
import lumos.brdf.library
import numpy as np

sat_height = 550 * 1000 # Geodetic height of satellite (meters)

angles_past_terminator = [5, 10, 15, 20]

lumos.plot.brightness_summary_satellite_frame(
    simple_sat.SURFACES,
    angles_past_terminator,
    sat_height,
    levels = (3, 10)
    )


# Angle of rotation of the surface
angle = np.deg2rad(20)

# Rotate the surface in the x direction
simple_sat.SURFACES[0].normal = np.array([np.sin(angle), 0, -np.cos(angle)])

lumos.plot.brightness_summary_satellite_frame(
    simple_sat.SURFACES,
    angles_past_terminator,
    sat_height,
    levels = (3, 10)
    )

# Angle of rotation of the surface
angle = np.deg2rad(45)

# Rotate the surface in the y direction
simple_sat.SURFACES[0].normal = np.array([0, np.sin(angle), -np.cos(angle)])

lumos.plot.brightness_summary_satellite_frame(
    simple_sat.SURFACES,
    angles_past_terminator,
    sat_height,
    levels = (3, 10)
    )


import lumos.geometry
observers = lumos.geometry.GroundObservers(
                sat_height,
                np.deg2rad(20),
                density = 200
            )

observers.calculate_intensity(simple_sat.SURFACES)

idx = np.argmax( observers.intensities )
dist_off_plane = observers.dists_off_plane.flatten()[idx]
dist_on_plane = observers.dists_on_plane.flatten()[idx]

print( f"Peak Intensity: {observers.intensities.max() * 1e9:0.0f} nW / m^2")
print( f"Distance off Plane: {dist_off_plane / 1000 :0.0f} km")
print( f"Distance on Plane: {dist_on_plane / 1000 :0.0f} km")