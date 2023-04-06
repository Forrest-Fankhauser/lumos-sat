import numpy as np
import matplotlib.pyplot as plt
import lumos.brdf.library as brdf_library
import lumos.conversions

def plot1D(model, incident_angles = (10, 30, 50, 70)):
    with plt.style.context(["dark_background"]):
        for angle in incident_angles:
            x = np.linspace(-90, 90, 400)
            ix, iy, iz = lumos.conversions.spherical_to_unit(np.deg2rad(angle), np.pi)
            ox, oy, oz = lumos.conversions.spherical_to_unit(np.deg2rad(x), 0)
            y = model((ix, iy, iz), (0, 0, 1), (ox, oy, oz))
            plt.semilogy(x, y, '--k')
        plt.legend()
        plt.show()

def plot2D(model, incident_angles = (10, 30, 50, 70)):
    fig, axs = plt.subplots(ncols = len(incident_angles), subplot_kw={'projection': 'polar'})
    with plt.style.context(["dark_background"]):
        for ax, aoi in zip(axs, incident_angles):
            phi_out = np.linspace(0, 90, 180)
            theta_out = np.linspace(0, 360, 360)
            phi_out, theta_out = np.meshgrid(phi_out, theta_out)

            ix, iy, iz = lumos.conversions.spherical_to_unit(np.deg2rad(aoi), np.pi)
            ox, oy, oz = lumos.conversions.spherical_to_unit(np.deg2rad(phi_out), np.deg2rad(theta_out))
            y = model((ix, iy, iz), (0, 0, 1), (ox, oy, oz))
            ax.contourf(np.deg2rad(theta_out), phi_out, np.log10(y))
            ax.set_title(r"$\phi_{in}$ = " + f"{aoi:0.1f}°")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_theta_zero_location('N')
        plt.show()