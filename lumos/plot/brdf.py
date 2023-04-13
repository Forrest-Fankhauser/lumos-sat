import numpy as np
import lumos.conversions

def plot1D(ax, brdf_model, incident_angles = (10, 30, 50, 70), log_space = True):
    """
    Plots the in-plane BRDF at given incident angles.

    :param ax: Plotting axis
    :type ax: :class:`matplotlib.pyplot.axes`
    :param brdf_model: BRDF function to plot
    :type brdf_model: callable
    :param incident_angles: Incident angles for which to plot BRDF model (degrees)
    :type incident_angles: tuple, optional
    :param log_space: Whether to plot in log space
    :type log_space: bool, optional
    """
    for angle in incident_angles:
        x = np.linspace(-90, 90, 400)
        ix, iy, iz = lumos.conversions.spherical_to_unit(np.deg2rad(angle), np.pi)
        ox, oy, oz = lumos.conversions.spherical_to_unit(np.deg2rad(x), 0)
        y = brdf_model((ix, iy, iz), (0, 0, 1), (ox, oy, oz))
        label = r"$\phi_{in}$ = " + f"{angle:0.1f}°"
        if log_space:
            ax.semilogy(x, y, label = label)
        else:
            ax.plot(x, y, label = label)
    ax.legend()

def plot2D(polar_ax, brdf_model, incident_angle):
    """
    Plots the BRDF at given incident angle.

    :param ax: Plotting axis. Must have polar projection.
    :type ax: :class:`matplotlib.axes.Axes`
    :param brdf_model: BRDF function to plot
    :type brdf_model: callable
    :param incident_angle: Incident angle for which to plot BRDF model (degrees)
    :type incident_angle: float
    """
    phi_out = np.linspace(0, 90, 180)
    theta_out = np.linspace(0, 360, 360)
    phi_out, theta_out = np.meshgrid(phi_out, theta_out)

    ix, iy, iz = lumos.conversions.spherical_to_unit(np.deg2rad(incident_angle), np.pi)
    ox, oy, oz = lumos.conversions.spherical_to_unit(np.deg2rad(phi_out), np.deg2rad(theta_out))
    y = brdf_model((ix, iy, iz), (0, 0, 1), (ox, oy, oz))
    polar_ax.contourf(np.deg2rad(theta_out), phi_out, np.log10(y))
    polar_ax.set_title(r"$\phi_{in}$ = " + f"{incident_angle:0.1f}°")
    polar_ax.set_xticks([])
    polar_ax.set_yticks([])
    polar_ax.set_theta_zero_location('N')