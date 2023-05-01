""" Plotting tools for calculations in the ground frame """

import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib.colors
import numpy as np

def plot_contour(
    ax,
    altitudes,
    azimuths,
    values,
    levels = None,
    cmap = 'plasma'
    ):

    """
    Creates contour plot in ground frame

    :param ax: Matplotlib axis for plotting on
    :type ax: :class:`matplotlib.pyplot.axes`
    :param altitudes: Altitudes in HCS frame (degrees)
    :type altitudes: :class:`np.ndarray`
    :param azimuths: Azimuths in HCS frame (degrees)
    :type azimuths: :class:`np.ndarray`
    :param values: Values to plot
    :type values: :class:`np.ndarray`
    :param levels: Minimum and maximum value to plot
    :type levels: tuple, optional
    :param cmap: Matplotlib colormap to use
    :type cmap: str
    """

    if levels is None:
        levels = (values.min(), values.max())

    cmap = matplotlib.colormaps[cmap]
    norm = matplotlib.colors.Normalize(levels[0], levels[1])

    ax.contourf(
        np.deg2rad(azimuths),
        90 - altitudes,
        values,
        cmap = cmap,
        norm = norm,
        levels = range(levels[0], levels[1] + 1),
        extend = 'both'
        )

    ax.set_rmax(90)
    ax.set_yticklabels([])
    ax.set_theta_zero_location('N')
    ax.set_rticks([10, 20, 30, 40, 50, 60, 70, 80, 90])
    ax.set_xticks(np.deg2rad([0, 90, 180, 270]))
    ax.set_xticklabels(['N', 'E', 'S', 'W'])
    ax.set_rlabel_position(-22.5)
    ax.grid(True)

def plot_colorbar(cax, levels, label):
    cmap = matplotlib.colormaps['plasma_r']
    norm = matplotlib.colors.Normalize(levels[0], levels[1])
    plt.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax = cax, extend = 'both')
    cax.invert_yaxis()
    cax.set_aspect(3)
    cax.set_ylabel(label)
    cax.yaxis.set_label_position("left")

def mark_sun_azimuth(
    ax,
    sun_azimuth,
    ):
    
    """
    Adds small sun on edge of plot

    :param ax: Axis for plotting
    :type ax: :class:`matplotlib.pyplot.axes`
    :param sun_azimuth: Azimuth of sun (degrees)
    :type sun_azimuth: float
    """
    ax.plot([np.deg2rad(sun_azimuth)], [95], 
            marker = (3, 0, 180 - sun_azimuth), 
            markersize = 6, 
            color = "white", 
            clip_on = False)
    
    ax.plot([np.deg2rad(sun_azimuth)], [101], 
            marker = '$\u2600$', 
            markersize = 10, 
            color = "orange", 
            clip_on = False)

def mark_sun_altitude(
    cax,
    sun_altitude
    ):

    """
    Adds colorbar with small sun to mark time of evening

    :param cax: Matplotlib axis for plotting
    :type cax: :class:`matplotlib.pyplot.axes`
    :param sun_altitude: Altitude of sun (degrees)
    :type sun_altitude: float
    """
    
    color_list = ('#15171c', '#20242d', '#25406a', '#4872bc', '#88a5d1', '#b5c9e6')
    evening_cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Evening', color_list)
    norm = matplotlib.colors.Normalize(vmin = -18, vmax = 0)

    plt.colorbar(matplotlib.cm.ScalarMappable(norm = norm, cmap = evening_cmap), cax = cax)
    cax.set_aspect('equal')
    cax.set_xlim(0, 1)
    cax.yaxis.set_tick_params(which = 'minor', length = 0, rotation = 0, pad = 35)
    cax.yaxis.set_tick_params(which = 'major', length = 5, pad = 5)
    cax.xaxis.set_label_position("top")
    cax.set_yticks(ticks = [-18, -12, -6, 0], labels = ['Night', '', '', 'Day'], fontsize = 12)
    cax.set_yticks(ticks = [-15.5, -9.5, -3.5], labels = ['Astronomical\nTwilight', 'Nautical\nTwilight', 'Civil\nTwilight'], minor = True, fontsize = 8, ma = 'center', ha = 'center')
    plot_alt = np.clip(sun_altitude, -18, 0)
    cax.plot([-0.4], [plot_alt], marker = (3, 0, 270), color = 'white', markersize = 6, clip_on = False)
    cax.plot([-1.0], [plot_alt], marker = '$\u2600$', color = 'orange', markersize = 10, clip_on = False)