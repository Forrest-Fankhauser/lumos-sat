""" Plotting tools for calculations in brightness frame """

import matplotlib.pyplot as plt
import matplotlib.patches
import matplotlib.cm

import numpy as np

import lumos.constants
import lumos.plot
import lumos.conversions
import lumos.geometry

def plot_contour(
    ax,
    observers,
    values,
    cmap = 'plasma',
    levels = None
    ):

    """
    Makes a contour plot

    :param ax: Plotting axis
    :type ax: :class:`matplotlib.pyplot.axes`
    :param observers: Mesh of observers
    :type observers: :class:`lumos.geometry.GroundObservers`
    :param values: 2D array of values to plot
    :type values: :class:`np.ndarray`
    :param cmap: Matplotlib colormap to make contour with
    :type cmap: str, optional
    :param levels: Minimum and maximum value to plot
    :type levels: tuple, optional
    """

    if levels is None:
        levels = (values.min(), values.max())

    max_dist = observers.max_angle * lumos.constants.EARTH_RADIUS / 1000

    ax.set_xlim(-max_dist, max_dist)
    ax.set_ylim(-max_dist, max_dist)
    ax.set_aspect("equal")
    
    ax.set_xlabel("Distance on Plane (km)")
    ax.set_ylabel("Distance off Plane (km)")

    circ = matplotlib.patches.Circle((0, 0), max_dist , transform=ax.transData, facecolor = (1, 1, 1))
    ax.add_patch(circ)

    cmap = matplotlib.colormaps[cmap]
    norm = matplotlib.colors.Normalize(levels[0], levels[1])

    cs = ax.contourf(observers.dists_on_plane / 1000,
                    observers.dists_off_plane / 1000,
                    values,
                    cmap = cmap,
                    norm = norm,
                    levels = range(levels[0], levels[1] + 1),
                    extend = 'both')
    
    for collection in cs.collections:
        collection.set_clip_path(circ)

def mark_angles_above_horizon(
    ax,
    observers,
    angles_above_horizon = (15, 30) 
    ):

    """
    Marks discrete angles above horizon using dashed circles

    :param ax: Matplotlib axis for plotting
    :type ax: :class:`matplotlib.pyplot.axes`
    :param observers: Mesh of observers
    :type observers: :class:`lumos.geometry.GroundObservers`
    :param angles_above_horizon: Angles to mark (degrees)
    :type angles_above_horizon: tuple, optional
    """

    theta = np.linspace(0, 2 * np.pi, 50)

    for angle in angles_above_horizon:
        angle = np.deg2rad(angle)
        R = (lumos.constants.EARTH_RADIUS / 1000
        * (np.pi/2 - angle - np.arcsin(np.cos(observers.max_angle) * np.cos(angle)))
        )

        x = R * np.cos(theta)
        y = R * np.sin(theta)

        annotation_loc = (0.707 * R, 0.707 * R)
        y = np.where( (x - annotation_loc[0])**2 + (y - annotation_loc[1]) ** 2 < 100**2, np.inf, y)

        ax.plot(x, y, '--k')
        ax.annotate(f'{np.rad2deg(angle):.0f}°', annotation_loc, fontsize = 10, color = 'black', 
                horizontalalignment = 'center', verticalalignment = 'center')