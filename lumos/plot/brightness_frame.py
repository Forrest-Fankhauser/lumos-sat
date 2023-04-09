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
    intensities,
    title = None,
    levels = (0, 8)
    ):

    ab_mags = lumos.conversions.intensity_to_ab_mag(intensities)
    ab_mags = np.clip(ab_mags, levels[0] + 0.01, levels[-1] - 0.01)

    max_dist = observers.max_angle * lumos.constants.EARTH_RADIUS / 1000

    ax.set_xlim(-max_dist, max_dist)
    ax.set_ylim(-max_dist, max_dist)
    ax.set_aspect("equal")
    
    ax.set_xlabel("Distance on Plane (km)")
    ax.set_ylabel("Distance off Plane (km)")

    if title != None:
        ax.set_title(title, pad = 15)

    circ = matplotlib.patches.Circle((0, 0), max_dist , transform=ax.transData, facecolor = (1, 1, 1))
    ax.add_patch(circ)

    cmap = matplotlib.cm.get_cmap('plasma_r')
    cs = ax.contourf(observers.dists_on_plane / 1000,
                    observers.dists_off_plane / 1000,
                    ab_mags,
                    cmap = cmap,
                    levels = range(levels[0], levels[-1] + 1))
    
    for collection in cs.collections:
        collection.set_clip_path(circ)

    angles_above_horizon = (15, 30)
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
        
def plot_AB_Mag(observers : lumos.geometry.GroundObservers, 
                intensities : np.ndarray, 
                title : str = None, 
                save_file : str = None, 
                levels : tuple = (0, 8)):
    """
    Creates AB Magnitude plot in brightness frame
        Parameters:
            observers (lumos.geometry.GroundObservers) : Mesh of observers on ground
            intensities (np.ndarray) : Calculated intensities (W / m^2)
            title (str) : Plot title
            save_file (str) : Path to output file
            levels (tuple[int, int]) : Maximum and minimum AB Magnitude to show on plot
    """
    ab_mags = lumos.conversions.intensity_to_ab_mag(intensities)

    with plt.style.context(['dark_background', 'lumos.plot.brightness_frame']):
        fig, ax = plt.subplots()

        max_dist = observers.max_angle * lumos.constants.EARTH_RADIUS / 1000

        ax.set_xlim(-max_dist, max_dist)
        ax.set_ylim(-max_dist, max_dist)
        ax.set_aspect("equal")
    
        ax.set_xlabel("Distance on Plane (km)")
        ax.set_ylabel("Distance off Plane (km)")

        if title != None:
            ax.set_title(title, pad = 15)

        circ = matplotlib.patches.Circle((0, 0), max_dist , transform=ax.transData, facecolor = (1, 1, 1))
        ax.add_patch(circ)

        cmap = matplotlib.cm.get_cmap('plasma_r')

        contour_levels = range(levels[0], levels[-1] + 1)
        ab_mags = np.clip(ab_mags, levels[0] + 0.01, levels[-1] - 0.01)

        cs = ax.contourf(observers.dists_on_plane / 1000,
                        observers.dists_off_plane / 1000,
                        ab_mags,
                        cmap = cmap,
                        levels = contour_levels)
        
        for collection in cs.collections:
            collection.set_clip_path(circ)

        ax.invert_yaxis()
        cbar = fig.colorbar(cs, label = "AB Magnitude")
        cbar.ax.invert_yaxis()

        angles_above_horizon = (15, 30)
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
    
        if save_file == None:
            plt.show()
        else:
            fig.savefig(save_file)
            plt.close()