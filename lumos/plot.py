""" Helpful plotting functions for Lumos """

import matplotlib.pyplot as plt
import matplotlib.patches
import matplotlib.cm

import numpy as np

import lumos.conversions
import lumos.constants
import lumos.geometry
import lumos.calculator

import os
import cv2

def BRDF_1D(ax, brdf_model, incident_angles = (10, 30, 50, 70), log_space = True):
    """
    Plots the in-plane BRDF at given incident angles.

    :param ax: Plotting axis
    :type ax: :class:`matplotlib.pyplot.axes`
    :param brdf_model: BRDF function to plot
    :type brdf_model: function
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

def BRDF_2D(polar_ax, brdf_model, incident_angle):
    """
    Plots the BRDF at given incident angle.

    :param ax: Plotting axis. Must have polar projection.
    :type ax: :class:`matplotlib.axes.Axes`
    :param brdf_model: BRDF function to plot
    :type brdf_model: function
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

def contour_satellite_frame(
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

    cs = ax.contourf(observers.dists_on_plane / 1000,
                    observers.dists_off_plane / 1000,
                    values,
                    cmap = matplotlib.colormaps[cmap],
                    norm = matplotlib.colors.Normalize(levels[0], levels[1]),
                    levels = range(levels[0], levels[1] + 1),
                    extend = 'both')
    
    for collection in cs.collections:
        collection.set_clip_path(circ)

def mark_angles_above_horizon_satellite_frame(
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

        annotation_loc = (0.707 * R, -0.707 * R)
        y = np.where( (x - annotation_loc[0])**2 + (y - annotation_loc[1]) ** 2 < 200**2, np.inf, y)

        ax.plot(x, y, '--', color = "grey", linewidth = 0.8)
        ax.annotate(f'{np.rad2deg(angle):.0f}°', annotation_loc, fontsize = 8, color = 'grey', 
                horizontalalignment = 'center', verticalalignment = 'center')

def contour_observer_frame(
    ax,
    altitudes,
    azimuths,
    values,
    levels = None,
    cmap = 'plasma'
    ):

    """
    Creates contour plot in observer frame

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

    ax.contourf(
        np.deg2rad(azimuths),
        90 - altitudes,
        values,
        cmap = matplotlib.colormaps[cmap],
        norm = matplotlib.colors.Normalize(levels[0], levels[1]),
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

def colorbar(cax, levels):
    cmap = matplotlib.colormaps['plasma_r']
    norm = matplotlib.colors.Normalize(levels[0], levels[1])
    plt.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax = cax, extend = 'both')
    cax.set_aspect(2)

def mark_sun_azimuth_observer_frame(
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

def mark_sun_altitude_observer_frame(
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

def plot_compass(ax):
    ax.arrow(0, 0, 0, 0.75, width = 0.05, color = "white", head_length = 0.2, head_width = 0.2)
    ax.arrow(0, 0, -0.75, 0, width = 0.05, color = "white", head_length = 0.2, head_width = 0.2)
    ax.scatter([0], [0], s = 30, c = "white", zorder = 1)
    ax.annotate("N", (0, 1.15), horizontalalignment = 'center', 
                    verticalalignment = 'center', fontsize = 14, annotation_clip = False,
                    fontweight = 'bold')
    ax.annotate("E", (-1.15, 0), horizontalalignment = 'center', 
                    verticalalignment = 'center', fontsize = 14, annotation_clip = False,
                    fontweight = 'bold')
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim((-1.5, 0.5))
    ax.set_ylim((-0.5, 1.5))
    ax.set_frame_on(False)
    ax.set_aspect("equal")

def create_video(image_folder_path, video_output_path, frame_rate):
    """
    Combines folder of .png images to create a .mp4 video

    :param image_folder_path: Path to folder containing images.
    :type image_folder_path: str
    :param video_output_path: Destination path for output video
    :type video_output_path: str
    :param frame_rate: Video frames per second
    :type frame_rate: int
    """
    
    images = [img for img in sorted(os.listdir(image_folder_path)) if img.endswith(".png")]
    frame = cv2.imread(os.path.join(image_folder_path, images[0]))
    
    height, width, _ = frame.shape
    
    video = cv2.VideoWriter(video_output_path, 0, frame_rate, (width, height))
    
    for image in images:
        video.write(cv2.imread(os.path.join(image_folder_path, image)))
    
    cv2.destroyAllWindows()
    video.release()

def brightness_summary_satellite_frame(
        surfaces, angles_past_terminator, sat_height,
        include_sun = True, include_earthshine = False,
        earth_panel_density = 151, earth_brdf = None,
        levels = (0, 10)):
    
    N_frames = len(angles_past_terminator)

    with plt.style.context("dark_background"):
        fig, axs = plt.subplots(1, N_frames + 1,
                                width_ratios = N_frames * [1] + [0.1])

        for ax, angle_past_terminator in zip(axs[:-1], angles_past_terminator):

            observers = lumos.geometry.GroundObservers(
                sat_height,
                np.deg2rad(angle_past_terminator),
                density = 50
            )

            observers.calculate_intensity(surfaces)
            
            # Convert intensity to AB Magnitude
            ab_magnitudes = lumos.conversions.intensity_to_ab_mag(observers.intensities)

            # Plots intensity
            contour_satellite_frame(
                ax,
                observers,
                ab_magnitudes,
                levels = levels,
                cmap = 'plasma_r'
                )
            
            ax.set_title(r"$\alpha = $" + f"{angle_past_terminator:0.2f}°")

            mark_angles_above_horizon_satellite_frame(ax, observers)
        
        for ax in axs[1:-1]:
            ax.set_ylabel("")
            ax.set_yticklabels("")

        colorbar(axs[-1], levels)
        axs[-1].invert_yaxis()
        axs[-1].set_ylabel("AB Magnitude")
        axs[-1].set_aspect( 15 / (levels[1] - levels[0]))

        plt.show()

def brightness_summary_observer_frame(
        surfaces, sat_height, sun_altitudes, sun_azimuths,
        include_sun = True, include_earthshine = False,
        earth_panel_density = 151, earth_brdf = None,
        levels = (0, 10)):
    
    N_frames = len(sun_altitudes)

    altitudes = np.linspace(0, 90, 45)
    azimuths = np.linspace(0, 360, 90)
    altitudes, azimuths = np.meshgrid(altitudes, azimuths)

    with plt.style.context("dark_background"):
        fig = plt.figure(figsize = (6.4 * 2, 4.8 * 2))

        cax = fig.add_axes([0, 0.25, 0.15, 0.3])

        h = 0.8
        w = h * 6 / 32

        ax1 = fig.add_axes([0.5 - w/2 - 2 * (w + 0.0075), 0.025, w, h], projection = 'polar')
        ax2 = fig.add_axes([0.5 - w/2 - 1 * (w + 0.0075), 0.025, w, h], projection = 'polar')
        ax3 = fig.add_axes([0.5 - w/2 + 0 * (w + 0.0075), 0.025, w, h], projection = 'polar')
        ax4 = fig.add_axes([0.5 - w/2 + 1 * (w + 0.0075), 0.025, w, h], projection = 'polar')
        ax5 = fig.add_axes([0.5 - w/2 + 2 * (w + 0.0075), 0.025, w, h], projection = 'polar')
        ax6 = fig.add_axes([0.85, 0.05, 0.45 * w, 0.45 * h])

        axs = (ax1, ax2, ax3, ax4, ax5)

        plot_compass(ax6)

        for ax, sun_altitude, sun_azimuth in zip(axs, sun_altitudes, sun_azimuths):

            intensities = lumos.calculator.get_intensity_observer_frame(
                surfaces, sat_height, altitudes, azimuths,
                sun_altitude, sun_azimuth, 
                include_sun, include_earthshine, earth_panel_density, earth_brdf
            )
            
            # Convert intensity to AB Magnitude
            ab_magnitudes = lumos.conversions.intensity_to_ab_mag(intensities)

            # Plots intensity
            contour_observer_frame(ax, altitudes, azimuths, ab_magnitudes, levels, cmap = "plasma_r")

            ax.set_xticklabels(["", "", "", ""])
            ax.set_title(f"Sun Alt. = {sun_altitude:0.0f}°\nSun Az. = {sun_azimuth:0.0f}°")
            ax.grid(linewidth = 0.5, alpha = 0.25)
            
        colorbar(cax, levels = levels)
        cax.tick_params(labelsize = 14)
        cax.set_ylabel("AB Magnitude", fontsize = 16)
        cax.invert_yaxis()
        cax.yaxis.set_label_position("left")

        plt.show()