""" Plotting AB Magnitude in ground frame """

import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib.gridspec
import numpy as np
import lumos.plot
import lumos.conversions as conversions
import scipy.ndimage

def plot_AB_Mag_contour(altitudes : np.ndarray, 
                        azimuths : np.ndarray, 
                        intensities : np.ndarray, 
                        sun_alt : float, 
                        sun_az : float,
                        levels : tuple[int, int] = (0, 8),
                        filter : bool = True, 
                        title : str = None, 
                        save_file : str = None):
    """
    Creates AB Magnitude plot in ground frame
        Parameters:
            altitudes (np.ndarray) : Altitudes (degrees)
            azimuths (np.ndarray) : Azimuths (degrees)
            intensities (np.ndarray) : Calculated intensities (W / m^2)
            sun_alt (float) : Altitude of sun (degrees)
            sun_az (float) : Azimuth of sun (degrees)
            levels (tuple[int, int]) : Maximum and minimum AB Magnitude to show on plot
            title (str) : Plot title
            save_file (str) : Path to output file
    """
    
    if filter:
        intensities = scipy.ndimage.gaussian_filter(intensities, 2, mode = ('wrap', 'reflect'))
    
    ab_mags = conversions.intensity_to_ab_mag(intensities)

    with plt.style.context(['dark_background', 'lumos.plot.ground_frame']):
        fig = plt.figure()

        gs = matplotlib.gridspec.GridSpec(1, 3, width_ratios = [1, 5, 1])
        cax1 = plt.subplot(gs[0])
        ax = plt.subplot(gs[1], projection = 'polar')
        cax2 = plt.subplot(gs[2])

        plt.subplots_adjust(top=0.75)

        # -------------- Central Heat Map ----------------------
        if title != None:
            fig.suptitle(title, x = 0.5, y = 0.95)
        
        ab_mags = np.clip(ab_mags, levels[0] + 0.01, levels[-1] - 0.01)

        cmap = matplotlib.colormaps['plasma_r']
        cs = ax.contourf(np.radians(azimuths), 90 - altitudes,
                         ab_mags, cmap = cmap, 
                         levels = range(levels[0], levels[-1] + 1))
        
        ax.plot([np.deg2rad(sun_az)], [95], marker = (3, 0, 180 - sun_az), markersize = 6, color = "white", clip_on = False)
        ax.plot([np.deg2rad(sun_az)], [101], marker = '$\u2600$', markersize = 10, color = "orange", clip_on = False)

        ax.set_rmax(90)
        ax.set_yticklabels([])
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.set_rticks([10, 20, 30, 40, 50, 60, 70, 80, 90])
        ax.set_xticks(np.deg2rad([0, 90, 180, 270]))
        ax.set_xticklabels(['N', 'E', 'S', 'W'])
        ax.set_rlabel_position(-22.5)
        ax.grid(True)

        plt.colorbar(cs, label = "AB Magnitude", cax = cax1, location = 'left')
        cax1.invert_yaxis()
        cax1.set_aspect(20)

        evening_cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Evening', ['#15171c', '#20242d', '#25406a', '#4872bc', '#88a5d1', '#b5c9e6'])
        norm = matplotlib.colors.Normalize(vmin=-18, vmax=0)
        plt.colorbar(matplotlib.cm.ScalarMappable(norm = norm, cmap = evening_cmap), cax = cax2)
        cax2.set_aspect('equal')
        cax2.set_xlim(0, 1)
        cax2.yaxis.set_tick_params(which = 'minor', length = 0, rotation = 0, pad = 35)
        cax2.yaxis.set_tick_params(which = 'major', length = 5, pad = 5)
        cax2.xaxis.set_label_position("top")
        cax2.set_yticks(ticks = [-18, -12, -6, 0], labels = ['Night', '', '', 'Day'], fontsize = 12)
        cax2.set_yticks(ticks = [-15.5, -9.5, -3.5], labels = ['Astronomical\nTwilight', 'Nautical\nTwilight', 'Civil\nTwilight'], minor = True, fontsize = 8, ma = 'center', ha = 'center')
        plot_alt = np.clip(sun_alt, -18, 0)
        cax2.plot([-0.4], [plot_alt], marker = (3, 0, 270), color = 'white', markersize = 6, clip_on = False)
        cax2.plot([-1.0], [plot_alt], marker = '$\u2600$', color = 'orange', markersize = 10, clip_on = False)
    
        if save_file == None:
            plt.show()
        else:
            fig.savefig(save_file)
            plt.close()

def plot_AB_Mag_scatter(altitudes : np.ndarray, 
                        azimuths : np.ndarray, 
                        intensities : np.ndarray, 
                        sun_alt : float, 
                        sun_az : float,
                        levels : tuple[int, int] = (0, 8),
                        title : str = None, 
                        save_file : str = None):
    """
    Creates AB Magnitude plot in ground frame
        Parameters:
            altitudes (np.ndarray) : Altitudes (degrees)
            azimuths (np.ndarray) : Azimuths (degrees)
            intensities (np.ndarray) : Calculated intensities (W / m^2)
            sun_alt (float) : Altitude of sun (degrees)
            sun_az (float) : Azimuth of sun (degrees)
            levels (tuple[int, int]) : Maximum and minimum AB Magnitude to show on plot
            title (str) : Plot title
            save_file (str) : Path to output file
    """
    
    ab_mags = conversions.intensity_to_ab_mag(intensities)

    with plt.style.context(['dark_background', 'lumos.plot.ground_frame']):
        fig = plt.figure()

        gs = matplotlib.gridspec.GridSpec(1, 3, width_ratios = [1, 5, 1])
        cax1 = plt.subplot(gs[0])
        ax = plt.subplot(gs[1], projection = 'polar')
        cax2 = plt.subplot(gs[2])

        plt.subplots_adjust(top=0.75)

        # -------------- Central Scatter Plot ----------------------
        if title != None:
            fig.suptitle(title, x = 0.5, y = 0.95)
        
        ab_mags = np.clip(ab_mags, levels[0] + 0.01, levels[-1] - 0.01)

        cmap = matplotlib.colormaps['plasma_r']
        cs = ax.scatter(np.radians(azimuths), 90 - altitudes,
            c = ab_mags, cmap = cmap, 
            norm = matplotlib.colors.Normalize(vmin = levels[0], vmax = levels[-1]))
        
        ax.plot([np.deg2rad(sun_az)], [95], marker = (3, 0, 180 - sun_az), markersize = 6, color = "white", clip_on = False)
        ax.plot([np.deg2rad(sun_az)], [101], marker = '$\u2600$', markersize = 10, color = "orange", clip_on = False)

        ax.set_rmax(90)
        ax.set_yticklabels([])
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.set_rticks([10, 20, 30, 40, 50, 60, 70, 80, 90])
        ax.set_xticks(np.deg2rad([0, 90, 180, 270]))
        ax.set_xticklabels(['N', 'E', 'S', 'W'])
        ax.set_rlabel_position(-22.5)
        ax.grid(True)

        plt.colorbar(cs, label = "AB Magnitude", cax = cax1, location = 'left')
        cax1.invert_yaxis()
        cax1.set_xlim((0, 1))
        cax1.set_aspect(18 / (levels[1] - levels[0]))

        evening_cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Evening', ['#15171c', '#20242d', '#25406a', '#4872bc', '#88a5d1', '#b5c9e6'])
        norm = matplotlib.colors.Normalize(vmin=-18, vmax=0)
        plt.colorbar(matplotlib.cm.ScalarMappable(norm = norm, cmap = evening_cmap), cax = cax2)
        cax2.set_aspect('equal')
        cax2.set_xlim(0, 1)
        cax2.yaxis.set_tick_params(which = 'minor', length = 0, rotation = 0, pad = 35)
        cax2.yaxis.set_tick_params(which = 'major', length = 5, pad = 5)
        cax2.xaxis.set_label_position("top")
        cax2.set_yticks(ticks = [-18, -12, -6, 0], labels = ['Night', '', '', 'Day'], fontsize = 12)
        cax2.set_yticks(ticks = [-15.5, -9.5, -3.5], labels = ['Astronomical\nTwilight', 'Nautical\nTwilight', 'Civil\nTwilight'], minor = True, fontsize = 8, ma = 'center', ha = 'center')
        plot_alt = np.clip(sun_alt, -18, 0)
        cax2.plot([-0.4], [plot_alt], marker = (3, 0, 270), color = 'white', markersize = 6, clip_on = False)
        cax2.plot([-1.0], [plot_alt], marker = '$\u2600$', color = 'orange', markersize = 10, clip_on = False)
    
        if save_file == None:
            plt.show()
        else:
            fig.savefig(save_file)
            plt.close()