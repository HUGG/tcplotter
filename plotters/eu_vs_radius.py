#!/usr/bin/env python3

# Import libraries we need
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
from scipy.interpolate import interp1d
import subprocess


# Define function for calculating effective uranium concentration
def calc_eu(uranium, thorium):
    """Calculates effective uranium concentration from U, Th inputs"""
    return uranium + 0.235 * thorium


# Define function for making contour plot of cooling ages and closure temperatures
def eu_vs_radius(num_points=21, cooling_hist_type=1, temp_max=250.0, rate=10.0, time_hist=[0.0, 10.0, 25.0],
                 temp_hist=[0.0, 200.0, 250.0], ap_u_min=1.0, ap_u_max=150.0, zr_u_min=1.0, zr_u_max=4000.0,
                 ap_rad_min=40.0, ap_rad_max=100.0, zr_rad_min=40.0, zr_rad_max=100.0, ap_thorium=0.0, zr_thorium=0.0,
                 plot_type=3, save_plot=False, plot_file_format='pdf', display_plot=True, tt_plot=False, verbose=False, use_widget=False):
    """
    Calculates thermochronometer ages and closure temperatures for different effective uranium concentrations and
    equivalent spherical radii.

    Parameters
    ----------
    num_points : int, default=21
        Number of points along x and y axes where ages/closure temperatures are
        calculated.
        NOTE: A value of num_points = 101 was used in the manuscript. It has been
        reduced here to make the plotting faster. Set this to 101 to reproduce
        the manuscript Figures 2 or 3.
    cooling_hist_type : int, default=1
        Cooling history type.
        1 = constant cooling rate (specify rate as parameter rate)
        2 = list of time-temperature points (fill in lists as parameters
        time_hist, temp_hist)
    # Options for cooling history type 1
    temp_max : float, default=250.0
        Max temperature for cooling history (in degrees C).
    rate : float, default=10.0
        Cooling rate in degrees C per Myr.
    # Options for cooling history type 2
    time_hist : list of floats or ints, default=[0.0, 10.0, 25.0]
        Time points defining cooling history in Ma (millions of years ago).
        NOTE: Present-day point should be first in list.
    temp_hist : list of floats or ints, default=[0.0, 200.0, 250.0]
        Temperature points defining cooling history in degrees C.
        NOTE: Present-day point should be first in list.
    # End cooling history options
    ap_u_min : float, default=1.0
        Minimum apatite uranium concentration in ppm.
    ap_u_max : float, default=150.0
        Maximum apatite uranium concentration in ppm.
    zr_u_min : float, default=1.0
        Minimum zircon uranium concentration in ppm.
    zr_u_max : float, default=1500.0
        Maximum zircon uranium concentration in ppm.
    ap_rad_min : float, default=40.0
        Minimum apatite equivalent spherical grain radius in micrometers.
    ap_rad_max : float, default=100.0
        Maximum apatite equivalent spherical grain radius in micrometers.
    zr_rad_min : float, default=40.0
        Minimum zircon equivalent spherical grain radius in micrometers.
    zr_rad_max : float, default=100.0
        Maximum zircon equivalent spherical grain radius in micrometers.
    ap_thorium : float, default=0.0
        Apatite thorium concentration in ppm.
    zr_thorium : float, default=0.0
        Zircon thorium concentration in ppm.
    plot_type : int, default=3
        eU versus radius plot type.
        1 = apatite, 2 = zircon, 3 = both
    save_plot : bool, default=False
        Flag for whether to save the plot to a file.
    plot_file_format : str, default='pdf'
        File format for saving plot(s) to file (examples: png, pdf, svg, eps).
    display_plot : bool, default=True
        Flag for whether to display the plot.
    tt_plot : bool, default=False
        Flag for whether to create/display the time-temperature history plot.
    verbose : bool, default=False
        Enable/disable verbose output.
    use_widget : bool, default=False
        Enable/disable IPython progress bar widget. Disabled for command-line usage.
    """

    # --- Plotting parameters ---------------------------------------------------- #
    # Plotting flags and options
    dpi = 300

    # Set file name prefix
    plot_filename = 'eu_vs_radius'

    # Set plot style
    plt.style.use('seaborn-colorblind')
    colormap = 'plasma'
    colormap_alpha = 1.0

    # Number of plot contour lines
    n_cont = 12

    # Number of plot contour colors
    n_fill_cont = 256

    # --- End of user-defined parameters ----------------------------------------- #
    #                                                                              #
    #  You probably don't need to modify anything below unless you know what you   #
    #  are doing :)                                                                #
    #                                                                              #
    # ---------------------------------------------------------------------------- #

    # Check to see whether ipywidgets and IPython are available for widget use
    # If not, disable widgets and display a warning
    try:
        import ipywidgets as widgets
    except ModuleNotFoundError:
        print("Warning: ipywidgets module not found. Disabling graphical progress bar.")
        use_widget = False
    try:
        from IPython.display import display
    except ModuleNotFoundError:
        print("Warning: IPython.display module not found. Disabling graphical progress bar.")
        use_widget = False

    # Ensure relative paths work by setting working dir to dir containing this script file
    wd_orig = os.getcwd()
    script_path = os.path.abspath(__file__)
    dir_name = os.path.dirname(script_path)
    os.chdir(dir_name)

    # Define cooling history using constant cooling rate
    if cooling_hist_type == 1:
        # Define time and temperature histories
        start_time = temp_max / rate
        time_hist = [0.0, start_time]
        temp_hist = [0.0, temp_max]

    # Option 2: Define time-temperature history using list of tT points
    elif cooling_hist_type == 2:
        pass

    # Raise error if an unsupported value is given for cooling_hist_type
    else:
        raise ValueError('Bad value for cooling_hist_type. Should be 1 or 2.')

    # Create arrays of U concentrations
    ap_u = np.linspace(ap_u_min, ap_u_max, num_points)
    zr_u = np.linspace(zr_u_min, zr_u_max, num_points)

    # Create grain radius arrays
    ap_rad = np.linspace(ap_rad_min, ap_rad_max, num_points)
    zr_rad = np.linspace(zr_rad_min, zr_rad_max, num_points)

    # Calculate effective uranium
    ap_eu = calc_eu(ap_u, ap_thorium)
    zr_eu = calc_eu(zr_u, zr_thorium)

    # Calculate total number of models
    total_models = len(ap_u) * len(ap_rad)

    # Screen output info
    if plot_type == 1:
        model_type = 'apatite age/Tc (eU vs. radius)'
    elif plot_type == 2:
        model_type = 'zircon age/Tc (eU vs. radius)'
    elif plot_type == 3:
        model_type = 'apatite/zircon age/Tc (eU vs. radius)'
    else:
        raise ValueError('Bad value for plot_type. Should be 1, 2, or 3.')

    # Define time-temperature history filename
    tt_file = 'simple_time_temp.txt'

    # Check to make sure necessary age calculation executable(s) exist
    if not Path('../bin/RDAAM_He').is_file():
        raise FileNotFoundError("Age calculation program bin/RDAAM_He not found. Did you compile and install it?")

    # Create figure
    if plot_type < 3:
        fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    else:
        fig, ax = plt.subplots(2, 2, figsize=(10, 10))

    # Set plot loop variables
    ap_x = ap_eu
    ap_y = ap_rad
    zr_x = zr_eu
    zr_y = zr_rad

    # Create lists for storing closure temperatures, ages
    ahe_tc_list = []
    ahe_age_list = []
    ap_x_list = []
    ap_y_list = []
    zhe_tc_list = []
    zhe_age_list = []
    zr_x_list = []
    zr_y_list = []

    # Write cooling history points to file
    with open(tt_file, 'w') as f:
        for i in range(len(time_hist)):
            f.write(f'{time_hist[i]:.4f},{temp_hist[i]:.1f}\n')

    # Echo total model run time and cooling rate
    if verbose and cooling_hist_type == 1:
        print(
            f'Cooling from {temp_max:.1f}°C at a rate of {rate:.1f} °C/Myr will require {start_time:.2f} million years')

    # Create visual progress bar, if enabled
    if use_widget:
        s = widgets.IntProgress(
            value=0,
            min=0,
            max=total_models,
            description='Calculating:',
            bar_style='', # 'success', 'info', 'warning', 'danger' or ''
            style={'bar_color': '#ff6666'},
            orientation='horizontal'
        )
        display(s)

    # Loop over plotables
    model_count = 0
    for i in range(len(ap_x)):
        for j in range(len(ap_y)):
            model_count += 1
            if not verbose:
                if use_widget:
                    s.value = model_count
                else:
                    print(
                        f'Calculating {model_type} - {int(round(100 * model_count / total_models)):3d}% ({model_count:5d} / {total_models:5d})\r',
                        end="")

            # Define parameters for this iteration
            ap_uranium = ap_u[i]
            zr_uranium = zr_u[i]
            ap_radius = ap_rad[j]
            zr_radius = zr_rad[j]
            ap_x_list.append(ap_uranium)
            zr_x_list.append(zr_uranium)
            ap_y_list.append(ap_radius)
            zr_y_list.append(zr_radius)

            # Calculate (U-Th)/He ages
            command = '../bin/RDAAM_He ' + tt_file + ' ' + str(ap_radius) + ' ' + str(ap_uranium) + ' ' + str(
                ap_thorium) + ' ' + str(zr_radius) + ' ' + str(zr_uranium) + ' ' + str(zr_thorium)
            p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

            # Parse output for ages
            stdout = p.stdout.readlines()
            corr_ahe_age = stdout[0].split()[7].decode('UTF-8')
            corr_zhe_age = stdout[1].split()[7].decode('UTF-8')

            # Find closure temperatures from cooling ages and thermal history
            tc_interp = interp1d(time_hist, temp_hist)
            ahe_tc = tc_interp(float(corr_ahe_age))
            zhe_tc = tc_interp(float(corr_zhe_age))

            # Add closure temperatures, ages to lists
            ahe_tc_list.append(ahe_tc)
            ahe_age_list.append(float(corr_ahe_age))
            zhe_tc_list.append(zhe_tc)
            zhe_age_list.append(float(corr_zhe_age))

            if verbose:
                print(
                    f'AHe: {float(corr_ahe_age):.2f} Ma (Tc: {ahe_tc:.1f}°C); ZHe: {float(corr_zhe_age):.2f} Ma (Tc: {zhe_tc:.1f}°C)')

    # Clean up Tt file
    os.remove(tt_file)

    # Apatite eU versus radius
    if plot_type == 1:
        # Create age contour lines
        ap_contours_age = ax[0].tricontour(ap_x_list, ap_y_list, ahe_age_list, n_cont, linewidths=0.5, colors='k')
        # Add age contour labels
        ax[0].clabel(ap_contours_age, fmt='%1.1f')
        # Create age contour fill
        ap_contourf_age = ax[0].tricontourf(ap_x_list, ap_y_list, ahe_age_list, n_fill_cont, cmap=colormap,
                                            alpha=colormap_alpha)

        # This is the fix for the white lines between contour levels
        for c in ap_contourf_age.collections:
            c.set_edgecolor("face")

        # Create closure temperature contour lines
        ap_contours_tc = ax[1].tricontour(ap_x_list, ap_y_list, ahe_tc_list, n_cont, linewidths=0.5, colors='k')
        # Add closure temperature contour labels
        ax[1].clabel(ap_contours_tc, fmt='%1.1f')
        # Create closure temperature contour fill
        ap_contourf_tc = ax[1].tricontourf(ap_x_list, ap_y_list, ahe_tc_list, n_fill_cont, cmap=colormap,
                                           alpha=colormap_alpha)

        # This is the fix for the white lines between contour levels
        for c in ap_contourf_tc.collections:
            c.set_edgecolor("face")

    # Zircon eU versus radius
    elif plot_type == 2:
        # Create age contour lines
        zr_contours_age = ax[0].tricontour(zr_x_list, zr_y_list, zhe_age_list, n_cont, linewidths=0.5, colors='k')
        # Add age contour labels
        ax[0].clabel(zr_contours_age, fmt='%1.1f')
        # Create age contour fill
        zr_contourf_age = ax[0].tricontourf(zr_x_list, zr_y_list, zhe_age_list, n_fill_cont, cmap=colormap,
                                            alpha=colormap_alpha)

        # This is the fix for the white lines between contour levels
        for c in zr_contourf_age.collections:
            c.set_edgecolor("face")

        # Create closure temperature contour lines
        zr_contours_tc = ax[1].tricontour(zr_x_list, zr_y_list, zhe_tc_list, n_cont, linewidths=0.5, colors='k')
        # Add closure temperature contour labels
        ax[1].clabel(zr_contours_tc, fmt='%1.1f')
        # Create closure temperature contour fill
        zr_contourf_tc = ax[1].tricontourf(zr_x_list, zr_y_list, zhe_tc_list, n_fill_cont, cmap=colormap,
                                           alpha=colormap_alpha)

        # This is the fix for the white lines between contour levels
        for c in zr_contourf_tc.collections:
            c.set_edgecolor("face")

    # Apatite and zircon eU versus radius
    else:
        # Create age contour lines
        ap_contours_age = ax[0][0].tricontour(ap_x_list, ap_y_list, ahe_age_list, n_cont, linewidths=0.5, colors='k')
        # Add age contour labels
        ax[0][0].clabel(ap_contours_age, fmt='%1.1f')
        # Create age contour fill
        ap_contourf_age = ax[0][0].tricontourf(ap_x_list, ap_y_list, ahe_age_list, n_fill_cont, cmap=colormap,
                                               alpha=colormap_alpha)

        # This is the fix for the white lines between contour levels
        for c in ap_contourf_age.collections:
            c.set_edgecolor("face")

        # Create closure temperature contour lines
        ap_contours_tc = ax[0][1].tricontour(ap_x_list, ap_y_list, ahe_tc_list, n_cont, linewidths=0.5, colors='k')
        # Add closure temperature contour labels
        ax[0][1].clabel(ap_contours_tc, fmt='%1.1f')
        # Create closure temperature contour fill
        ap_contourf_tc = ax[0][1].tricontourf(ap_x_list, ap_y_list, ahe_tc_list, n_fill_cont, cmap=colormap,
                                              alpha=colormap_alpha)

        # This is the fix for the white lines between contour levels
        for c in ap_contourf_tc.collections:
            c.set_edgecolor("face")

        # Create age contour lines
        zr_contours_age = ax[1][0].tricontour(zr_x_list, zr_y_list, zhe_age_list, n_cont, linewidths=0.5, colors='k')
        # Add age contour labels
        ax[1][0].clabel(zr_contours_age, fmt='%1.1f')
        # Create age contour fill
        zr_contourf_age = ax[1][0].tricontourf(zr_x_list, zr_y_list, zhe_age_list, n_fill_cont, cmap=colormap,
                                               alpha=colormap_alpha)

        # This is the fix for the white lines between contour levels
        for c in zr_contourf_age.collections:
            c.set_edgecolor("face")

        # Create closure temperature contour lines
        zr_contours_tc = ax[1][1].tricontour(zr_x_list, zr_y_list, zhe_tc_list, n_cont, linewidths=0.5, colors='k')
        # Add closure temperature contour labels
        ax[1][1].clabel(zr_contours_tc, fmt='%1.1f')
        # Create closure temperature contour fill
        zr_contourf_tc = ax[1][1].tricontourf(zr_x_list, zr_y_list, zhe_tc_list, n_fill_cont, cmap=colormap,
                                              alpha=colormap_alpha)

        # This is the fix for the white lines between contour levels
        for c in zr_contourf_tc.collections:
            c.set_edgecolor("face")

    # Format plot

    # Apatite eU versus radius
    if plot_type == 1:
        ax[0].set_title('Apatite (U-Th)/He age [Ma]')
        ax[1].set_title('Apatite (U-Th)/He closure temperature [°C]')

    # Zircon eU versus radius
    elif plot_type == 2:
        ax[0].set_title('Zircon (U-Th)/He age [Ma]')
        ax[1].set_title('Zircon (U-Th)/He closure temperature [°C]')

    # Apatite and zircon eU versus radius
    else:
        ax[0][0].set_title('Apatite (U-Th)/He age [Ma]')
        ax[0][1].set_title('Apatite (U-Th)/He closure temperature [°C]')
        ax[1][0].set_title('Zircon (U-Th)/He age [Ma]')
        ax[1][1].set_title('Zircon (U-Th)/He closure temperature [°C]')

    # Apatite or Zircon eU versus radius
    if plot_type < 3:
        ax[0].set_xlabel('Effective uranium (ppm)')
        ax[1].set_xlabel('Effective uranium (ppm)')
        ax[0].set_ylabel('Equivalent spherical radius (µm)')
        ax[1].set_ylabel('Equivalent spherical radius (µm)')

    # Apatite and zircon eU versus radius
    else:
        ax[0][0].set_xlabel('Effective uranium (ppm)')
        ax[0][1].set_xlabel('Effective uranium (ppm)')
        ax[0][0].set_ylabel('Equivalent spherical radius (µm)')
        ax[0][1].set_ylabel('Equivalent spherical radius (µm)')
        ax[1][0].set_xlabel('Effective uranium (ppm)')
        ax[1][1].set_xlabel('Effective uranium (ppm)')
        ax[1][0].set_ylabel('Equivalent spherical radius (µm)')
        ax[1][1].set_ylabel('Equivalent spherical radius (µm)')

    # Use tight layout for subplots
    plt.tight_layout()

    # Save plot if desired
    if save_plot:
        # Create plots directory if it does not already exist
        plots_exists = os.path.exists('../plots')
        if not plots_exists:
            # Create a new directory because it does not exist
            os.makedirs('../plots')

        # Define plot filename based on type of plot and save plot
        if plot_type == 1:
            plot_savename = plot_filename + '_apatite_' + str(dpi) + 'dpi.' + plot_file_format
        elif plot_type == 2:
            plot_savename = plot_filename + '_zircon_' + str(dpi) + 'dpi.' + plot_file_format
        else:
            plot_savename = plot_filename + '_apatite_zircon_' + str(dpi) + 'dpi.' + plot_file_format
        plt.savefig('../plots/' + plot_savename, dpi=dpi)

    # Display plot if desired
    if display_plot:
        plt.show()

    # Create tT history plot if requested
    if tt_plot:
        # Create figure 2
        fig2, ax2 = plt.subplots(1, 1, figsize=(6, 5))

        # Plot tT history
        ax2.plot(time_hist, temp_hist, color='black')

        # Set plot x and y range
        ax2.set_xlim([0.0, max(time_hist)])
        ax2.set_ylim([0.0, max(temp_hist)])

        # Add axis labels
        ax2.set_xlabel('Time (Ma)')
        ax2.set_ylabel('Temperature (°C)')

        # Add title
        ax2.set_title('Time-temperature history')

        # Flip axis directions
        plt.gca().invert_xaxis()
        plt.gca().invert_yaxis()

        # Use tight layout
        plt.tight_layout()

        # Save plot if desired
        if save_plot:
            # Define plot filename and save plot
            plot_savename2 = plot_filename + '_tT_history_' + str(dpi) + 'dpi.' + plot_file_format
            plt.savefig('../plots/' + plot_savename2, dpi=dpi)

        # Display plot if desired
        if display_plot:
            plt.show()

    # Revert to original working directory
    os.chdir(wd_orig)

    return None


def main():
    parser = argparse.ArgumentParser(description='Calculates thermochronometer ages and closure temperatures for '
                                                 'different effective uranium concentrations and equivalent spherical '
                                                 'radii.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--num-points', dest='num_points', help='Number of points along x and y axes where ages/closure temperatures are '
                                        'calculated. NOTE: A value of num_points = 101 was used in the manuscript. It has '
                                        'been reduced here to make the plotting faster. Set this to 101 to reproduce '
                                        'the manuscript Figures 2 or 3.', default=21, type=int)
    parser.add_argument('--cooling-hist-type', dest='cooling_hist_type',
                        help='Cooling history type. 1 = constant cooling rate (specify rate as parameter rate). 2 = '
                             'list of time-temperature points (fill in lists as parameters time_hist, temp_hist).',
                        default=1, type=int)
    parser.add_argument('--temp-max', dest='temp_max', help='Max temperature for cooling history (in degrees C). Option for cooling '
                                           'history type 1.', default=250.0,
                        type=float)
    parser.add_argument('--rate', help='Cooling rate in degrees C per Myr. Option for cooling history type 1.', default=10.0, type=float)
    parser.add_argument('--time-hist', dest='time_hist', help='Time points defining cooling history in Ma (millions of years ago). '
                                            'NOTE: Present-day point should be first in list. Option for cooling '
                                            'history type 2.', nargs='+', default=[0.0, 10.0, 25.0],
                        type=float)
    parser.add_argument('--temp-hist', dest='temp_hist', help='Temperature points defining cooling history in degrees C. NOTE: '
                                            'Present-day point should be first in list. Option for cooling history '
                                            'type 2.', nargs='+', default=[0.0, 200.0, 250.0],
                        type=float)
    parser.add_argument('--ap-u-min', dest='ap_u_min', help='Minimum apatite uranium concentration in ppm', default=1.0, type=float)
    parser.add_argument('--ap-u-max', dest='ap_u_max', help='Maximum apatite uranium concentration in ppm', default=150.0, type=float)
    parser.add_argument('--zr-u-min', dest='zr_u_min', help='Minimum zircon uranium concentration in ppm', default=1.0, type=float)
    parser.add_argument('--zr-u-max', dest='zr_u_max', help='Maximum zircon uranium concentration in ppm', default=1500.0, type=float)
    parser.add_argument('--ap-rad-min', dest='ap_rad_min', help='Minimum apatite equivalent spherical grain radius in micrometers', default=40.0, type=float)
    parser.add_argument('--ap-rad-max', dest='ap_rad_max', help='Maximum apatite equivalent spherical grain radius in micrometers', default=100.0, type=float)
    parser.add_argument('--zr-rad-min', dest='zr_rad_min', help='Minimum zircon equivalent spherical grain radius in micrometers', default=40.0, type=float)
    parser.add_argument('--zr-rad-max', dest='zr_rad_max', help='Maximum zircon equivalent spherical grain radius in micrometers', default=100.0, type=float)
    parser.add_argument('--ap-thorium', dest='ap_thorium', help='Apatite thorium concentration in ppm', default=0.0, type=float)
    parser.add_argument('--zr-thorium', dest='zr_thorium', help='Zircon thorium concentration in ppm', default=0.0, type=float)
    parser.add_argument('--plot-type', dest='plot_type', help='eU versus radius plot type. 1 = apatite, 2 = zircon, 3 = both.', default=3, type=int)
    parser.add_argument('--save-plot', dest='save_plot', help='Save plot(s) to file', action='store_true',
                        default=False)
    parser.add_argument('--plot-file-format', dest='plot_file_format',
                        help='File format for saving plot(s) to file (examples: png, pdf, svg, eps)', default='pdf',
                        type=str)
    parser.add_argument('--no-display-plot', dest='no_display_plot', help='Do not display plots on the screen',
                        action='store_true', default=False)
    parser.add_argument('--tt-plot', dest='tt_plot', help='Create (and display) time-temperature history plot',
                        action='store_true', default=False)
    parser.add_argument('-v', '--verbose', help='Enable verbose output', action='store_true', default=False)

    args = parser.parse_args()

    # Flip command-line flag to be opposite for function call
    # Function call expects display_plot = True for plot to be displayed
    display_plot = not args.no_display_plot

    eu_vs_radius(num_points=args.num_points, cooling_hist_type=args.cooling_hist_type, temp_max=args.temp_max,
                 rate=args.rate, time_hist=args.time_hist, temp_hist=args.temp_hist, ap_u_min=args.ap_u_min,
                 ap_u_max=args.ap_u_max, zr_u_min=args.zr_u_min, zr_u_max=args.zr_u_max, ap_rad_min=args.ap_rad_min,
                 ap_rad_max=args.ap_rad_max, zr_rad_min=args.zr_rad_min, zr_rad_max=args.zr_rad_max,
                 ap_thorium=args.ap_thorium, zr_thorium=args.zr_thorium, plot_type=args.plot_type,
                 save_plot=args.save_plot, plot_file_format=args.plot_file_format, display_plot=display_plot,
                 tt_plot=args.tt_plot, verbose=args.verbose, use_widget=False)


if __name__ == "__main__":
    # execute only if run as a script
    main()
