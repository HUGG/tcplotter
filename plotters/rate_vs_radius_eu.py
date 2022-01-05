#!/usr/bin/env python3

# Import libraries we need
import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np
import os
from pathlib import Path
from scipy.interpolate import interp1d
import subprocess


# Define function for calculating effective uranium concentration
def calc_eu(uranium, thorium):
    """Calculates effective uranium concentration from U, Th inputs"""
    return uranium + 0.235 * thorium


# Define function for creating plot of cooling rates
def rate_vs_radius_eu(n_inc=21, rate_min=0.1, rate_max=100.0, temp_max=250.0, ap_u_min=1.0, ap_u_max=150.0, zr_u_min=1.0, zr_u_max=4000.0,
                 ap_rad_min=40.0, ap_rad_max=100.0, zr_rad_min=40.0, zr_rad_max=100.0, ap_thorium=0.0, zr_thorium=0.0,
                 plot_type=3, save_plot=False, show_plot=True, verbose=False, use_widget=False):
    """
    A script for calculating thermochronometer ages and closure temperatures for different cooling rates, effective
    uranium concentrations, and equivalent spherical radii.

    Parameters
    ----------
    n_inc : int, default=21
        Number of points along x and y axes where ages/closure temperatures are
        calculated.
        NOTE: A value of n_inc = 101 was used in the manuscript. It has been
        reduced here to make the plotting faster. Set this to 101 to reproduce
        the manuscript Figure 4.
    rate_min : float, default=0.1
        Minimum cooling rate in degrees C per Myr.
    rate_max : float, default=100.0
        Maximum cooling rate in degrees C per Myr.
    temp_max : float, default=250.0
        Max temperature for cooling history (in degrees C).
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
        Cooling rate versus eU/radius.
        1 = apatite, 2 = zircon, 3 = both
    save_plot : bool, default=False
        Flag for whether to save the plot to a file.
    show_plot : bool, default=True
        Flag for whether to display the plot.
    verbose : bool, default=False
        Enable/disable verbose output.
    use_widget : bool, default=False
        Enable/disable IPython progress bar widget. Disabled for command-line usage.

    Returns
    -------
    None

    """
    # --- General model parameters ----------------------------------------------- #

    # Set base directory (use empty string for local working directory)
    fp = ''

    # --- Plotting parameters ---------------------------------------------------- #

    # Plotting flags and options
    dpi = 300
    out_fmt = 'pdf'

    # Set file name prefix
    plot_filename = 'rate_vs_radius_eu'

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

    # Create arrays of U concentrations
    ap_u = np.linspace(ap_u_min, ap_u_max, n_inc)
    zr_u = np.linspace(zr_u_min, zr_u_max, n_inc)

    # Create grain radius arrays
    ap_rad = np.linspace(ap_rad_min, ap_rad_max, n_inc)
    zr_rad = np.linspace(zr_rad_min, zr_rad_max, n_inc)

    # Create cooling rate array
    rates = np.logspace(start=np.log10(rate_min), stop=np.log10(rate_max), num=n_inc)

    # Calculate effective uranium
    ap_eu = calc_eu(ap_u, ap_thorium)
    zr_eu = calc_eu(zr_u, zr_thorium)

    # Total number of models
    total_models = len(ap_u) * len(rates) + len(ap_rad) * len(rates)

    # Screen output info
    if plot_type == 1:
        model_type = 'apatite age/Tc (cooling rate vs. radius/eU)'
    elif plot_type == 2:
        model_type = 'zircon age/Tc (cooling rate vs. radius/eU)'
    elif plot_type == 3:
        model_type = 'apatite/zircon age/Tc (cooling rate vs. radius/eU)'
    else:
        raise ValueError('Bad value for parameter plot_type. Must be 1, 2, or 3.')

    # Define time-temperature history filename
    tt_file = 'simple_time_temp.txt'

    # Check to make sure necessary age calculation executable(s) exist
    if not Path(fp + '../bin/RDAAM_He').is_file():
        raise FileNotFoundError("Age calculation program bin/RDAAM_He not found. Did you compile and install it?")

    # Create figure
    if plot_type < 3:
        fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    else:
        fig, ax = plt.subplots(2, 2, figsize=(10, 10))

    # Set plot loop variables
    ap_x1 = rates
    ap_y1 = ap_rad
    zr_x1 = rates
    zr_y1 = zr_rad
    ap_x2 = rates
    ap_y2 = ap_eu
    zr_x2 = rates
    zr_y2 = zr_eu

    # Create lists for storing closure temperatures, ages
    ahe_tc_list1 = []
    ahe_tc_list2 = []
    ap_x_list1 = []
    ap_y_list1 = []
    ap_x_list2 = []
    ap_y_list2 = []
    zhe_tc_list1 = []
    zhe_tc_list2 = []
    zr_x_list1 = []
    zr_y_list1 = []
    zr_x_list2 = []
    zr_y_list2 = []

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

    # Loop over plotables - loop 1: rate versus radius
    model_count = 0
    for i in range(len(ap_x1)):
        for j in range(len(ap_y1)):
            model_count += 1
            if not verbose:
                if use_widget:
                    s.value = model_count
                else:
                    print(
                        f'Calculating {model_type} - {int(round(100 * model_count / total_models)):3d}% ({model_count:5d} / {total_models:5d})\r',
                        end="")

            # Define parameters for this iteration
            rate = rates[i]
            ap_radius = ap_rad[j]
            zr_radius = zr_rad[j]
            ap_uranium = 10.0
            zr_uranium = 100.0
            ap_x_list1.append(rate)
            zr_x_list1.append(rate)
            ap_y_list1.append(ap_radius)
            zr_y_list1.append(zr_radius)

            # Write synthetic cooling history points to file
            start_time = temp_max / rate
            with open(tt_file, 'w') as f:
                f.write('0.0,0.0\n')
                f.write('{0:.4f},{1:.1f}'.format(start_time, temp_max))

            # Screen output
            if verbose:
                print(
                    f'Cooling from {temp_max:.1f}°C at a rate of {rate:.1f} °C/Myr will require {start_time:.2f} million years')

            # Calculate (U-Th)/He ages
            command = fp + '../bin/RDAAM_He ' + tt_file + ' ' + str(ap_radius) + ' ' + str(ap_uranium) + ' ' + str(
                ap_thorium) + ' ' + str(zr_radius) + ' ' + str(zr_uranium) + ' ' + str(zr_thorium)
            p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

            # Parse output for ages
            stdout = p.stdout.readlines()
            corr_ahe_age = stdout[0].split()[7].decode('UTF-8')
            corr_zhe_age = stdout[1].split()[7].decode('UTF-8')

            # Find closure temperatures from cooling ages and thermal history
            tc_interp = interp1d([0.0, start_time], [0.0, temp_max])
            ahe_tc = tc_interp(float(corr_ahe_age))
            zhe_tc = tc_interp(float(corr_zhe_age))

            # Add closure temperatures to lists
            ahe_tc_list1.append(ahe_tc)
            zhe_tc_list1.append(zhe_tc)

            if verbose:
                print(
                    f'AHe: {float(corr_ahe_age):.2f} Ma (Tc: {ahe_tc:.1f}°C); ZHe: {float(corr_zhe_age):.2f} Ma (Tc: {zhe_tc:.1f}°C)')

    # Loop over plotables - loop 2: rate versus eU
    for i in range(len(ap_x2)):
        for j in range(len(ap_y2)):
            model_count += 1
            if not verbose:
                if use_widget:
                    s.value = model_count
                else:
                    print(
                        f'Calculating {model_type} - {int(round(100 * (model_count) / total_models)):3d}% ({model_count:5d} / {total_models:5d})\r',
                        end="")

            # Define parameters for this iteration
            rate = rates[i]
            ap_radius = 45.0
            zr_radius = 60.0
            ap_uranium = ap_u[j]
            zr_uranium = zr_u[j]
            ap_x_list2.append(rate)
            zr_x_list2.append(rate)
            ap_y_list2.append(ap_uranium)
            zr_y_list2.append(zr_uranium)

            # Write synthetic cooling history points to file
            start_time = temp_max / rate
            with open(tt_file, 'w') as f:
                f.write('0.0,0.0\n')
                f.write('{0:.4f},{1:.1f}'.format(start_time, temp_max))

            # Screen output
            if verbose:
                print(
                    f'Cooling from {temp_max:.1f}°C at a rate of {rate:.1f} °C/Myr will require {start_time:.2f} million years')

            # Calculate (U-Th)/He ages
            command = fp + '../bin/RDAAM_He ' + tt_file + ' ' + str(ap_radius) + ' ' + str(ap_uranium) + ' ' + str(
                ap_thorium) + ' ' + str(zr_radius) + ' ' + str(zr_uranium) + ' ' + str(zr_thorium)
            p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

            # Parse output for ages
            stdout = p.stdout.readlines()
            corr_ahe_age = stdout[0].split()[7].decode('UTF-8')
            corr_zhe_age = stdout[1].split()[7].decode('UTF-8')

            # Find closure temperatures from cooling ages and thermal history
            tc_interp = interp1d([0.0, start_time], [0.0, temp_max])
            ahe_tc = tc_interp(float(corr_ahe_age))
            zhe_tc = tc_interp(float(corr_zhe_age))

            # Add closure temperatures to lists
            ahe_tc_list2.append(ahe_tc)
            zhe_tc_list2.append(zhe_tc)

            if verbose:
                print(
                    f'AHe: {float(corr_ahe_age):.2f} Ma (Tc: {ahe_tc:.1f}°C); ZHe: {float(corr_zhe_age):.2f} Ma (Tc: {zhe_tc:.1f}°C)')

    # Clean up temporary tt file
    os.remove(tt_file)

    # Plot only values for apatite (U-Th)/He
    if plot_type == 1:
        # --- Apatite cooling rate versus radius ---
        # Create closure temperature contour lines
        ap_contours_tc = ax[0].tricontour(ap_x_list1, ap_y_list1, ahe_tc_list1, n_cont, linewidths=0.5, colors='k')
        # Use log x-axis scaling
        ax[0].set_xscale('log')
        # Add closure temperature contour labels
        ax[0].clabel(ap_contours_tc, fmt='%1.1f')
        # Create closure temperature contour fill
        ap_contourf_tc1 = ax[0].tricontourf(ap_x_list1, ap_y_list1, ahe_tc_list1, n_fill_cont, cmap=colormap,
                                               alpha=colormap_alpha)

        # This is the fix for the white lines between contour levels
        for c in ap_contourf_tc1.collections:
            c.set_edgecolor("face")

        # --- Apatite cooling rate versus eU plot ---
        # Create closure temperature contour lines
        ap_contours_tc = ax[1].tricontour(ap_x_list2, ap_y_list2, ahe_tc_list2, n_cont, linewidths=0.5, colors='k')
        # Use log x-axis scaling
        ax[1].set_xscale('log')
        # Add closure temperature contour labels
        ax[1].clabel(ap_contours_tc, fmt='%1.1f')
        # Create closure temperature contour fill
        ap_contourf_tc2 = ax[1].tricontourf(ap_x_list2, ap_y_list2, ahe_tc_list2, n_fill_cont, cmap=colormap,
                                               alpha=colormap_alpha)

        # This is the fix for the white lines between contour levels
        for c in ap_contourf_tc2.collections:
            c.set_edgecolor("face")

    # Plot only values for zircon (U-Th)/He
    elif plot_type == 2:
        # --- Zircon cooling rate versus radius ---
        # Create closure temperature contour lines
        zr_contours_tc = ax[0].tricontour(zr_x_list1, zr_y_list1, zhe_tc_list1, n_cont, linewidths=0.5, colors='k')
        # Use log x-axis scaling
        ax[0].set_xscale('log')
        # Add closure temperature contour labels
        ax[0].clabel(zr_contours_tc, fmt='%1.1f')
        # Create closure temperature contour fill
        zr_contourf_tc1 = ax[0].tricontourf(zr_x_list1, zr_y_list1, zhe_tc_list1, n_fill_cont, cmap=colormap,
                                               alpha=colormap_alpha)

        # This is the fix for the white lines between contour levels
        for c in zr_contourf_tc1.collections:
            c.set_edgecolor("face")

        # --- Zircon cooling rate versus eU plot ---
        # Create closure temperature contour lines
        zr_contours_tc = ax[1].tricontour(zr_x_list2, zr_y_list2, zhe_tc_list2, n_cont, linewidths=0.5, colors='k')
        # Use log x-axis scaling
        ax[1].set_xscale('log')
        # Add closure temperature contour labels
        ax[1].clabel(zr_contours_tc, fmt='%1.1f')
        # Create closure temperature contour fill
        zr_contourf_tc2 = ax[1].tricontourf(zr_x_list2, zr_y_list2, zhe_tc_list2, n_fill_cont, cmap=colormap,
                                               alpha=colormap_alpha)

        # This is the fix for the white lines between contour levels
        for c in zr_contourf_tc2.collections:
            c.set_edgecolor("face")

    # Plot values for apatite and zircon (U-Th)/He
    else:
        # --- Apatite cooling rate versus radius ---
        # Create closure temperature contour lines
        ap_contours_tc = ax[0][0].tricontour(ap_x_list1, ap_y_list1, ahe_tc_list1, n_cont, linewidths=0.5, colors='k')
        # Use log x-axis scaling
        ax[0][0].set_xscale('log')
        # Add closure temperature contour labels
        ax[0][0].clabel(ap_contours_tc, fmt='%1.1f')
        # Create closure temperature contour fill
        ap_contourf_tc1 = ax[0][0].tricontourf(ap_x_list1, ap_y_list1, ahe_tc_list1, n_fill_cont, cmap=colormap,
                                               alpha=colormap_alpha)

        # This is the fix for the white lines between contour levels
        for c in ap_contourf_tc1.collections:
            c.set_edgecolor("face")

        # --- Apatite cooling rate versus eU plot ---
        # Create closure temperature contour lines
        ap_contours_tc = ax[0][1].tricontour(ap_x_list2, ap_y_list2, ahe_tc_list2, n_cont, linewidths=0.5, colors='k')
        # Use log x-axis scaling
        ax[0][1].set_xscale('log')
        # Add closure temperature contour labels
        ax[0][1].clabel(ap_contours_tc, fmt='%1.1f')
        # Create closure temperature contour fill
        ap_contourf_tc2 = ax[0][1].tricontourf(ap_x_list2, ap_y_list2, ahe_tc_list2, n_fill_cont, cmap=colormap,
                                               alpha=colormap_alpha)

        # This is the fix for the white lines between contour levels
        for c in ap_contourf_tc2.collections:
            c.set_edgecolor("face")

        # --- Zircon cooling rate versus radius plot ---
        # Create closure temperature contour lines
        zr_contours_tc = ax[1][0].tricontour(zr_x_list1, zr_y_list1, zhe_tc_list1, n_cont, linewidths=0.5, colors='k')
        # Use log x-axis scaling
        ax[1][0].set_xscale('log')
        # Add closure temperature contour labels
        ax[1][0].clabel(zr_contours_tc, fmt='%1.1f')
        # Create closure temperature contour fill
        zr_contourf_tc1 = ax[1][0].tricontourf(zr_x_list1, zr_y_list1, zhe_tc_list1, n_fill_cont, cmap=colormap,
                                               alpha=colormap_alpha)

        # This is the fix for the white lines between contour levels
        for c in zr_contourf_tc1.collections:
            c.set_edgecolor("face")

        # --- Zircon cooling rate versus eU plot ---
        # Create closure temperature contour lines
        zr_contours_tc = ax[1][1].tricontour(zr_x_list2, zr_y_list2, zhe_tc_list2, n_cont, linewidths=0.5, colors='k')
        # Use log x-axis scaling
        ax[1][1].set_xscale('log')
        # Add closure temperature contour labels
        ax[1][1].clabel(zr_contours_tc, fmt='%1.1f')
        # Create closure temperature contour fill
        zr_contourf_tc2 = ax[1][1].tricontourf(zr_x_list2, zr_y_list2, zhe_tc_list2, n_fill_cont, cmap=colormap,
                                               alpha=colormap_alpha)

        # This is the fix for the white lines between contour levels
        for c in zr_contourf_tc2.collections:
            c.set_edgecolor("face")

    # Format plot

    # Apatite only
    if plot_type == 1:
        ax[0].set_title('Apatite (U-Th)/He closure temperature [°C]')
        ax[1].set_title('Apatite (U-Th)/He closure temperature [°C]')
    # Zircon only
    elif plot_type == 2:
        ax[0].set_title('Zircon (U-Th)/He closure temperature [°C]')
        ax[1].set_title('Zircon (U-Th)/He closure temperature [°C]')
    # Apatite and zircon
    else:
        ax[0][0].set_title('Apatite (U-Th)/He closure temperature [°C]')
        ax[0][1].set_title('Apatite (U-Th)/He closure temperature [°C]')
        ax[1][0].set_title('Zircon (U-Th)/He closure temperature [°C]')
        ax[1][1].set_title('Zircon (U-Th)/He closure temperature [°C]')

    # Apatite or Zircon
    if plot_type < 3:
        ax[0].set_xlabel('Cooling rate [°C/Myr]')
        ax[1].set_xlabel('Cooling rate [°C/Myr]')
        ax[0].set_ylabel('Equivalent spherical radius (µm)')
        ax[1].set_ylabel('Effective uranium (ppm)')

    # Apatite and zircon eU versus radius
    else:
        ax[0][0].set_xlabel('Cooling rate [°C/Myr]')
        ax[0][1].set_xlabel('Cooling rate [°C/Myr]')
        ax[0][0].set_ylabel('Equivalent spherical radius (µm)')
        ax[0][1].set_ylabel('Effective uranium (ppm)')
        ax[1][0].set_xlabel('Cooling rate [°C/Myr]')
        ax[1][1].set_xlabel('Cooling rate [°C/Myr]')
        ax[1][0].set_ylabel('Equivalent spherical radius (µm)')
        ax[1][1].set_ylabel('Effective uranium (ppm)')

    # Don't use scientific notation for x-axis
    if plot_type < 3:
        ax[0].xaxis.set_major_formatter(ScalarFormatter())
        ax[1].xaxis.set_major_formatter(ScalarFormatter())
    else:
        ax[0][0].xaxis.set_major_formatter(ScalarFormatter())
        ax[0][1].xaxis.set_major_formatter(ScalarFormatter())
        ax[1][0].xaxis.set_major_formatter(ScalarFormatter())
        ax[1][1].xaxis.set_major_formatter(ScalarFormatter())

    # Use tight layout for subplots
    plt.tight_layout()

    # Save plot if requested
    if save_plot:
        if plot_type == 1:
            plot_savename = plot_filename + '_apatite_' + str(dpi) + 'dpi.' + out_fmt
        elif plot_type == 2:
            plot_savename = plot_filename + '_zircon_' + str(dpi) + 'dpi.' + out_fmt
        else:
            plot_savename = plot_filename + '_apatite_zircon_' + str(dpi) + 'dpi.' + out_fmt
        plt.savefig(fp + '../' + plot_savename, dpi=dpi)

    # Save plot if requested
    if show_plot:
        plt.show()

    return None


def main():
    parser = argparse.ArgumentParser(description='Calculates thermochronometer ages and closure temperatures for '
                                                 'different cooling rates, effective uranium concentrations, '
                                                 'and equivalent spherical radii.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--n_inc', help='Number of points along x and y axes where ages/closure temperatures are '
                                        'calculated. NOTE: A value of n_inc = 101 was used in the manuscript. It has '
                                        'been reduced here to make the plotting faster. Set this to 101 to reproduce '
                                        'the manuscript Figure 4.', default=21, type=int)
    parser.add_argument('--rate_min', help='Minimum cooling rate in degrees C per Myr.', default=0.1, type=float)
    parser.add_argument('--rate_max', help='Maximum cooling rate in degrees C per Myr.', default=100.0, type=float)
    parser.add_argument('--temp_max', help='Max temperature for cooling history (in degrees C).', default=250.0,
                        type=float)
    parser.add_argument('--ap_u_min', help='Minimum apatite uranium concentration in ppm', default=1.0, type=float)
    parser.add_argument('--ap_u_max', help='Maximum apatite uranium concentration in ppm', default=150.0, type=float)
    parser.add_argument('--zr_u_min', help='Minimum zircon uranium concentration in ppm', default=1.0, type=float)
    parser.add_argument('--zr_u_max', help='Maximum zircon uranium concentration in ppm', default=1500.0, type=float)
    parser.add_argument('--ap_rad_min', help='Minimum apatite equivalent spherical grain radius in micrometers', default=40.0, type=float)
    parser.add_argument('--ap_rad_max', help='Maximum apatite equivalent spherical grain radius in micrometers', default=100.0, type=float)
    parser.add_argument('--zr_rad_min', help='Minimum zircon equivalent spherical grain radius in micrometers', default=40.0, type=float)
    parser.add_argument('--zr_rad_max', help='Maximum zircon equivalent spherical grain radius in micrometers', default=100.0, type=float)
    parser.add_argument('--ap_thorium', help='Apatite thorium concentration in ppm', default=0.0, type=float)
    parser.add_argument('--zr_thorium', help='Zircon thorium concentration in ppm', default=0.0, type=float)
    parser.add_argument('--plot_type', help='Cooling rate versus radius/eU plot type. 1 = apatite, 2 = zircon, 3 = both.', default=3, type=int)
    parser.add_argument('--save_plot', help='Save plot to file?', default=False, type=bool)
    parser.add_argument('--show_plot', help='Display plot on the screen?', default=True, type=bool)
    parser.add_argument('--verbose', help='Enable/disable verbose output', default=False, type=bool)

    args = parser.parse_args()

    rate_vs_radius_eu(args.n_inc, args.rate_min, args.rate_max, args.temp_max, args.ap_u_min, args.ap_u_max, args.zr_u_min, args.zr_u_max, args.ap_rad_min, args.ap_rad_max, args.zr_rad_min, args.zr_rad_max, args.ap_thorium, args.zr_thorium, args.plot_type, args.save_plot, args.show_plot, args.verbose, use_widget=False)


if __name__ == "__main__":
    # execute only if run as a script
    main()