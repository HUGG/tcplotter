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
def rate_vs_age_tc(num_points=101, rate_min=0.1, rate_max=100.0, temp_max=250.0, ap_u1=1.0, ap_u2=20.0, ap_u3=150.0,
                   zr_u1=10.0, zr_u2=200.0, zr_u3=4000.0, ap_rad=45.0, zr_rad=60.0, ap_thorium=0.0, zr_thorium=0.0,
                   ahe_uncertainty=0.1, aft_uncertainty=0.2, zhe_uncertainty=0.1, plot_type=3, save_plot=False,
                   plot_file_format='pdf', display_plot=True, verbose=False, use_widget=False):
    """
    Calculates thermochronometer ages and closure temperatures for different cooling rates and effective uranium
    concentrations.

    Parameters
    ----------
    num_points : int, default=101
        Number of points along x and y axes where ages/closure temperatures are
        calculated.
    rate_min : float, default=0.1
        Minimum cooling rate in degrees C per Myr.
    rate_max : float, default=100.0
        Maximum cooling rate in degrees C per Myr.
    temp_max : float, default=250.0
        Max temperature for cooling history (in degrees C).
    ap_u1 : float, default=1.0
        Apatite uranium concentration in ppm for upper plot panel.
    ap_u2 : float, default=10.0
        Apatite uranium concentration in ppm for middle plot panel.
    ap_u3 : float, default=150.0
        Apatite uranium concentration in ppm for lower plot panel.
    zr_u1 : float, default=10.0
        Zircon uranium concentration in ppm for upper plot panel.
    zr_u2 : float, default=200.0
        Zircon uranium concentration in ppm for middle plot panel.
    zr_u3 : float, default=4000.0
        Zircon uranium concentration in ppm for lower plot panel.
    ap_rad : float, default=45.0
        Apatite equivalent spherical grain radius in micrometers.
    zr_rad : float, default=60.0
        Zircon equivalent spherical grain radius in micrometers.
    ap_thorium : float, default=0.0
        Apatite thorium concentration in ppm.
    zr_thorium : float, default=0.0
        Zircon thorium concentration in ppm.
    ahe_uncertainty : float, default=0.1
        Apatite (U-Th)/He age uncertainty fraction (0.1 = 10%)
    aft_uncertainty : float, default=0.2
        Apatite fission-track age uncertainty fraction (0.2 = 20%)
    zhe_uncertainty : float, default=0.1
        Zircon (U-Th)/He age uncertainty fraction (0.1 = 10%)
    plot_type : int, default=3
        1 = Cooling rate versus closure temperature
        2 = Cooling rate versus age
        3 = Cooling rate versus age and closure temperature
    save_plot : bool, default=False
        Flag for whether to save the plot to a file.
    plot_file_format : str, default='pdf'
        File format for saving plot to file (examples: png, pdf, svg, eps).
    display_plot : bool, default=True
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

    # Clean up temporary files?
    clean_up_files = True

    # --- Plotting parameters ---------------------------------------------------- #

    # Plotting flags and options
    dpi = 300

    # Set plot style
    plt.style.use('seaborn-darkgrid')

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

    # Make lists for apatite and zircon uranium concentrations
    ap_u_list = [ap_u1, ap_u2, ap_u3]
    zr_u_list = [zr_u1, zr_u2, zr_u3]

    # Set plot file name prefix
    if plot_type == 1:
        plot_filename = 'rate_vs_tc'
    elif plot_type == 2:
        plot_filename = 'rate_vs_age'
    elif plot_type == 3:
        plot_filename = 'rate_vs_age_tc'
    else:
        raise ValueError('Bad value for plot_type. Must be 1, 2, or 3.')

    # Define cooling rates to consider
    rates = np.logspace(start=np.log10(rate_min), stop=np.log10(rate_max), num=num_points)

    # Plot titles
    title_list = [f'Low eU (ap={ap_u_list[0]:.1f}, zr={zr_u_list[0]:.1f} ppm)',
                  f'Intermediate eU (ap={ap_u_list[1]:.1f}, zr={zr_u_list[1]:.1f} ppm)',
                  f'High eU (ap={ap_u_list[2]:.1f}, zr={zr_u_list[2]:.1f} ppm)']

    # Define time-temperature history filename
    tt_file = 'simple_time_temp.txt'

    # Check to make sure necessary age calculation executable(s) exist
    if not Path('../bin/RDAAM_He').is_file():
        raise FileNotFoundError("Age calculation program bin/RDAAM_He not found. Did you compile and install it?")
    if not Path('../bin/ketch_aft').is_file():
        raise FileNotFoundError("Age calculation program bin/ketch_aft not found. Did you compile and install it?")

    # Calculate total number of models that will be run
    total_models = len(ap_u_list) * len(rates)

    # Set model type string
    if plot_type == 1:
        model_type = 'cooling rate versus closure temperature'
    elif plot_type == 2:
        model_type = 'cooling rate versus age'
    elif plot_type == 3:
        model_type = 'cooling rate versus age and closure temperature'

    # Create figure
    if plot_type == 3:
        fig, ax = plt.subplots(3, 2, figsize=(12, 10))
    else:
        fig, ax = plt.subplots(3, 1, figsize=(6, 10))

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

    # Loop over plots/plot pairs
    model_count = 0
    for i in range(len(ap_u_list)):
        ap_uranium = ap_u_list[i]
        zr_uranium = zr_u_list[i]

        # Create lists for plotables
        rate_list = []
        ahe_tc_list = []
        aft_tc_list = []
        zhe_tc_list = []
        ahe_age_list = []
        aft_age_list = []
        zhe_age_list = []
        for rate in rates:
            model_count += 1
            if not verbose:
                if use_widget:
                    s.value = model_count
                else:
                    print(
                        f'Calculating {model_type} - {int(round(100 * model_count / total_models)):3d}% ({model_count:5d} / {total_models:5d})\r',
                        end="")

            # Define thermal history
            start_time = temp_max / rate
            with open(tt_file, 'w') as f:
                f.write('0.0,0.0\n')
                f.write('{0:.4f},{1:.1f}'.format(start_time, temp_max))

            # Calculate He ages
            command = '../bin/RDAAM_He ' + tt_file + ' ' + str(ap_rad) + ' ' + str(ap_uranium) + ' ' + str(
                ap_thorium) + ' ' + str(zr_rad) + ' ' + str(zr_uranium) + ' ' + str(zr_thorium)
            p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

            # Parse output for ages
            stdout = p.stdout.readlines()
            corr_ahe_age = stdout[0].split()[7].decode('UTF-8')
            corr_zhe_age = stdout[1].split()[7].decode('UTF-8')

            # Calculate AFT age
            command = '../bin/ketch_aft ' + tt_file
            p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

            # Parse output for AFT age
            stdout = p.stdout.readlines()
            aft_age = stdout[0].split()[4][:-1].decode('UTF-8')

            # Use predicted ages to get closure temperature
            tc_interp = interp1d([0.0, start_time], [0.0, temp_max])
            ahe_tc = tc_interp(float(corr_ahe_age))
            aft_tc = tc_interp(float(aft_age))
            zhe_tc = tc_interp(float(corr_zhe_age))

            # Add current iteration values to plotable lists
            rate_list.append(rate)
            ahe_tc_list.append(ahe_tc)
            aft_tc_list.append(aft_tc)
            zhe_tc_list.append(zhe_tc)
            ahe_age_list.append(float(corr_ahe_age))
            aft_age_list.append(float(aft_age))
            zhe_age_list.append(float(corr_zhe_age))

            # Echo ages for this iteration
            if verbose:
                print(
                    f'AHe: {float(corr_ahe_age):.2f} Ma (Tc: {ahe_tc:.1f}°C); AFT: {float(aft_age):.2f} Ma (Tc: {aft_tc:.1f}°C); ZHe: {float(corr_zhe_age):.2f} Ma (Tc: {zhe_tc:.1f}°C) -- total time: {start_time:.1f} Myr')

        # Assign uncertainties if plotting ages
        if plot_type != 1:
            # Calculate age min and max values using given uncertainties
            ahe_age_min = np.array(ahe_age_list) * (1.0 - ahe_uncertainty)
            ahe_age_max = np.array(ahe_age_list) * (1.0 + ahe_uncertainty)
            aft_age_min = np.array(aft_age_list) * (1.0 - aft_uncertainty)
            aft_age_max = np.array(aft_age_list) * (1.0 + aft_uncertainty)
            zhe_age_min = np.array(zhe_age_list) * (1.0 - zhe_uncertainty)
            zhe_age_max = np.array(zhe_age_list) * (1.0 + zhe_uncertainty)

        # Create plots for rate versus closure temperature
        if plot_type == 1:
            ax[i].semilogx(rate_list, ahe_tc_list, label='Apatite (U-Th)/He')
            ax[i].semilogx(rate_list, aft_tc_list, label='Apatite FT')
            ax[i].semilogx(rate_list, zhe_tc_list, label='Zircon (U-Th)/He')

        # Create plots for rate versus age
        if plot_type == 2:
            ax[i].fill_between(rate_list, ahe_age_min, ahe_age_max, alpha=0.5,
                               label=f'Apatite (U-Th)/He age ± {ahe_uncertainty * 100:.0f}%')
            ax[i].fill_between(rate_list, aft_age_min, aft_age_max, alpha=0.5,
                               label=f'Apatite FT age ± {aft_uncertainty * 100:.0f}%')
            ax[i].fill_between(rate_list, zhe_age_min, zhe_age_max, alpha=0.5,
                               label=f'Zircon (U-Th)/He age ± {zhe_uncertainty * 100:.0f}%')

            # Scale axes
            ax[i].set_xscale('log')
            ax[i].set_yscale('log')

        # Create plots for rate versus age and closure temperature
        if plot_type == 3:
            # Plot ages and closure temperatures (low eU)
            ax[i][0].fill_between(rate_list, ahe_age_min, ahe_age_max, alpha=0.5,
                                  label=f'Apatite (U-Th)/He age ± {ahe_uncertainty * 100:.0f}%')
            ax[i][1].plot(rate_list, ahe_tc_list, label='Apatite (U-Th)/He')

            # Plot ages and closure temperatures (intermediate eU)
            ax[i][0].fill_between(rate_list, aft_age_min, aft_age_max, alpha=0.5,
                                  label=f'Apatite FT age ± {aft_uncertainty * 100:.0f}%')
            ax[i][1].plot(rate_list, aft_tc_list, label='Apatite FT')

            # Plot ages and closure temperatures (high eU)
            ax[i][0].fill_between(rate_list, zhe_age_min, zhe_age_max, alpha=0.5,
                                  label=f'Zircon (U-Th)/He age ± {zhe_uncertainty * 100:.0f}%')
            ax[i][1].plot(rate_list, zhe_tc_list, label='Zircon (U-Th)/He')

            # Set axis scalings
            ax[i][0].set_xscale('log')
            ax[i][0].set_yscale('log')
            ax[i][1].set_xscale('log')

        # Format axis tick labels
        if plot_type == 3:
            ax[i][0].xaxis.set_major_formatter(ScalarFormatter())
            ax[i][1].xaxis.set_major_formatter(ScalarFormatter())
            ax[i][0].yaxis.set_major_formatter(ScalarFormatter())
        else:
            ax[i].xaxis.set_major_formatter(ScalarFormatter())
            ax[i].yaxis.set_major_formatter(ScalarFormatter())

        # Set axis range and add axis labels
        if plot_type == 1:
            ax[i].set_xlim([rate_min, rate_max])
            ax[i].set_ylabel('Closure temperature (°C)')
            if i == 2:
                ax[i].set_xlabel('Cooling rate (°C/Myr)')

        # Set axis range and add axis labels
        if plot_type == 2:
            ax[i].set_xlim([rate_min, rate_max])
            ax[i].set_ylabel('Age (Ma)')
            if i == 2:
                ax[i].set_xlabel('Cooling rate (°C/Myr)')

        # Set axis ranges and add axis labels
        if plot_type == 3:
            ax[i][0].set_xlim([rate_min, rate_max])
            ax[i][0].set_ylim([0.5, 1800.0])
            ax[i][1].set_xlim([rate_min, rate_max])
            ax[i][1].set_ylim([0.0, 200.0])
            ax[i][0].set_ylabel('Age (Ma)')
            ax[i][1].set_ylabel('Closure temperature (°C)')
            if i == 2:
                ax[i][0].set_xlabel('Cooling rate (°C/Myr)')
                ax[i][1].set_xlabel('Cooling rate (°C/Myr)')

        # Add subplot titles
        if plot_type == 3:
            ax[i][0].set_title(title_list[i])
            ax[i][1].set_title(title_list[i])
        else:
            ax[i].set_title(title_list[i])

        if plot_type == 3:
            ax[i][0].legend()
            ax[i][1].legend()
        else:
            ax[i].legend()

    # Delete temporary tt file
    if clean_up_files:
        os.remove(tt_file)
        os.remove('ft_length.csv')

    # Use tight layout
    plt.tight_layout()

    # Save plot if requested
    if save_plot:
        # Create plots directory if it does not already exist
        plots_exists = os.path.exists('../plots')
        if not plots_exists:
            # Create a new directory because it does not exist
            os.makedirs('../plots')

        # Define plot filename and save plot
        plot_filename = plot_filename + '_' + str(dpi) + 'dpi.' + plot_file_format
        plt.savefig('../plots/' + plot_filename, dpi=dpi)

    # Show plot if requested
    if display_plot:
        plt.show()

    # Revert to original working directory
    os.chdir(wd_orig)

    return None


def main():
    parser = argparse.ArgumentParser(description='Calculates thermochronometer ages and closure temperatures for '
                                                 'different cooling rates and eU concentrations.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--num-points', dest='num_points', help='Number of points along x and y axes where ages/closure temperatures are '
                                        'calculated.', default=101, type=int)
    parser.add_argument('--rate-min', dest='rate_min', help='Minimum cooling rate in degrees C per Myr.', default=0.1, type=float)
    parser.add_argument('--rate-max', dest='rate_max', help='Maximum cooling rate in degrees C per Myr.', default=100.0, type=float)
    parser.add_argument('--temp-max', dest='temp_max', help='Max temperature for cooling history (in degrees C).', default=250.0,
                        type=float)
    parser.add_argument('--ap-u1', dest='ap_u1', help='Apatite uranium concentration in ppm for upper plot panel', default=1.0,
                        type=float)
    parser.add_argument('--ap-u2', dest='ap_u2', help='Apatite uranium concentration in ppm for middle plot panel', default=20.0,
                        type=float)
    parser.add_argument('--ap-u3', dest='ap_u3', help='Apatite uranium concentration in ppm for lower plot panel', default=150.0,
                        type=float)
    parser.add_argument('--zr-u1', dest='zr_u1', help='Zircon uranium concentration in ppm for upper plot panel', default=10.0,
                        type=float)
    parser.add_argument('--zr-u2', dest='zr_u2', help='Zircon uranium concentration in ppm for middle plot panel', default=200.0,
                        type=float)
    parser.add_argument('--zr-u3', dest='zr_u3', help='Zircon uranium concentration in ppm for lower plot panel', default=4000.0,
                        type=float)
    parser.add_argument('--ap-rad', dest='ap_rad', help='Apatite equivalent spherical grain radius in micrometers',
                        default=45.0, type=float)
    parser.add_argument('--zr-rad', dest='zr_rad', help='Zircon equivalent spherical grain radius in micrometers',
                        default=60.0, type=float)
    parser.add_argument('--ap-thorium', dest='ap_thorium', help='Apatite thorium concentration in ppm', default=0.0, type=float)
    parser.add_argument('--zr-thorium', dest='zr_thorium', help='Zircon thorium concentration in ppm', default=0.0, type=float)
    parser.add_argument('--ahe-uncertainty', dest='ahe_uncertainty', help='Apatite U-Th/He age uncertainty fraction. Enter 0.1 for 10 percent '
                                                  'uncertainty.', default=0.1,
                        type=float)
    parser.add_argument('--aft-uncertainty', dest='aft_uncertainty', help='Apatite fission-track age uncertainty fraction. Enter 0.2 for 20 '
                                                  'percent uncertainty.',
                        default=0.2,
                        type=float)
    parser.add_argument('--zhe-uncertainty', dest='zhe_uncertainty', help='Zircon U-Th/He age uncertainty fraction. Enter 0.1 for 10 percent '
                                                  'uncertainty.', default=0.1,
                        type=float)
    parser.add_argument('--plot-type', dest='plot_type',
                        help='1 = Cooling rate versus closure temperature. 2 = Cooling rate versus age. 3 = Cooling '
                             'rate versus age and closure temperature',
                        default=3, type=int)
    parser.add_argument('--save-plot', dest='save_plot', help='Save plot to file', action='store_true',
                        default=False)
    parser.add_argument('--plot-file-format', dest='plot_file_format',
                        help='File format for saving plot to file (examples: png, pdf, svg, eps)', default='pdf',
                        type=str)
    parser.add_argument('--no-display-plot', dest='no_display_plot', help='Do not display plot on the screen',
                        action='store_true', default=False)
    parser.add_argument('-v', '--verbose', help='Enable verbose output', action='store_true', default=False)

    args = parser.parse_args()

    # Flip command-line flag to be opposite for function call
    # Function call expects display_plot = True for plot to be displayed
    display_plot = not args.no_display_plot

    rate_vs_age_tc(num_points=args.num_points, rate_min=args.rate_min, rate_max=args.rate_max, temp_max=args.temp_max,
                   ap_u1=args.ap_u1, ap_u2=args.ap_u2, ap_u3=args.ap_u3, zr_u1=args.zr_u1, zr_u2=args.zr_u2,
                   zr_u3=args.zr_u3, ap_rad=args.ap_rad, zr_rad=args.zr_rad, ap_thorium=args.ap_thorium,
                   zr_thorium=args.zr_thorium, ahe_uncertainty=args.ahe_uncertainty,
                   aft_uncertainty=args.aft_uncertainty, zhe_uncertainty=args.zhe_uncertainty, plot_type=args.plot_type,
                   save_plot=args.save_plot, plot_file_format=args.plot_file_format, display_plot=display_plot, verbose=args.verbose, use_widget=False)


if __name__ == "__main__":
    # execute only if run as a script
    main()
