#!/usr/bin/env python3

# Import libraries we need
import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np
import os


# Define function for creating plot of cooling rates
def time_vs_temp(rate_min=0.1, rate_slow=1.0, rate_avg=10.0, rate_max=100.0, temp_max=250.0, time_max=50.0,
                 save_plot=False, plot_file_format='pdf', display_plot=True,
                 ):
    """
    Plots cooling rate lines for different input rates.

    Parameters
    ----------
    rate_min : float or int, default=0.1
        Minimum cooling rate to plot in degrees C / Myr.
    rate_slow : float or int, default=1.0
        "Slow" cooling rate to plot in degrees C / Myr.
    rate_avg : float or int, default=10.0
        "Average" cooling rate to plot in degrees C / Myr.
    rate_max : float or int, default=100.0
        Maximum cooling rate to plot in degrees C / Myr.
    temp_max : float or int, default=250.0
        Maximum temperature for cooling history in degrees C.
    time_max : float or int, default=50.0
        Maximum value for time on x-axis of plot in millions of years ago (Ma).
    save_plot : bool, default=False
        Flag for whether to save the plot to a file.
    plot_file_format : str, default='pdf'
        File format for saving plot to file (examples: png, pdf, svg, eps).
    display_plot : bool, default=True
        Flag for whether to display the plot.

    Returns
    -------
    None
    """

    # --- Plotting parameters ---------------------------------------------------- #
    # Plotting flags and options
    dpi = 300

    # Set plot style
    plt.style.use('seaborn-whitegrid')

    # --- End of user-defined parameters ----------------------------------------- #
    #                                                                              #
    #  You probably don't need to modify anything below unless you know what you   #
    #  are doing :)                                                                #
    #                                                                              #
    # ---------------------------------------------------------------------------- #

    # Ensure relative paths work by setting working dir to dir containing this script file
    wd_orig = os.getcwd()
    script_path = os.path.abspath(__file__)
    dir_name = os.path.dirname(script_path)
    os.chdir(dir_name)

    # Create arrays of points to plot
    min_rate_x = np.array([temp_max / rate_min, 0.0])
    min_rate_y = np.array([temp_max, 0.0])
    slow_rate_x = np.array([temp_max / rate_slow, 0.0])
    slow_rate_y = np.array([temp_max, 0.0])
    avg_rate_x = np.array([temp_max / rate_avg, 0.0])
    avg_rate_y = np.array([temp_max, 0.0])
    max_rate_x = np.array([temp_max / rate_max, 0.0])
    max_rate_y = np.array([temp_max, 0.0])

    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))

    # Plot lines and fill
    ax.fill_betweenx(min_rate_y, min_rate_x, max_rate_x, color='black', alpha=0.15,
                     label='Range of model cooling rates')
    ax.plot(min_rate_x, min_rate_y, color='black')
    ax.plot(slow_rate_x, slow_rate_y, color='black')
    ax.plot(avg_rate_x, avg_rate_y, color='black')
    ax.plot(max_rate_x, max_rate_y, color='black')

    # Set axis tick label format
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_major_formatter(ScalarFormatter())

    # Set plot x and y range
    ax.set_xlim([0.0, time_max])
    ax.set_ylim([0.0, temp_max])

    # Add axis labels
    ax.set_xlabel('Time (Ma)')
    ax.set_ylabel('Temperature (°C)')

    # Flip axis directions
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()

    # Use tight layout
    plt.tight_layout()

    # Save plot if requested
    if save_plot:
        # Create plots directory if it does not already exist
        plots_exists = os.path.exists('../plots')
        if not plots_exists:
            # Create a new directory because it does not exist
            os.makedirs('../plots')

        # Set plot filename and save plot
        plot_filename = 'time_vs_temp_' + str(dpi) + 'dpi.' + plot_file_format
        plt.savefig('../plots/' + plot_filename, dpi=dpi)

    # Display plot if requested
    if display_plot:
        plt.show()

    # Revert to original working directory
    os.chdir(wd_orig)

    return None


# Define main function to support command-line usage with argparse
def main():
    parser = argparse.ArgumentParser(description='Plots cooling rate lines for different input rates',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--rate-min', dest='rate_min', help='Minimum cooling rate to plot in degrees C / Myr', default=0.1, type=float)
    parser.add_argument('--rate-slow', dest='rate_slow', help='"Slow" cooling rate to plot in degrees C / Myr', default=1.0, type=float)
    parser.add_argument('--rate-avg', dest='rate_avg', help='"Average" cooling rate to plot in degrees C / Myr', default=10.0,
                        type=float)
    parser.add_argument('--rate-max', dest='rate_max', help='Maximum cooling rate to plot in degrees C / Myr', default=100.0, type=float)
    parser.add_argument('--temp-max', dest='temp_max', help='Maximum temperature for cooling history in degrees C', default=250.0,
                        type=float)
    parser.add_argument('--time-max', dest='time_max', help='Maximum value for time on x-axis of plot in millions of years ago (Ma)',
                        default=50.0,
                        type=float)
    parser.add_argument('--save-plot', dest='save_plot', help='Save plot to file', action='store_true', default=False)
    parser.add_argument('--plot-file-format', dest='plot_file_format', help='File format for saving plot to file (examples: png, pdf, svg, eps)', default='pdf', type=str)
    parser.add_argument('--no-display-plot', dest='no_display_plot', help='Do not display plot on the screen', action='store_true', default=False)

    args = parser.parse_args()

    # Flip command-line flag to be opposite for function call
    # Function call expects display_plot = True for plot to be displayed
    display_plot = not args.no_display_plot

    time_vs_temp(rate_min=args.rate_min, rate_slow=args.rate_slow, rate_avg=args.rate_avg, rate_max=args.rate_max, temp_max=args.temp_max, time_max=args.time_max,
                 save_plot=args.save_plot, plot_file_format=args.plot_file_format, display_plot=display_plot
                 )


if __name__ == "__main__":
    # Execute only if run as a script from the command line
    main()
