#!/usr/bin/env python3

# Import libraries we need
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np
import argparse


# Define function for creating plot of cooling rates
def time_vs_temp(rate_min=0.1, rate_slow=1.0, rate_avg=10.0, rate_max=100.0, temp_max=250.0, time_max=50.0,
                 save_plot=False, show_plot=True,
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
    show_plot : bool, default=True
        Flag for whether to display the plot.

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

    # Set plot style
    plt.style.use('seaborn-whitegrid')

    # --- End of user-defined parameters ----------------------------------------- #
    #                                                                              #
    #  You probably don't need to modify anything below unless you know what you   #
    #  are doing :)                                                                #
    #                                                                              #
    # ---------------------------------------------------------------------------- #

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
        plot_filename = 'cooling_rates_figure1_' + str(dpi) + 'dpi.' + out_fmt
        # Save plots in base directory (not py subdirectory)
        plt.savefig(fp + '../' + plot_filename, dpi=dpi)

    # Show plot if requested
    if show_plot:
        plt.show()

    return None


# Define main function to support command-line usage with argparse
def main():
    parser = argparse.ArgumentParser(description='Plots cooling rate lines for different input rates',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--rate_min', help='Minimum cooling rate to plot in degrees C / Myr', default=0.1, type=float)
    parser.add_argument('--rate_slow', help='"Slow" cooling rate to plot in degrees C / Myr', default=1.0, type=float)
    parser.add_argument('--rate_avg', help='"Average" cooling rate to plot in degrees C / Myr', default=10.0,
                        type=float)
    parser.add_argument('--rate_max', help='Maximum cooling rate to plot in degrees C / Myr', default=100.0, type=float)
    parser.add_argument('--temp_max', help='Maximum temperature for cooling history in degrees C', default=250.0,
                        type=float)
    parser.add_argument('--time_max', help='Maximum value for time on x-axis of plot in millions of years ago (Ma)',
                        default=50.0,
                        type=float)
    parser.add_argument('--save_plot', help='Save plot to file?', default=False, type=bool)
    parser.add_argument('--show_plot', help='Display plot on the screen?', default=True, type=bool)

    args = parser.parse_args()

    time_vs_temp(args.rate_min, args.rate_slow, args.rate_avg, args.rate_max, args.temp_max, args.time_max,
                 args.save_plot, args.show_plot
                 )


if __name__ == "__main__":
    # Execute only if run as a script from the command line
    main()
