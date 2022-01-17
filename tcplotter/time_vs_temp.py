#!/usr/bin/env python3

# Import libraries we need
import argparse
from tcplotter import time_vs_temp


# Define main function to support command-line usage with argparse
def main():
    """Command-line interface for the time_vs_temp function in tcplotter."""
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
    parser.add_argument('--plot-dpi', dest='plot_dpi', help='Saved plot resolution in dots per inch', default=300, type=int)
    parser.add_argument('--plot-style', dest='plot_style', help='Style sheet used for plotting. See https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html.', default='seaborn-whitegrid',
                        type=str)
    parser.add_argument('--no-display-plot', dest='no_display_plot', help='Do not display plot on the screen', action='store_true', default=False)

    args = parser.parse_args()

    # Flip command-line flag to be opposite for function call
    # Function call expects display_plot = True for plot to be displayed
    display_plot = not args.no_display_plot

    time_vs_temp(rate_min=args.rate_min, rate_slow=args.rate_slow, rate_avg=args.rate_avg, rate_max=args.rate_max, temp_max=args.temp_max, time_max=args.time_max,
                 save_plot=args.save_plot, plot_file_format=args.plot_file_format, plot_dpi=args.plot_dpi, plot_style=args.plot_style, display_plot=display_plot
                 )


if __name__ == "__main__":
    # Execute only if run as a script from the command line
    main()
