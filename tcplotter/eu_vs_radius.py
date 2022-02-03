#!/usr/bin/env python3

# Import libraries we need
import argparse
from tcplotter import eu_vs_radius


def main():
    """Command-line interface for the eu_vs_radius function in tcplotter."""
    parser = argparse.ArgumentParser(
        description="Calculates thermochronometer ages and closure temperatures for "
        "different effective uranium concentrations and equivalent spherical "
        "radii.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--num-points",
        dest="num_points",
        help="Number of points along x and y axes where ages/closure temperatures are "
        "calculated. NOTE: A value of num_points = 101 was used in the manuscript. It has "
        "been reduced here to make the plotting faster. Set this to 101 to reproduce "
        "the manuscript Figures 2 or 3.",
        default=21,
        type=int,
    )
    parser.add_argument(
        "--cooling-hist-type",
        dest="cooling_hist_type",
        help="Cooling history type. 1 = constant cooling rate (specify rate as parameter rate). 2 = "
        "list of time-temperature points (fill in lists as parameters time_hist, temp_hist).",
        default=1,
        type=int,
    )
    parser.add_argument(
        "--temp-max",
        dest="temp_max",
        help="Max temperature for cooling history (in degrees C). Option for cooling "
        "history type 1.",
        default=350.0,
        type=float,
    )
    parser.add_argument(
        "--cooling-rate",
        dest="cooling_rate",
        help="Cooling rate in degrees C per Myr. Option for cooling history type 1.",
        default=10.0,
        type=float,
    )
    parser.add_argument(
        "--time-hist",
        dest="time_hist",
        help="Time points defining cooling history in Ma (millions of years ago). "
        "NOTE: Present-day point should be first in list. Option for cooling "
        "history type 2.",
        nargs="+",
        default=[0.0, 10.0, 25.0, 35.0],
        type=float,
    )
    parser.add_argument(
        "--temp-hist",
        dest="temp_hist",
        help="Temperature points defining cooling history in degrees C. NOTE: "
        "Present-day point should be first in list. Option for cooling history "
        "type 2.",
        nargs="+",
        default=[0.0, 75.0, 50.0, 350.0],
        type=float,
    )
    parser.add_argument(
        "--ap-u-min",
        dest="ap_u_min",
        help="Minimum apatite uranium concentration in ppm",
        default=1.0,
        type=float,
    )
    parser.add_argument(
        "--ap-u-max",
        dest="ap_u_max",
        help="Maximum apatite uranium concentration in ppm",
        default=150.0,
        type=float,
    )
    parser.add_argument(
        "--zr-u-min",
        dest="zr_u_min",
        help="Minimum zircon uranium concentration in ppm",
        default=1.0,
        type=float,
    )
    parser.add_argument(
        "--zr-u-max",
        dest="zr_u_max",
        help="Maximum zircon uranium concentration in ppm",
        default=4000.0,
        type=float,
    )
    parser.add_argument(
        "--ap-rad-min",
        dest="ap_rad_min",
        help="Minimum apatite equivalent spherical grain radius in micrometers",
        default=40.0,
        type=float,
    )
    parser.add_argument(
        "--ap-rad-max",
        dest="ap_rad_max",
        help="Maximum apatite equivalent spherical grain radius in micrometers",
        default=100.0,
        type=float,
    )
    parser.add_argument(
        "--zr-rad-min",
        dest="zr_rad_min",
        help="Minimum zircon equivalent spherical grain radius in micrometers",
        default=40.0,
        type=float,
    )
    parser.add_argument(
        "--zr-rad-max",
        dest="zr_rad_max",
        help="Maximum zircon equivalent spherical grain radius in micrometers",
        default=100.0,
        type=float,
    )
    parser.add_argument(
        "--ap-thorium",
        dest="ap_thorium",
        help="Apatite thorium concentration in ppm",
        default=0.0,
        type=float,
    )
    parser.add_argument(
        "--zr-thorium",
        dest="zr_thorium",
        help="Zircon thorium concentration in ppm",
        default=0.0,
        type=float,
    )
    parser.add_argument(
        "--plot-type",
        dest="plot_type",
        help="eU versus radius plot type. 1 = apatite, 2 = zircon, 3 = both.",
        default=3,
        type=int,
    )
    parser.add_argument(
        "--save-plot",
        dest="save_plot",
        help="Save plot(s) to file",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--plot-file-format",
        dest="plot_file_format",
        help="File format for saving plot(s) to file (examples: png, pdf, svg, eps)",
        default="pdf",
        type=str,
    )
    parser.add_argument(
        "--plot-dpi",
        dest="plot_dpi",
        help="Saved plot resolution in dots per inch",
        default=300,
        type=int,
    )
    parser.add_argument(
        "--plot-style",
        dest="plot_style",
        help="Style sheet used for plotting. See https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html.",
        default="seaborn-colorblind",
        type=str,
    )
    parser.add_argument(
        "--plot-colormap",
        dest="plot_colormap",
        help="Colormap used for plotting. See https://matplotlib.org/stable/tutorials/colors/colormaps.html.",
        default="plasma",
        type=str,
    )
    parser.add_argument(
        "--plot-alpha",
        dest="plot_alpha",
        help="Transparency used for plotting fill colors",
        default=1.0,
        type=float,
    )
    parser.add_argument(
        "--plot-contour-lines",
        dest="plot_contour_lines",
        help="Number of contour lines used for plotting",
        default=12,
        type=int,
    )
    parser.add_argument(
        "--plot-contour-fills",
        dest="plot_contour_fills",
        help="Number of contour fill colors from the selected colormap",
        default=256,
        type=int,
    )
    parser.add_argument(
        "--no-display-plot",
        dest="no_display_plot",
        help="Do not display plots on the screen",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--tt-plot",
        dest="tt_plot",
        help="Create (and display) time-temperature history plot",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-v",
        "--verbose",
        help="Enable verbose output",
        action="store_true",
        default=False,
    )

    args = parser.parse_args()

    # Flip command-line flag to be opposite for function call
    # Function call expects display_plot = True for plot to be displayed
    display_plot = not args.no_display_plot

    eu_vs_radius(
        num_points=args.num_points,
        cooling_hist_type=args.cooling_hist_type,
        temp_max=args.temp_max,
        cooling_rate=args.cooling_rate,
        time_hist=args.time_hist,
        temp_hist=args.temp_hist,
        ap_u_min=args.ap_u_min,
        ap_u_max=args.ap_u_max,
        zr_u_min=args.zr_u_min,
        zr_u_max=args.zr_u_max,
        ap_rad_min=args.ap_rad_min,
        ap_rad_max=args.ap_rad_max,
        zr_rad_min=args.zr_rad_min,
        zr_rad_max=args.zr_rad_max,
        ap_thorium=args.ap_thorium,
        zr_thorium=args.zr_thorium,
        plot_type=args.plot_type,
        save_plot=args.save_plot,
        plot_file_format=args.plot_file_format,
        plot_dpi=args.plot_dpi,
        plot_style=args.plot_style,
        plot_colormap=args.plot_colormap,
        plot_alpha=args.plot_alpha,
        plot_contour_lines=args.plot_contour_lines,
        plot_contour_fills=args.plot_contour_fills,
        display_plot=display_plot,
        tt_plot=args.tt_plot,
        verbose=args.verbose,
        use_widget=False,
    )


if __name__ == "__main__":
    # execute only if run as a script
    main()
