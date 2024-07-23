#!/usr/bin/env python3

# Import libraries we need
import argparse
from tcplotter import rate_vs_age_tc


def main():
    """Command-line interface for the rate_vs_age_tc function in tcplotter."""
    parser = argparse.ArgumentParser(
        description="Calculates thermochronometer ages and closure temperatures for "
        "different cooling rates and eU concentrations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--num-points",
        dest="num_points",
        help="Number of points along x and y axes where ages/closure temperatures are "
        "calculated.",
        default=101,
        type=int,
    )
    parser.add_argument(
        "--cooling-rate-min",
        dest="cooling_rate_min",
        help="Minimum cooling rate in degrees C per Myr.",
        default=0.1,
        type=float,
    )
    parser.add_argument(
        "--cooling-rate-max",
        dest="cooling_rate_max",
        help="Maximum cooling rate in degrees C per Myr.",
        default=100.0,
        type=float,
    )
    parser.add_argument(
        "--temp-max",
        dest="temp_max",
        help="Max temperature for cooling history (in degrees C).",
        default=350.0,
        type=float,
    )
    parser.add_argument(
        "--ap-u1",
        dest="ap_u1",
        help="Apatite uranium concentration in ppm for upper plot panel",
        default=1.0,
        type=float,
    )
    parser.add_argument(
        "--ap-u2",
        dest="ap_u2",
        help="Apatite uranium concentration in ppm for middle plot panel",
        default=20.0,
        type=float,
    )
    parser.add_argument(
        "--ap-u3",
        dest="ap_u3",
        help="Apatite uranium concentration in ppm for lower plot panel",
        default=150.0,
        type=float,
    )
    parser.add_argument(
        "--zr-u1",
        dest="zr_u1",
        help="Zircon uranium concentration in ppm for upper plot panel",
        default=10.0,
        type=float,
    )
    parser.add_argument(
        "--zr-u2",
        dest="zr_u2",
        help="Zircon uranium concentration in ppm for middle plot panel",
        default=200.0,
        type=float,
    )
    parser.add_argument(
        "--zr-u3",
        dest="zr_u3",
        help="Zircon uranium concentration in ppm for lower plot panel",
        default=4000.0,
        type=float,
    )
    parser.add_argument(
        "--ap-rad",
        dest="ap_rad",
        help="Apatite equivalent spherical grain radius in micrometers",
        default=45.0,
        type=float,
    )
    parser.add_argument(
        "--zr-rad",
        dest="zr_rad",
        help="Zircon equivalent spherical grain radius in micrometers",
        default=60.0,
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
        "--ahe-uncertainty",
        dest="ahe_uncertainty",
        help="Apatite U-Th/He age uncertainty fraction. Enter 0.1 for 10 percent "
        "uncertainty.",
        default=0.1,
        type=float,
    )
    parser.add_argument(
        "--aft-uncertainty",
        dest="aft_uncertainty",
        help="Apatite fission-track age uncertainty fraction. Enter 0.2 for 20 "
        "percent uncertainty.",
        default=0.2,
        type=float,
    )
    parser.add_argument(
        "--zhe-uncertainty",
        dest="zhe_uncertainty",
        help="Zircon U-Th/He age uncertainty fraction. Enter 0.1 for 10 percent "
        "uncertainty.",
        default=0.1,
        type=float,
    )
    parser.add_argument(
        "--plot-type",
        dest="plot_type",
        help="1 = Cooling rate versus closure temperature. 2 = Cooling rate versus age. 3 = Cooling "
        "rate versus age and closure temperature",
        default=3,
        type=int,
    )
    parser.add_argument(
        "--plot-age-min",
        dest="plot_age_min",
        help="Minimum age value in Ma for plots of cooling rate versus age. Only applies to plot_type 2 and 3.",
        default=0.5,
        type=float,
    )
    parser.add_argument(
        "--plot-age-max",
        dest="plot_age_max",
        help="Maximum age value in Ma for plots of cooling rate versus age. Only applies to plot_type 2 and 3.",
        default=1800.0,
        type=float,
    )
    parser.add_argument(
        "--plot-tc-min",
        dest="plot_tc_min",
        help="Minimum closure temperature value in deg. C for plots of cooling rate versus closure temperature. Only applies to plot_type 1 and 3.",
        default=0.0,
        type=float,
    )
    parser.add_argument(
        "--plot-tc-max",
        dest="plot_tc_max",
        help="Maximum closure temperature value in deg. C for plots of cooling rate versus closure temperature. Only applies to plot_type 1 and 3.",
        default=200.0,
        type=float,
    )
    parser.add_argument(
        "--save-plot",
        dest="save_plot",
        help="Save plot to file",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--plot-file-format",
        dest="plot_file_format",
        help="File format for saving plot to file (examples: png, pdf, svg, eps)",
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
        default="seaborn-v0_8-colorblind",
        type=str,
    )
    parser.add_argument(
        "--plot-alpha",
        dest="plot_alpha",
        help="Transparency used for plotting fill colors for age swath plots",
        default=0.6,
        type=float,
    )
    parser.add_argument(
        "--no-plot-grid",
        dest="no_plot_grid",
        help="Do not display grid lines on plot",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--no-display-plot",
        dest="no_display_plot",
        help="Do not display plot on the screen",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--no-clean-up-files",
        dest="no_clean_up_files",
        help="Do not delete temporary files",
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

    # Flip command-line flags to be opposite for function call
    # Function call expects
    # - plot_grid = True for grid lines to be added to plot
    # - display_plot = True for plot to be displayed
    # - clean_up_files = True for temporary files to be deleted
    plot_grid = not args.no_plot_grid
    display_plot = not args.no_display_plot
    clean_up_files = not args.no_clean_up_files

    rate_vs_age_tc(
        num_points=args.num_points,
        cooling_rate_min=args.cooling_rate_min,
        cooling_rate_max=args.cooling_rate_max,
        temp_max=args.temp_max,
        ap_u1=args.ap_u1,
        ap_u2=args.ap_u2,
        ap_u3=args.ap_u3,
        zr_u1=args.zr_u1,
        zr_u2=args.zr_u2,
        zr_u3=args.zr_u3,
        ap_rad=args.ap_rad,
        zr_rad=args.zr_rad,
        ap_thorium=args.ap_thorium,
        zr_thorium=args.zr_thorium,
        ahe_uncertainty=args.ahe_uncertainty,
        aft_uncertainty=args.aft_uncertainty,
        zhe_uncertainty=args.zhe_uncertainty,
        plot_type=args.plot_type,
        plot_age_min=args.plot_age_min,
        plot_age_max=args.plot_age_max,
        plot_tc_min=args.plot_tc_min,
        plot_tc_max=args.plot_tc_max,
        save_plot=args.save_plot,
        plot_file_format=args.plot_file_format,
        plot_dpi=args.plot_dpi,
        plot_style=args.plot_style,
        plot_alpha=args.plot_alpha,
        plot_grid=plot_grid,
        display_plot=display_plot,
        clean_up_files=clean_up_files,
        verbose=args.verbose,
        use_widget=False,
    )


if __name__ == "__main__":
    # execute only if run as a script
    main()
