# Import libraries we need
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np
import os
from pathlib import Path
from scipy.interpolate import interp1d, RegularGridInterpolator
import shutil
import subprocess


# Define function for calculating effective uranium concentration
def calc_eu(uranium, thorium):
    """Calculates effective uranium concentration from U, Th inputs"""
    return uranium + 0.238 * thorium


# Define function to find which version of the RDAAM_He/ketch_aft to use
def get_tc_exec(command):
    """Returns the location of the RDAAM_He or ketch_aft executable"""
    if shutil.which(command) is not None:
        tc_exec = command
    elif Path("bin/" + command).is_file():
        tc_exec = "bin/" + command
    else:
        raise FileNotFoundError(
            f"Age calculation program {command} not found. See Troubleshooting in tcplotter docs online."
        )
    return tc_exec


# Function for reading age data file
def read_age_data(file):
    """Reads in age data from a csv file"""
    # Make empty lists for column values
    ahe_age = []
    ahe_uncertainty = []
    ahe_eu = []
    ahe_radius = []
    zhe_age = []
    zhe_uncertainty = []
    zhe_eu = []
    zhe_radius = []

    # Read in data file and create nested lists of values
    with open(file, "r") as file:
        data = file.read().splitlines()
        for i in range(1, len(data)):
            # Split lines by commas
            data[i] = data[i].split(",")
            # Strip whitespace
            data[i] = [line.strip() for line in data[i]]
            # Use values only if the eU and radius were provided
            if (len(data[i][3]) > 0) and (len(data[i][4]) > 0):
                # Append AHe data if the age type is AHe
                if data[i][0].lower() == "ahe":
                    ahe_age.append(float(data[i][1]))
                    ahe_uncertainty.append(float(data[i][2]))
                    ahe_eu.append(float(data[i][3]))
                    ahe_radius.append(float(data[i][4]))
                # Append ZHe data if the age type is ZHe
                elif data[i][0].lower() == "zhe":
                    zhe_age.append(float(data[i][1]))
                    zhe_uncertainty.append(float(data[i][2]))
                    zhe_eu.append(float(data[i][3]))
                    zhe_radius.append(float(data[i][4]))
        # Create new lists with data file values
        ahe_data = [ahe_age, ahe_uncertainty, ahe_eu, ahe_radius]
        zhe_data = [zhe_age, zhe_uncertainty, zhe_eu, zhe_radius]

    return ahe_data, zhe_data


def chi_squared(observed, expected, std):
    """Returns the reduced chi-squared value for input array data."""
    misfit = 0
    for i in range(len(observed)):
        misfit += (observed[i]-expected[i])**2 / std[i]**2
    # Scale goodness-of-fit by number of ages
    return misfit / len(observed)


def calculate_misfit(age_data, age_type, age_list, param_x, param_y):
    """Calculates misfit between measured and predicted ages."""
    predicted_ages = []
    measured_ages = []
    std_dev = []
    for i in range(len(age_data[0])):
        # Create interpolation function
        age_grid = np.array(age_list).reshape((len(param_x), len(param_y)))
        age_interp = RegularGridInterpolator((param_x, param_y), age_grid)
        # Append interpolated age if within the eU and radius ranges
        # Check if value is not within eU range
        if not (min(param_x) <= age_data[2][i] <= max(param_x)):
            print(f"Warning: eU concentration for {age_type} age {i + 1} not within modelled range.")
            print(f"         This age will be excluded from the misfit calculation.")
        elif not (min(param_y) <= age_data[3][i] <= max(param_y)):
            print(f"Warning: Grain radius for {age_type} age {i + 1} not within modelled range.")
            print(f"         This age will be excluded from the misfit calculation.")
        else:
            predicted_ages.append(age_interp([age_data[2][i], age_data[3][i]]))
            # Include only measured ages within the eU/radius ranges
            measured_ages.append(age_data[0][i])
            std_dev.append(age_data[1][i])
    # Calculate misfit
    misfit = chi_squared(measured_ages, predicted_ages, std_dev)
    n_ages = len(measured_ages)
    return misfit[0], n_ages


# Define function for creating plot of cooling rates
def time_vs_temp(
    cooling_rate_min=0.1,
    cooling_rate_slow=1.0,
    cooling_rate_avg=10.0,
    cooling_rate_max=100.0,
    temp_max=350.0,
    time_max=50.0,
    save_plot=False,
    plot_file_format="pdf",
    plot_dpi=300,
    plot_style="seaborn-whitegrid",
    fill_between=True,
    display_plot=True,
):
    """
    Plots cooling rate lines for different input rates.

    Parameters
    ----------
    cooling_rate_min : float or int, default=0.1
        Minimum cooling rate to plot in degrees C / Myr.
    cooling_rate_slow : float or int, default=1.0
        "Slow" cooling rate to plot in degrees C / Myr.
    cooling_rate_avg : float or int, default=10.0
        "Average" cooling rate to plot in degrees C / Myr.
    cooling_rate_max : float or int, default=100.0
        Maximum cooling rate to plot in degrees C / Myr.
    temp_max : float or int, default=350.0
        Maximum temperature for cooling history in degrees C.
    time_max : float or int, default=50.0
        Maximum value for time on x-axis of plot in millions of years ago (Ma).
    save_plot : bool, default=False
        Flag for whether to save the plot to a file.
    plot_file_format : str, default='pdf'
        File format for saving plot to file (examples: png, pdf, svg, eps).
    plot_dpi : int, default=300
        Saved plot resolution in dots per inch.
    plot_style : str, default='seaborn-whitegrid'
        Style sheet used for plotting. See https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html.
    fill_between : bool, default=True
        Flag for whether to fill area between min, max cooling rates.
    display_plot : bool, default=True
        Flag for whether to display the plot.

    Returns
    -------
    None
    """

    # Ensure relative paths work by setting working dir to dir containing this script file
    wd_orig = os.getcwd()
    script_path = os.path.abspath(__file__)
    dir_name = os.path.dirname(script_path)
    os.chdir(dir_name)

    # Find time and temperature bounds for plot
    time_plot_min = min(time_max, temp_max / cooling_rate_min)
    temp_plot_min = min(temp_max, cooling_rate_min * time_plot_min)
    time_plot_slow = min(time_max, temp_max / cooling_rate_slow)
    temp_plot_slow = min(temp_max, cooling_rate_slow * time_plot_slow)
    time_plot_avg = min(time_max, temp_max / cooling_rate_avg)
    temp_plot_avg = min(temp_max, cooling_rate_avg * time_plot_avg)
    time_plot_max = min(time_max, temp_max / cooling_rate_max)
    temp_plot_max = min(temp_max, cooling_rate_max * time_plot_max)

    # Create arrays of points to plot
    min_rate_x = np.array([time_plot_min, 0.0])
    min_rate_y = np.array([temp_plot_min, 0.0])
    slow_rate_x = np.array([time_plot_slow, 0.0])
    slow_rate_y = np.array([temp_plot_slow, 0.0])
    avg_rate_x = np.array([time_plot_avg, 0.0])
    avg_rate_y = np.array([temp_plot_avg, 0.0])
    max_rate_x = np.array([time_plot_max, 0.0])
    max_rate_y = np.array([temp_plot_max, 0.0])

    # Set plot style
    plt.style.use(plot_style)

    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))

    if fill_between:
        # Define fill ranges
        min_rate_filly = np.array([temp_max, 0.0])
        min_rate_fillx = np.array([temp_max / cooling_rate_min, 0.0])
        max_rate_fillx = np.array([temp_max / cooling_rate_max, 0.0])

        # Plot fill
        ax.fill_betweenx(
            min_rate_filly,
            min_rate_fillx,
            max_rate_fillx,
            color="black",
            alpha=0.15,
            label="Range of model cooling rates",
        )

    # Plot lines
    ax.plot(min_rate_x, min_rate_y, color="black")
    ax.plot(slow_rate_x, slow_rate_y, color="black")
    ax.plot(avg_rate_x, avg_rate_y, color="black")
    ax.plot(max_rate_x, max_rate_y, color="black")

    # Set axis tick label format
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_major_formatter(ScalarFormatter())

    # Set plot x and y range
    ax.set_xlim([0.0, time_max])
    ax.set_ylim([0.0, temp_max])

    # Add axis labels
    ax.set_xlabel("Time (Ma)")
    ax.set_ylabel("Temperature (°C)")

    # Flip axis directions
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()

    # Use tight layout
    plt.tight_layout()

    # Save plot if requested
    if save_plot:
        # Set plot filename and save plot
        plot_filename = "time_vs_temp_" + str(plot_dpi) + "dpi." + plot_file_format
        plt.savefig(wd_orig + "/" + plot_filename, dpi=plot_dpi)

    # Display plot if requested
    if display_plot:
        plt.show()

    # Revert to original working directory
    os.chdir(wd_orig)

    return None


# Define function for making contour plot of cooling ages and closure temperatures
def eu_vs_radius(
    num_points=21,
    cooling_hist_type=1,
    temp_max=350.0,
    cooling_rate=10.0,
    time_hist=[0.0, 10.0, 25.0, 35.0],
    temp_hist=[0.0, 75.0, 50.0, 350.0],
    ap_u_min=1.0,
    ap_u_max=150.0,
    zr_u_min=1.0,
    zr_u_max=4000.0,
    ap_rad_min=40.0,
    ap_rad_max=100.0,
    zr_rad_min=40.0,
    zr_rad_max=100.0,
    ap_thorium=0.0,
    zr_thorium=0.0,
    plot_type=3,
    save_plot=False,
    plot_file_format="pdf",
    plot_dpi=300,
    plot_style="seaborn-colorblind",
    plot_colormap="plasma",
    plot_alpha=1.0,
    plot_contour_lines=12,
    plot_contour_fills=256,
    age_data_file="",
    calc_misfit=False,
    display_plot=True,
    tt_plot=False,
    verbose=False,
    use_widget=False,
):
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
    temp_max : float, default=350.0
        Max temperature for cooling history (in degrees C). Option only for cooling history type 1.
    cooling_rate : float, default=10.0
        Cooling rate in degrees C per Myr. Option only for cooling history type 1.
    time_hist : list of floats or ints, default=[0.0, 10.0, 25.0, 35.0]
        Time points defining cooling history in Ma (millions of years ago).
        NOTE: Present-day point should be first in list.
        Option only for cooling history type 2.
    temp_hist : list of floats or ints, default=[0.0, 75.0, 50.0, 350.0]
        Temperature points defining cooling history in degrees C.
        NOTE: Present-day point should be first in list.
        Option only for cooling history type 2.
    ap_u_min : float, default=1.0
        Minimum apatite uranium concentration in ppm.
    ap_u_max : float, default=150.0
        Maximum apatite uranium concentration in ppm.
    zr_u_min : float, default=1.0
        Minimum zircon uranium concentration in ppm.
    zr_u_max : float, default=4000.0
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
    plot_dpi : int, default=300
        Saved plot resolution in dots per inch.
    plot_style : str, default='seaborn-colorblind'
        Style sheet used for plotting. See https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html.
    plot_colormap : str, default='plasma'
        Colormap used for plotting. See https://matplotlib.org/stable/tutorials/colors/colormaps.html.
    plot_alpha : float, default=1.0
        Transparency used for plotting fill colors.
    plot_contour_lines : int, default=12
        Number of contour lines used for plotting.
    plot_contour_fills : int, default=256
        Number of contour fill colors from the selected colormap.
    age_data_file : str, default=''
        Filename for file containing measured thermochronometer ages.
    calc_misfit : bool, default=False
        Flag for whether a misfit should be calculated for measured and predicted ages.
    display_plot : bool, default=True
        Flag for whether to display the plot.
    tt_plot : bool, default=False
        Flag for whether to create/display the time-temperature history plot.
    verbose : bool, default=False
        Enable/disable verbose output.
    use_widget : bool, default=False
        Enable/disable IPython progress bar widget. Disabled for command-line usage.
    """

    # Check to see whether ipywidgets and IPython are available for widget use
    # If not, disable widgets and display a warning
    if use_widget:
        try:
            import ipywidgets as widgets
        except ModuleNotFoundError:
            print("Warning: ipywidgets module not found. Disabling graphical progress bar.")
            use_widget = False
    if use_widget:
        try:
            from IPython.display import display
        except ModuleNotFoundError:
            print(
                "Warning: IPython.display module not found. Disabling graphical progress bar."
            )
            use_widget = False

    # Read in measured ages from file, if age_data_file is defined
    if len(age_data_file) > 0:
        ahe_age_data, zhe_age_data = read_age_data(age_data_file)
    else:
        ahe_age_data = []
        zhe_age_data = []

    # Ensure relative paths work by setting working dir to dir containing this script file
    wd_orig = os.getcwd()
    script_path = os.path.abspath(__file__)
    dir_name = os.path.dirname(script_path)
    os.chdir(dir_name)

    # Define cooling history using constant cooling rate
    if cooling_hist_type == 1:
        # Define time and temperature histories
        start_time = temp_max / cooling_rate
        time_hist = [0.0, start_time]
        temp_hist = [0.0, temp_max]

    # Option 2: Define time-temperature history using list of tT points
    elif cooling_hist_type == 2:
        pass

    # Raise error if an unsupported value is given for cooling_hist_type
    else:
        raise ValueError("Bad value for cooling_hist_type. Should be 1 or 2.")

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
        model_type = "apatite age/Tc (eU vs. radius)"
    elif plot_type == 2:
        model_type = "zircon age/Tc (eU vs. radius)"
    elif plot_type == 3:
        model_type = "apatite/zircon age/Tc (eU vs. radius)"
    else:
        raise ValueError("Bad value for plot_type. Should be 1, 2, or 3.")

    # Define time-temperature history filename
    tt_file = "simple_time_temp.txt"

    # Get age calculation executable(s) to use
    rdaam_command = get_tc_exec("RDAAM_He")

    # Set plot style
    plt.style.use(plot_style)

    # Define plot size and number of subplots
    fig_width = 10
    if plot_type < 3:
        fig_height = 5
        # Make plot longer if plotting data
        if plot_type == 1:
            if len(ahe_age_data) > 0: fig_height += 1.5
        else:
            if len(zhe_age_data) > 0: fig_height += 1.5
        # Create figure and axes
        fig, ax = plt.subplots(1, 2, figsize=(fig_width, fig_height))
    else:
        fig_height = 10
        # Make plot longer if plotting data
        if (len(ahe_age_data) > 0): fig_height += 1.5
        if (len(zhe_age_data) > 0): fig_height += 1.5
        # Create figure and axes
        fig, ax = plt.subplots(2, 2, figsize=(fig_width, fig_height))

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
    with open(tt_file, "w") as f:
        for i in range(len(time_hist)):
            f.write(f"{time_hist[i]:.4f},{temp_hist[i]:.1f}\n")

    # Echo total model run time and cooling rate
    if verbose and cooling_hist_type == 1:
        print(
            f"Cooling from {temp_max:.1f}°C at a rate of {cooling_rate:.1f} °C/Myr will require {start_time:.2f} million years"
        )

    # Create visual progress bar, if enabled
    if use_widget and not verbose:
        s = widgets.IntProgress(
            value=0,
            min=0,
            max=total_models,
            description="Calculating:",
            bar_style="",  # 'success', 'info', 'warning', 'danger' or ''
            style={"bar_color": "#ff6666"},
            orientation="horizontal",
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
                        f"Calculating {model_type} - {int(round(100 * model_count / total_models)):3d}% ({model_count:5d} / {total_models:5d})\r",
                        end="",
                    )

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
            command = (
                rdaam_command
                + " "
                + tt_file
                + " "
                + str(ap_radius)
                + " "
                + str(ap_uranium)
                + " "
                + str(ap_thorium)
                + " "
                + str(zr_radius)
                + " "
                + str(zr_uranium)
                + " "
                + str(zr_thorium)
            )
            p = subprocess.Popen(
                command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
            )

            # Parse output for ages
            stdout = p.stdout.readlines()
            corr_ahe_age = stdout[0].split()[7].decode("UTF-8")
            corr_zhe_age = stdout[1].split()[7].decode("UTF-8")

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
                    f"AHe: {float(corr_ahe_age):.2f} Ma (Tc: {ahe_tc:.1f}°C); ZHe: {float(corr_zhe_age):.2f} Ma (Tc: {zhe_tc:.1f}°C)"
                )

    # Clean up Tt file
    os.remove(tt_file)

    # Calculate age misfits if age data file is used and option is selected
    if calc_misfit and len(age_data_file) > 0:
        total_misfit = 0
        if len(ahe_age_data) > 0:
            ahe_misfit, ahe_n_ages = calculate_misfit(ahe_age_data, "AHe", ahe_age_list, ap_x, ap_y)
            print(f"AHe misfit (n = {ahe_n_ages} ages): {ahe_misfit}")
        if len(zhe_age_data) > 0:
            zhe_misfit, zhe_n_ages = calculate_misfit(zhe_age_data, "ZHe", zhe_age_list, zr_x, zr_y)
            print(f"ZHe misfit (n = {zhe_n_ages} ages): {zhe_misfit}")
        total_misfit += ahe_misfit + zhe_misfit
        print(f"Total misfit (n = {ahe_n_ages + zhe_n_ages} ages): {total_misfit}")

    # Apatite eU versus radius
    if plot_type == 1:
        # Create age contour lines
        ap_contours_age = ax[0].tricontour(
            ap_x_list,
            ap_y_list,
            ahe_age_list,
            plot_contour_lines,
            linewidths=0.5,
            colors="black",
        )
        # Add age contour labels
        ax[0].clabel(ap_contours_age)
        # Determine bounds for contour colors if plotting age data
        if len(ahe_age_data) > 0:
            age_min = min(min(ahe_age_list), min(ahe_age_data[0]))
            age_max = max(max(ahe_age_list), max(ahe_age_data[0]))
        else:
            age_min = min(ahe_age_list)
            age_max = max(ahe_age_list)
        # Create age contour fill
        ap_contourf_age = ax[0].tricontourf(
            ap_x_list,
            ap_y_list,
            ahe_age_list,
            plot_contour_fills,
            cmap=plot_colormap,
            vmin=age_min,
            vmax=age_max,
            alpha=plot_alpha,
        )
        if len(ahe_age_data) > 0:
            ahe_label = f"AHe data (n = {ahe_n_ages}"
            if calc_misfit:
                ahe_label += f"; misfit: {ahe_misfit:.2f}"
            ahe_label += ")"
            age_data_plot = ax[0].scatter(
                x=ahe_age_data[2],
                y=ahe_age_data[3],
                c=ahe_age_data[0],
                edgecolors="black",
                cmap=plot_colormap,
                vmin=age_min,
                vmax=age_max,
                label=ahe_label
            )
            ax[0].set_xlim([min(ap_x_list), max(ap_x_list)])
            ax[0].set_ylim([min(ap_y_list), max(ap_y_list)])
            ax[0].legend()

            # Plot the colorbar if plotting data
            fig.colorbar(age_data_plot, ax=ax[0], orientation="horizontal", label="Apatite (U-Th)/He age (Ma)")

        # This is the fix for the white lines between contour levels
        for c in ap_contourf_age.collections:
            c.set_edgecolor("face")

        # Create closure temperature contour lines
        ap_contours_tc = ax[1].tricontour(
            ap_x_list,
            ap_y_list,
            ahe_tc_list,
            plot_contour_lines,
            linewidths=0.5,
            colors="black",
        )
        # Add closure temperature contour labels
        ax[1].clabel(ap_contours_tc)
        # Create closure temperature contour fill
        ap_contourf_tc = ax[1].tricontourf(
            ap_x_list,
            ap_y_list,
            ahe_tc_list,
            plot_contour_fills,
            cmap=plot_colormap,
            alpha=plot_alpha,
        )

        if len(ahe_age_data) > 0:
            # Plot the colorbar if plotting data
            fig.colorbar(ap_contourf_tc, ax=ax[1], orientation="horizontal", label="Apatite (U-Th)/He closure temperature (°C)")

        # This is the fix for the white lines between contour levels
        for c in ap_contourf_tc.collections:
            c.set_edgecolor("face")

    # Zircon eU versus radius
    elif plot_type == 2:
        # Create age contour lines
        zr_contours_age = ax[0].tricontour(
            zr_x_list,
            zr_y_list,
            zhe_age_list,
            plot_contour_lines,
            linewidths=0.5,
            colors="k",
        )
        # Add age contour labels
        ax[0].clabel(zr_contours_age)
        # Determine bounds for contour colors if plotting age data
        if len(zhe_age_data) > 0:
            age_min = min(min(zhe_age_list), min(zhe_age_data[0]))
            age_max = max(max(zhe_age_list), max(zhe_age_data[0]))
        else:
            age_min = min(zhe_age_list)
            age_max = max(zhe_age_list)
        # Create age contour fill
        zr_contourf_age = ax[0].tricontourf(
            zr_x_list,
            zr_y_list,
            zhe_age_list,
            plot_contour_fills,
            cmap=plot_colormap,
            vmin=age_min,
            vmax=age_max,
            alpha=plot_alpha,
        )
        if len(zhe_age_data) > 0:
            zhe_label = f"ZHe data (n = {zhe_n_ages}"
            if calc_misfit:
                zhe_label += f"; misfit: {zhe_misfit:.2f}"
            zhe_label += ")"
            age_data_plot = ax[0].scatter(
                x=zhe_age_data[2],
                y=zhe_age_data[3],
                c=zhe_age_data[0],
                edgecolors="black",
                cmap=plot_colormap,
                vmin=age_min,
                vmax=age_max,
                label=zhe_label
            )
            ax[0].set_xlim([min(zr_x_list), max(zr_x_list)])
            ax[0].set_ylim([min(zr_y_list), max(zr_y_list)])
            ax[0].legend()

            # Plot the colorbar if plotting data
            fig.colorbar(age_data_plot, ax=ax[0], orientation="horizontal", label="Zircon (U-Th)/He age (Ma)")

        # This is the fix for the white lines between contour levels
        for c in zr_contourf_age.collections:
            c.set_edgecolor("face")

        # Create closure temperature contour lines
        zr_contours_tc = ax[1].tricontour(
            zr_x_list,
            zr_y_list,
            zhe_tc_list,
            plot_contour_lines,
            linewidths=0.5,
            colors="k",
        )
        # Add closure temperature contour labels
        ax[1].clabel(zr_contours_tc)
        # Create closure temperature contour fill
        zr_contourf_tc = ax[1].tricontourf(
            zr_x_list,
            zr_y_list,
            zhe_tc_list,
            plot_contour_fills,
            cmap=plot_colormap,
            alpha=plot_alpha,
        )

        if len(zhe_age_data) > 0:
            # Plot the colorbar if plotting data
            fig.colorbar(zr_contourf_tc, ax=ax[1], orientation="horizontal", label="Zircon (U-Th)/He closure temperature (°C)")

        # This is the fix for the white lines between contour levels
        for c in zr_contourf_tc.collections:
            c.set_edgecolor("face")

    # Apatite and zircon eU versus radius
    else:
        # Create age contour lines
        ap_contours_age = ax[0][0].tricontour(
            ap_x_list,
            ap_y_list,
            ahe_age_list,
            plot_contour_lines,
            linewidths=0.5,
            colors="k",
        )
        # Add age contour labels
        ax[0][0].clabel(ap_contours_age)
        # Determine bounds for contour colors if plotting age data
        if len(ahe_age_data) > 0:
            age_min = min(min(ahe_age_list), min(ahe_age_data[0]))
            age_max = max(max(ahe_age_list), max(ahe_age_data[0]))
        else:
            age_min = min(ahe_age_list)
            age_max = max(ahe_age_list)
        # Create age contour fill
        ap_contourf_age = ax[0][0].tricontourf(
            ap_x_list,
            ap_y_list,
            ahe_age_list,
            plot_contour_fills,
            cmap=plot_colormap,
            vmin=age_min,
            vmax=age_max,
            alpha=plot_alpha,
        )
        if len(ahe_age_data) > 0:
            ahe_label = f"AHe data (n = {ahe_n_ages}"
            if calc_misfit:
                ahe_label += f"; misfit: {ahe_misfit:.2f}"
            ahe_label += ")"
            age_data_plot = ax[0][0].scatter(
                x=ahe_age_data[2],
                y=ahe_age_data[3],
                c=ahe_age_data[0],
                edgecolors="black",
                cmap=plot_colormap,
                vmin=age_min,
                vmax=age_max,
                label=ahe_label
            )
            ax[0][0].set_xlim([min(ap_x_list), max(ap_x_list)])
            ax[0][0].set_ylim([min(ap_y_list), max(ap_y_list)])
            ax[0][0].legend()

            # Plot the colorbar if plotting data
            fig.colorbar(age_data_plot, ax=ax[0][0], orientation="horizontal", label="Apatite (U-Th)/He age (Ma)")

        # This is the fix for the white lines between contour levels
        for c in ap_contourf_age.collections:
            c.set_edgecolor("face")

        # Create closure temperature contour lines
        ap_contours_tc = ax[0][1].tricontour(
            ap_x_list,
            ap_y_list,
            ahe_tc_list,
            plot_contour_lines,
            linewidths=0.5,
            colors="k",
        )
        # Add closure temperature contour labels
        ax[0][1].clabel(ap_contours_tc)
        # Create closure temperature contour fill
        ap_contourf_tc = ax[0][1].tricontourf(
            ap_x_list,
            ap_y_list,
            ahe_tc_list,
            plot_contour_fills,
            cmap=plot_colormap,
            alpha=plot_alpha,
        )

        if len(ahe_age_data) > 0:
            # Plot the colorbar if plotting data
            fig.colorbar(ap_contourf_tc, ax=ax[0][1], orientation="horizontal", label="Apatite (U-Th)/He closure temperature (°C)")

        # This is the fix for the white lines between contour levels
        for c in ap_contourf_tc.collections:
            c.set_edgecolor("face")

        # Create age contour lines
        zr_contours_age = ax[1][0].tricontour(
            zr_x_list,
            zr_y_list,
            zhe_age_list,
            plot_contour_lines,
            linewidths=0.5,
            colors="k",
        )
        # Add age contour labels
        ax[1][0].clabel(zr_contours_age)
        # Determine bounds for contour colors if plotting age data
        if len(zhe_age_data) > 0:
            age_min = min(min(zhe_age_list), min(zhe_age_data[0]))
            age_max = max(max(zhe_age_list), max(zhe_age_data[0]))
        else:
            age_min = min(zhe_age_list)
            age_max = max(zhe_age_list)
        # Create age contour fill
        zr_contourf_age = ax[1][0].tricontourf(
            zr_x_list,
            zr_y_list,
            zhe_age_list,
            plot_contour_fills,
            cmap=plot_colormap,
            vmin=age_min,
            vmax=age_max,
            alpha=plot_alpha,
        )
        if len(zhe_age_data) > 0:
            zhe_label = f"ZHe data  (n = {zhe_n_ages}"
            if calc_misfit:
                zhe_label += f"; misfit: {zhe_misfit:.2f}"
            zhe_label += ")"
            age_data_plot = ax[1][0].scatter(
                x=zhe_age_data[2],
                y=zhe_age_data[3],
                c=zhe_age_data[0],
                edgecolors="black",
                cmap=plot_colormap,
                vmin=age_min,
                vmax=age_max,
                label=zhe_label
            )
            ax[1][0].set_xlim([min(zr_x_list), max(zr_x_list)])
            ax[1][0].set_ylim([min(zr_y_list), max(zr_y_list)])
            ax[1][0].legend()

            # Plot the colorbar if plotting data
            fig.colorbar(age_data_plot, ax=ax[1][0], orientation="horizontal", label="Zircon (U-Th)/He age (Ma)")

        # This is the fix for the white lines between contour levels
        for c in zr_contourf_age.collections:
            c.set_edgecolor("face")

        # Create closure temperature contour lines
        zr_contours_tc = ax[1][1].tricontour(
            zr_x_list,
            zr_y_list,
            zhe_tc_list,
            plot_contour_lines,
            linewidths=0.5,
            colors="k",
        )
        # Add closure temperature contour labels
        ax[1][1].clabel(zr_contours_tc)
        # Create closure temperature contour fill
        zr_contourf_tc = ax[1][1].tricontourf(
            zr_x_list,
            zr_y_list,
            zhe_tc_list,
            plot_contour_fills,
            cmap=plot_colormap,
            alpha=plot_alpha,
        )

        if len(zhe_age_data) > 0:
            # Plot the colorbar if plotting data
            fig.colorbar(zr_contourf_tc, ax=ax[1][1], orientation="horizontal", label="Zircon (U-Th)/He closure temperature (°C)")

        # This is the fix for the white lines between contour levels
        for c in zr_contourf_tc.collections:
            c.set_edgecolor("face")

    # Format plot labels

    # Apatite eU versus radius
    if plot_type == 1:
        ax[0].set_title("Apatite (U-Th)/He age [Ma]")
        ax[1].set_title("Apatite (U-Th)/He closure temperature [°C]")

    # Zircon eU versus radius
    elif plot_type == 2:
        ax[0].set_title("Zircon (U-Th)/He age [Ma]")
        ax[1].set_title("Zircon (U-Th)/He closure temperature [°C]")

    # Apatite and zircon eU versus radius
    else:
        ax[0][0].set_title("Apatite (U-Th)/He age [Ma]")
        ax[0][1].set_title("Apatite (U-Th)/He closure temperature [°C]")
        ax[1][0].set_title("Zircon (U-Th)/He age [Ma]")
        ax[1][1].set_title("Zircon (U-Th)/He closure temperature [°C]")

    # Apatite or Zircon eU versus radius
    if plot_type < 3:
        ax[0].set_xlabel("Effective uranium (ppm)")
        ax[1].set_xlabel("Effective uranium (ppm)")
        ax[0].set_ylabel("Equivalent spherical radius (µm)")
        ax[1].set_ylabel("Equivalent spherical radius (µm)")

    # Apatite and zircon eU versus radius
    else:
        ax[0][0].set_xlabel("Effective uranium (ppm)")
        ax[0][1].set_xlabel("Effective uranium (ppm)")
        ax[0][0].set_ylabel("Equivalent spherical radius (µm)")
        ax[0][1].set_ylabel("Equivalent spherical radius (µm)")
        ax[1][0].set_xlabel("Effective uranium (ppm)")
        ax[1][1].set_xlabel("Effective uranium (ppm)")
        ax[1][0].set_ylabel("Equivalent spherical radius (µm)")
        ax[1][1].set_ylabel("Equivalent spherical radius (µm)")

    # Use tight layout for subplots
    plt.tight_layout()

    # Save plot if desired
    if save_plot:
        # Set file name prefix
        plot_filename = "eu_vs_radius"

        # Define plot filename based on type of plot and save plot
        if plot_type == 1:
            plot_savename = (
                plot_filename + "_apatite_" + str(plot_dpi) + "dpi." + plot_file_format
            )
        elif plot_type == 2:
            plot_savename = (
                plot_filename + "_zircon_" + str(plot_dpi) + "dpi." + plot_file_format
            )
        else:
            plot_savename = (
                plot_filename
                + "_apatite_zircon_"
                + str(plot_dpi)
                + "dpi."
                + plot_file_format
            )
        plt.savefig(wd_orig + "/" + plot_savename, dpi=plot_dpi)

    # Display plot if desired
    if display_plot:
        plt.show()

    # Create tT history plot if requested
    if tt_plot:
        # Create figure 2
        fig2, ax2 = plt.subplots(1, 1, figsize=(6, 5))

        # Plot tT history
        ax2.plot(time_hist, temp_hist, color="black")

        # Set plot x and y range
        ax2.set_xlim([0.0, max(time_hist)])
        ax2.set_ylim([0.0, max(temp_hist)])

        # Add axis labels
        ax2.set_xlabel("Time (Ma)")
        ax2.set_ylabel("Temperature (°C)")

        # Add title
        ax2.set_title("Time-temperature history")

        # Flip axis directions
        plt.gca().invert_xaxis()
        plt.gca().invert_yaxis()

        # Use tight layout
        plt.tight_layout()

        # Save plot if desired
        if save_plot:
            # Define plot filename and save plot
            plot_savename2 = (
                plot_filename
                + "_tT_history_"
                + str(plot_dpi)
                + "dpi."
                + plot_file_format
            )
            plt.savefig(wd_orig + "/" + plot_savename2, dpi=plot_dpi)

        # Display plot if desired
        if display_plot:
            plt.show()

    # Revert to original working directory
    os.chdir(wd_orig)

    return None


# Define function for creating plot of cooling rates
def rate_vs_radius_eu(
    num_points=21,
    cooling_rate_min=0.1,
    cooling_rate_max=100.0,
    temp_max=350.0,
    ap_u_min=1.0,
    ap_u_max=150.0,
    ap_u_ref=10.0,
    zr_u_min=1.0,
    zr_u_max=4000.0,
    zr_u_ref=100.0,
    ap_rad_min=40.0,
    ap_rad_max=100.0,
    ap_rad_ref=45.0,
    zr_rad_min=40.0,
    zr_rad_max=100.0,
    zr_rad_ref=60.0,
    ap_thorium=0.0,
    zr_thorium=0.0,
    plot_type=3,
    save_plot=False,
    plot_file_format="pdf",
    plot_dpi=300,
    plot_style="seaborn-colorblind",
    plot_colormap="plasma",
    plot_alpha=1.0,
    plot_contour_lines=12,
    plot_contour_fills=256,
    display_plot=True,
    verbose=False,
    use_widget=False,
):
    """
    Calculates thermochronometer ages and closure temperatures for different cooling rates, effective uranium
    concentrations, and equivalent spherical radii.

    Parameters
    ----------
    num_points : int, default=21
        Number of points along x and y axes where ages/closure temperatures are
        calculated.
        NOTE: A value of num_points = 101 was used in the manuscript. It has been
        reduced here to make the plotting faster. Set this to 101 to reproduce
        the manuscript Figure 4.
    cooling_rate_min : float, default=0.1
        Minimum cooling rate in degrees C per Myr.
    cooling_rate_max : float, default=100.0
        Maximum cooling rate in degrees C per Myr.
    temp_max : float, default=350.0
        Max temperature for cooling history (in degrees C).
    ap_u_min : float, default=1.0
        Minimum apatite uranium concentration in ppm.
    ap_u_max : float, default=150.0
        Maximum apatite uranium concentration in ppm.
    ap_u_ref : float, default=10.0
        Apatite uranium concentration in ppm for rate versus radius plot.
    zr_u_min : float, default=1.0
        Minimum zircon uranium concentration in ppm.
    zr_u_max : float, default=4000.0
        Maximum zircon uranium concentration in ppm.
    zr_u_ref : float, default=100.0
        Zircon uranium concentration in ppm for rate versus radius plot.
    ap_rad_min : float, default=40.0
        Minimum apatite equivalent spherical grain radius in micrometers.
    ap_rad_max : float, default=100.0
        Maximum apatite equivalent spherical grain radius in micrometers.
    ap_rad_ref : float, default=45.0
        Apatite equivalent spherical grain radius in micrometers for rate versus eU plot.
    zr_rad_min : float, default=40.0
        Minimum zircon equivalent spherical grain radius in micrometers.
    zr_rad_max : float, default=100.0
        Maximum zircon equivalent spherical grain radius in micrometers.
    zr_rad_ref : float, default=60.0
        Zircon equivalent spherical grain radius in micrometers for rate versus eU plot.
    ap_thorium : float, default=0.0
        Apatite thorium concentration in ppm.
    zr_thorium : float, default=0.0
        Zircon thorium concentration in ppm.
    plot_type : int, default=3
        Cooling rate versus eU/radius.
        1 = apatite, 2 = zircon, 3 = both
    save_plot : bool, default=False
        Flag for whether to save the plot to a file.
    plot_file_format : str, default='pdf'
        File format for saving plot to file (examples: png, pdf, svg, eps).
    plot_dpi : int, default=300
        Saved plot resolution in dots per inch.
    plot_style : str, default='seaborn-colorblind'
        Style sheet used for plotting. See https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html.
    plot_colormap : str, default='plasma'
        Colormap used for plotting. See https://matplotlib.org/stable/tutorials/colors/colormaps.html.
    plot_alpha : float, default=1.0
        Transparency used for plotting fill colors.
    plot_contour_lines : int, default=12
        Number of contour lines used for plotting.
    plot_contour_fills : int, default=256
        Number of contour fill colors from the selected colormap.
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

    # Check to see whether ipywidgets and IPython are available for widget use
    # If not, disable widgets and display a warning
    if use_widget:
        try:
            import ipywidgets as widgets
        except ModuleNotFoundError:
            print("Warning: ipywidgets module not found. Disabling graphical progress bar.")
            use_widget = False
    if use_widget:
        try:
            from IPython.display import display
        except ModuleNotFoundError:
            print(
                "Warning: IPython.display module not found. Disabling graphical progress bar."
            )
            use_widget = False

    # Ensure relative paths work by setting working dir to dir containing this script file
    wd_orig = os.getcwd()
    script_path = os.path.abspath(__file__)
    dir_name = os.path.dirname(script_path)
    os.chdir(dir_name)

    # Create arrays of U concentrations
    ap_u = np.linspace(ap_u_min, ap_u_max, num_points)
    zr_u = np.linspace(zr_u_min, zr_u_max, num_points)

    # Create grain radius arrays
    ap_rad = np.linspace(ap_rad_min, ap_rad_max, num_points)
    zr_rad = np.linspace(zr_rad_min, zr_rad_max, num_points)

    # Create cooling rate array
    rates = np.logspace(
        start=np.log10(cooling_rate_min),
        stop=np.log10(cooling_rate_max),
        num=num_points,
    )

    # Calculate effective uranium
    ap_eu = calc_eu(ap_u, ap_thorium)
    zr_eu = calc_eu(zr_u, zr_thorium)

    # Total number of models
    total_models = len(ap_u) * len(rates) + len(ap_rad) * len(rates)

    # Screen output info
    if plot_type == 1:
        model_type = "apatite age/Tc (cooling rate vs. radius/eU)"
    elif plot_type == 2:
        model_type = "zircon age/Tc (cooling rate vs. radius/eU)"
    elif plot_type == 3:
        model_type = "apatite/zircon age/Tc (cooling rate vs. radius/eU)"
    else:
        raise ValueError("Bad value for parameter plot_type. Must be 1, 2, or 3.")

    # Define time-temperature history filename
    tt_file = "simple_time_temp.txt"

    # Get age calculation executable(s) to use
    rdaam_command = get_tc_exec("RDAAM_He")

    # Set plot style
    plt.style.use(plot_style)

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
    if use_widget and not verbose:
        s = widgets.IntProgress(
            value=0,
            min=0,
            max=total_models,
            description="Calculating:",
            bar_style="",  # 'success', 'info', 'warning', 'danger' or ''
            style={"bar_color": "#ff6666"},
            orientation="horizontal",
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
                        f"Calculating {model_type} - {int(round(100 * model_count / total_models)):3d}% ({model_count:5d} / {total_models:5d})\r",
                        end="",
                    )

            # Define parameters for this iteration
            rate = rates[i]
            ap_radius = ap_rad[j]
            zr_radius = zr_rad[j]
            ap_uranium = ap_u_ref
            zr_uranium = zr_u_ref
            ap_x_list1.append(rate)
            zr_x_list1.append(rate)
            ap_y_list1.append(ap_radius)
            zr_y_list1.append(zr_radius)

            # Write synthetic cooling history points to file
            start_time = temp_max / rate
            with open(tt_file, "w") as f:
                f.write("0.0,0.0\n")
                f.write("{0:.4f},{1:.1f}".format(start_time, temp_max))

            # Screen output
            if verbose:
                print(
                    f"Cooling from {temp_max:.1f}°C at a rate of {rate:.1f} °C/Myr will require {start_time:.2f} million years"
                )

            # Calculate (U-Th)/He ages
            command = (
                rdaam_command
                + " "
                + tt_file
                + " "
                + str(ap_radius)
                + " "
                + str(ap_uranium)
                + " "
                + str(ap_thorium)
                + " "
                + str(zr_radius)
                + " "
                + str(zr_uranium)
                + " "
                + str(zr_thorium)
            )
            p = subprocess.Popen(
                command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
            )

            # Parse output for ages
            stdout = p.stdout.readlines()
            corr_ahe_age = stdout[0].split()[7].decode("UTF-8")
            corr_zhe_age = stdout[1].split()[7].decode("UTF-8")

            # Find closure temperatures from cooling ages and thermal history
            tc_interp = interp1d([0.0, start_time], [0.0, temp_max])
            ahe_tc = tc_interp(float(corr_ahe_age))
            zhe_tc = tc_interp(float(corr_zhe_age))

            # Add closure temperatures to lists
            ahe_tc_list1.append(ahe_tc)
            zhe_tc_list1.append(zhe_tc)

            if verbose:
                print(
                    f"AHe: {float(corr_ahe_age):.2f} Ma (Tc: {ahe_tc:.1f}°C); ZHe: {float(corr_zhe_age):.2f} Ma (Tc: {zhe_tc:.1f}°C)"
                )

    # Loop over plotables - loop 2: rate versus eU
    for i in range(len(ap_x2)):
        for j in range(len(ap_y2)):
            model_count += 1
            if not verbose:
                if use_widget:
                    s.value = model_count
                else:
                    print(
                        f"Calculating {model_type} - {int(round(100 * (model_count) / total_models)):3d}% ({model_count:5d} / {total_models:5d})\r",
                        end="",
                    )

            # Define parameters for this iteration
            rate = rates[i]
            ap_radius = ap_rad_ref
            zr_radius = zr_rad_ref
            ap_uranium = ap_u[j]
            zr_uranium = zr_u[j]
            ap_x_list2.append(rate)
            zr_x_list2.append(rate)
            ap_y_list2.append(ap_uranium)
            zr_y_list2.append(zr_uranium)

            # Write synthetic cooling history points to file
            start_time = temp_max / rate
            with open(tt_file, "w") as f:
                f.write("0.0,0.0\n")
                f.write("{0:.4f},{1:.1f}".format(start_time, temp_max))

            # Screen output
            if verbose:
                print(
                    f"Cooling from {temp_max:.1f}°C at a rate of {rate:.1f} °C/Myr will require {start_time:.2f} million years"
                )

            # Calculate (U-Th)/He ages
            command = (
                rdaam_command
                + " "
                + tt_file
                + " "
                + str(ap_radius)
                + " "
                + str(ap_uranium)
                + " "
                + str(ap_thorium)
                + " "
                + str(zr_radius)
                + " "
                + str(zr_uranium)
                + " "
                + str(zr_thorium)
            )
            p = subprocess.Popen(
                command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
            )

            # Parse output for ages
            stdout = p.stdout.readlines()
            corr_ahe_age = stdout[0].split()[7].decode("UTF-8")
            corr_zhe_age = stdout[1].split()[7].decode("UTF-8")

            # Find closure temperatures from cooling ages and thermal history
            tc_interp = interp1d([0.0, start_time], [0.0, temp_max])
            ahe_tc = tc_interp(float(corr_ahe_age))
            zhe_tc = tc_interp(float(corr_zhe_age))

            # Add closure temperatures to lists
            ahe_tc_list2.append(ahe_tc)
            zhe_tc_list2.append(zhe_tc)

            if verbose:
                print(
                    f"AHe: {float(corr_ahe_age):.2f} Ma (Tc: {ahe_tc:.1f}°C); ZHe: {float(corr_zhe_age):.2f} Ma (Tc: {zhe_tc:.1f}°C)"
                )

    # Clean up temporary tt file
    os.remove(tt_file)

    # Plot only values for apatite (U-Th)/He
    if plot_type == 1:
        # --- Apatite cooling rate versus radius ---
        # Create closure temperature contour lines
        ap_contours_tc = ax[0].tricontour(
            ap_x_list1,
            ap_y_list1,
            ahe_tc_list1,
            plot_contour_lines,
            linewidths=0.5,
            colors="k",
        )
        # Use log x-axis scaling
        ax[0].set_xscale("log")
        # Add closure temperature contour labels
        ax[0].clabel(ap_contours_tc)
        # Create closure temperature contour fill
        ap_contourf_tc1 = ax[0].tricontourf(
            ap_x_list1,
            ap_y_list1,
            ahe_tc_list1,
            plot_contour_fills,
            cmap=plot_colormap,
            alpha=plot_alpha,
        )

        # This is the fix for the white lines between contour levels
        for c in ap_contourf_tc1.collections:
            c.set_edgecolor("face")

        # --- Apatite cooling rate versus eU plot ---
        # Create closure temperature contour lines
        ap_contours_tc = ax[1].tricontour(
            ap_x_list2,
            ap_y_list2,
            ahe_tc_list2,
            plot_contour_lines,
            linewidths=0.5,
            colors="k",
        )
        # Use log x-axis scaling
        ax[1].set_xscale("log")
        # Add closure temperature contour labels
        ax[1].clabel(ap_contours_tc)
        # Create closure temperature contour fill
        ap_contourf_tc2 = ax[1].tricontourf(
            ap_x_list2,
            ap_y_list2,
            ahe_tc_list2,
            plot_contour_fills,
            cmap=plot_colormap,
            alpha=plot_alpha,
        )

        # This is the fix for the white lines between contour levels
        for c in ap_contourf_tc2.collections:
            c.set_edgecolor("face")

    # Plot only values for zircon (U-Th)/He
    elif plot_type == 2:
        # --- Zircon cooling rate versus radius ---
        # Create closure temperature contour lines
        zr_contours_tc = ax[0].tricontour(
            zr_x_list1,
            zr_y_list1,
            zhe_tc_list1,
            plot_contour_lines,
            linewidths=0.5,
            colors="k",
        )
        # Use log x-axis scaling
        ax[0].set_xscale("log")
        # Add closure temperature contour labels
        ax[0].clabel(zr_contours_tc)
        # Create closure temperature contour fill
        zr_contourf_tc1 = ax[0].tricontourf(
            zr_x_list1,
            zr_y_list1,
            zhe_tc_list1,
            plot_contour_fills,
            cmap=plot_colormap,
            alpha=plot_alpha,
        )

        # This is the fix for the white lines between contour levels
        for c in zr_contourf_tc1.collections:
            c.set_edgecolor("face")

        # --- Zircon cooling rate versus eU plot ---
        # Create closure temperature contour lines
        zr_contours_tc = ax[1].tricontour(
            zr_x_list2,
            zr_y_list2,
            zhe_tc_list2,
            plot_contour_lines,
            linewidths=0.5,
            colors="k",
        )
        # Use log x-axis scaling
        ax[1].set_xscale("log")
        # Add closure temperature contour labels
        ax[1].clabel(zr_contours_tc)
        # Create closure temperature contour fill
        zr_contourf_tc2 = ax[1].tricontourf(
            zr_x_list2,
            zr_y_list2,
            zhe_tc_list2,
            plot_contour_fills,
            cmap=plot_colormap,
            alpha=plot_alpha,
        )

        # This is the fix for the white lines between contour levels
        for c in zr_contourf_tc2.collections:
            c.set_edgecolor("face")

    # Plot values for apatite and zircon (U-Th)/He
    else:
        # --- Apatite cooling rate versus radius ---
        # Create closure temperature contour lines
        ap_contours_tc = ax[0][0].tricontour(
            ap_x_list1,
            ap_y_list1,
            ahe_tc_list1,
            plot_contour_lines,
            linewidths=0.5,
            colors="k",
        )
        # Use log x-axis scaling
        ax[0][0].set_xscale("log")
        # Add closure temperature contour labels
        ax[0][0].clabel(ap_contours_tc)
        # Create closure temperature contour fill
        ap_contourf_tc1 = ax[0][0].tricontourf(
            ap_x_list1,
            ap_y_list1,
            ahe_tc_list1,
            plot_contour_fills,
            cmap=plot_colormap,
            alpha=plot_alpha,
        )

        # This is the fix for the white lines between contour levels
        for c in ap_contourf_tc1.collections:
            c.set_edgecolor("face")

        # --- Apatite cooling rate versus eU plot ---
        # Create closure temperature contour lines
        ap_contours_tc = ax[0][1].tricontour(
            ap_x_list2,
            ap_y_list2,
            ahe_tc_list2,
            plot_contour_lines,
            linewidths=0.5,
            colors="k",
        )
        # Use log x-axis scaling
        ax[0][1].set_xscale("log")
        # Add closure temperature contour labels
        ax[0][1].clabel(ap_contours_tc)
        # Create closure temperature contour fill
        ap_contourf_tc2 = ax[0][1].tricontourf(
            ap_x_list2,
            ap_y_list2,
            ahe_tc_list2,
            plot_contour_fills,
            cmap=plot_colormap,
            alpha=plot_alpha,
        )

        # This is the fix for the white lines between contour levels
        for c in ap_contourf_tc2.collections:
            c.set_edgecolor("face")

        # --- Zircon cooling rate versus radius plot ---
        # Create closure temperature contour lines
        zr_contours_tc = ax[1][0].tricontour(
            zr_x_list1,
            zr_y_list1,
            zhe_tc_list1,
            plot_contour_lines,
            linewidths=0.5,
            colors="k",
        )
        # Use log x-axis scaling
        ax[1][0].set_xscale("log")
        # Add closure temperature contour labels
        ax[1][0].clabel(zr_contours_tc)
        # Create closure temperature contour fill
        zr_contourf_tc1 = ax[1][0].tricontourf(
            zr_x_list1,
            zr_y_list1,
            zhe_tc_list1,
            plot_contour_fills,
            cmap=plot_colormap,
            alpha=plot_alpha,
        )

        # This is the fix for the white lines between contour levels
        for c in zr_contourf_tc1.collections:
            c.set_edgecolor("face")

        # --- Zircon cooling rate versus eU plot ---
        # Create closure temperature contour lines
        zr_contours_tc = ax[1][1].tricontour(
            zr_x_list2,
            zr_y_list2,
            zhe_tc_list2,
            plot_contour_lines,
            linewidths=0.5,
            colors="k",
        )
        # Use log x-axis scaling
        ax[1][1].set_xscale("log")
        # Add closure temperature contour labels
        ax[1][1].clabel(zr_contours_tc)
        # Create closure temperature contour fill
        zr_contourf_tc2 = ax[1][1].tricontourf(
            zr_x_list2,
            zr_y_list2,
            zhe_tc_list2,
            plot_contour_fills,
            cmap=plot_colormap,
            alpha=plot_alpha,
        )

        # This is the fix for the white lines between contour levels
        for c in zr_contourf_tc2.collections:
            c.set_edgecolor("face")

    # Format plot

    # Apatite only
    if plot_type == 1:
        ax[0].set_title("Apatite (U-Th)/He closure temperature [°C]")
        ax[1].set_title("Apatite (U-Th)/He closure temperature [°C]")
    # Zircon only
    elif plot_type == 2:
        ax[0].set_title("Zircon (U-Th)/He closure temperature [°C]")
        ax[1].set_title("Zircon (U-Th)/He closure temperature [°C]")
    # Apatite and zircon
    else:
        ax[0][0].set_title("Apatite (U-Th)/He closure temperature [°C]")
        ax[0][1].set_title("Apatite (U-Th)/He closure temperature [°C]")
        ax[1][0].set_title("Zircon (U-Th)/He closure temperature [°C]")
        ax[1][1].set_title("Zircon (U-Th)/He closure temperature [°C]")

    # Apatite or Zircon
    if plot_type < 3:
        ax[0].set_xlabel("Cooling rate [°C/Myr]")
        ax[1].set_xlabel("Cooling rate [°C/Myr]")
        ax[0].set_ylabel("Equivalent spherical radius (µm)")
        ax[1].set_ylabel("Effective uranium (ppm)")

    # Apatite and zircon eU versus radius
    else:
        ax[0][0].set_xlabel("Cooling rate [°C/Myr]")
        ax[0][1].set_xlabel("Cooling rate [°C/Myr]")
        ax[0][0].set_ylabel("Equivalent spherical radius (µm)")
        ax[0][1].set_ylabel("Effective uranium (ppm)")
        ax[1][0].set_xlabel("Cooling rate [°C/Myr]")
        ax[1][1].set_xlabel("Cooling rate [°C/Myr]")
        ax[1][0].set_ylabel("Equivalent spherical radius (µm)")
        ax[1][1].set_ylabel("Effective uranium (ppm)")

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
        # Set file name prefix
        plot_filename = "rate_vs_radius_eu"

        # Define plot filename based on type of plot and save plot
        if plot_type == 1:
            plot_savename = (
                plot_filename + "_apatite_" + str(plot_dpi) + "dpi." + plot_file_format
            )
        elif plot_type == 2:
            plot_savename = (
                plot_filename + "_zircon_" + str(plot_dpi) + "dpi." + plot_file_format
            )
        else:
            plot_savename = (
                plot_filename
                + "_apatite_zircon_"
                + str(plot_dpi)
                + "dpi."
                + plot_file_format
            )
        plt.savefig(wd_orig + "/" + plot_savename, dpi=plot_dpi)

    # Save plot if requested
    if display_plot:
        plt.show()

    # Revert to original working directory
    os.chdir(wd_orig)

    return None


# Define function for creating plot of cooling rates
def rate_vs_age_tc(
    num_points=101,
    cooling_rate_min=0.1,
    cooling_rate_max=100.0,
    temp_max=350.0,
    ap_u1=1.0,
    ap_u2=20.0,
    ap_u3=150.0,
    zr_u1=10.0,
    zr_u2=200.0,
    zr_u3=4000.0,
    ap_rad=45.0,
    zr_rad=60.0,
    ap_thorium=0.0,
    zr_thorium=0.0,
    ahe_uncertainty=0.1,
    aft_uncertainty=0.2,
    zhe_uncertainty=0.1,
    plot_type=3,
    plot_age_min=0.5,
    plot_age_max=1800.0,
    plot_tc_min=0.0,
    plot_tc_max=200.0,
    save_plot=False,
    plot_file_format="pdf",
    plot_dpi=300,
    plot_style="seaborn-colorblind",
    plot_alpha=0.6,
    plot_grid=True,
    display_plot=True,
    clean_up_files=True,
    verbose=False,
    use_widget=False,
):
    """
    Calculates thermochronometer ages and closure temperatures for different cooling rates and effective uranium
    concentrations.

    Parameters
    ----------
    num_points : int, default=101
        Number of points along x and y axes where ages/closure temperatures are
        calculated.
    cooling_rate_min : float, default=0.1
        Minimum cooling rate in degrees C per Myr.
    cooling_rate_max : float, default=100.0
        Maximum cooling rate in degrees C per Myr.
    temp_max : float, default=350.0
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
    plot_age_min : float, default=0.5
        Minimum age value in Ma for plots of cooling rate versus age. Only applies to plot_type 2 and 3.
    plot_age_max : float, default=1800.0
        Maximum age value in Ma for plots of cooling rate versus age. Only applies to plot_type 2 and 3.
    plot_tc_min : float, default=0.0
        Minimum closure temperature value in deg. C for plots of cooling rate versus closure temperature.
        Only applies to plot_type 1 and 3.
    plot_tc_max : float, default=200.0
        Maximum closure temperature value in deg. C for plots of cooling rate versus closure temperature.
        Only applies to plot_type 1 and 3.
    save_plot : bool, default=False
        Flag for whether to save the plot to a file.
    plot_file_format : str, default='pdf'
        File format for saving plot to file (examples: png, pdf, svg, eps).
    plot_dpi : int, default=300
        Saved plot resolution in dots per inch.
    plot_style : str, default='seaborn-colorblind'
        Style sheet used for plotting. See https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html.
    plot_alpha : float, default=0.6
        Transparency used for plotting fill colors for age swath plots.
    plot_grid : bool, default=True
        Flag for whether or not to display the plot grid lines.
    display_plot : bool, default=True
        Flag for whether to display the plot.
    clean_up_files : bool, default=True
        Flag for whether to delete temporary output files after the code has run.
    verbose : bool, default=False
        Enable/disable verbose output.
    use_widget : bool, default=False
        Enable/disable IPython progress bar widget. Disabled for command-line usage.

    Returns
    -------
    None

    """

    # Check to see whether ipywidgets and IPython are available for widget use
    # If not, disable widgets and display a warning
    if use_widget:
        try:
            import ipywidgets as widgets
        except ModuleNotFoundError:
            print("Warning: ipywidgets module not found. Disabling graphical progress bar.")
            use_widget = False
    if use_widget:
        try:
            from IPython.display import display
        except ModuleNotFoundError:
            print(
                "Warning: IPython.display module not found. Disabling graphical progress bar."
            )
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
        plot_filename = "rate_vs_tc"
    elif plot_type == 2:
        plot_filename = "rate_vs_age"
    elif plot_type == 3:
        plot_filename = "rate_vs_age_tc"
    else:
        raise ValueError("Bad value for plot_type. Must be 1, 2, or 3.")

    # Define cooling rates to consider
    rates = np.logspace(
        start=np.log10(cooling_rate_min),
        stop=np.log10(cooling_rate_max),
        num=num_points,
    )

    # Plot titles
    title_list = [
        f"Low eU (ap={ap_u_list[0]:.1f}, zr={zr_u_list[0]:.1f} ppm)",
        f"Intermediate eU (ap={ap_u_list[1]:.1f}, zr={zr_u_list[1]:.1f} ppm)",
        f"High eU (ap={ap_u_list[2]:.1f}, zr={zr_u_list[2]:.1f} ppm)",
    ]

    # Define time-temperature history filename
    tt_file = "simple_time_temp.txt"

    # Get age calculation executable(s) to use
    rdaam_command = get_tc_exec("RDAAM_He")
    ketch_command = get_tc_exec("ketch_aft")

    # Calculate total number of models that will be run
    total_models = len(ap_u_list) * len(rates)

    # Set model type string
    if plot_type == 1:
        model_type = "cooling rate versus closure temperature"
    elif plot_type == 2:
        model_type = "cooling rate versus age"
    elif plot_type == 3:
        model_type = "cooling rate versus age and closure temperature"

    # Set plot style
    plt.style.use(plot_style)

    # Create figure
    if plot_type == 3:
        fig, ax = plt.subplots(3, 2, figsize=(12, 10))
    else:
        fig, ax = plt.subplots(3, 1, figsize=(6, 10))

    # Create visual progress bar, if enabled
    if use_widget and not verbose:
        s = widgets.IntProgress(
            value=0,
            min=0,
            max=total_models,
            description="Calculating:",
            bar_style="",  # 'success', 'info', 'warning', 'danger' or ''
            style={"bar_color": "#ff6666"},
            orientation="horizontal",
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
                        f"Calculating {model_type} - {int(round(100 * model_count / total_models)):3d}% ({model_count:5d} / {total_models:5d})\r",
                        end="",
                    )

            # Define thermal history
            start_time = temp_max / rate
            with open(tt_file, "w") as f:
                f.write("0.0,0.0\n")
                f.write("{0:.4f},{1:.1f}".format(start_time, temp_max))

            # Calculate He ages
            command = (
                rdaam_command
                + " "
                + tt_file
                + " "
                + str(ap_rad)
                + " "
                + str(ap_uranium)
                + " "
                + str(ap_thorium)
                + " "
                + str(zr_rad)
                + " "
                + str(zr_uranium)
                + " "
                + str(zr_thorium)
            )
            p = subprocess.Popen(
                command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
            )

            # Parse output for ages
            stdout = p.stdout.readlines()
            corr_ahe_age = stdout[0].split()[7].decode("UTF-8")
            corr_zhe_age = stdout[1].split()[7].decode("UTF-8")

            # Calculate AFT age
            command = ketch_command + " " + tt_file
            p = subprocess.Popen(
                command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
            )

            # Parse output for AFT age
            stdout = p.stdout.readlines()
            aft_age = stdout[0].split()[4][:-1].decode("UTF-8")

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
                    f"AHe: {float(corr_ahe_age):.2f} Ma (Tc: {ahe_tc:.1f}°C); AFT: {float(aft_age):.2f} Ma (Tc: {aft_tc:.1f}°C); ZHe: {float(corr_zhe_age):.2f} Ma (Tc: {zhe_tc:.1f}°C) -- total time: {start_time:.1f} Myr"
                )

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
            ax[i].semilogx(rate_list, ahe_tc_list, label="Apatite (U-Th)/He")
            ax[i].semilogx(rate_list, aft_tc_list, label="Apatite FT")
            ax[i].semilogx(rate_list, zhe_tc_list, label="Zircon (U-Th)/He")

        # Create plots for rate versus age
        if plot_type == 2:
            ax[i].fill_between(
                rate_list,
                ahe_age_min,
                ahe_age_max,
                alpha=plot_alpha,
                label=f"Apatite (U-Th)/He age ± {ahe_uncertainty * 100:.0f}%",
            )
            ax[i].fill_between(
                rate_list,
                aft_age_min,
                aft_age_max,
                alpha=plot_alpha,
                label=f"Apatite FT age ± {aft_uncertainty * 100:.0f}%",
            )
            ax[i].fill_between(
                rate_list,
                zhe_age_min,
                zhe_age_max,
                alpha=plot_alpha,
                label=f"Zircon (U-Th)/He age ± {zhe_uncertainty * 100:.0f}%",
            )

            # Scale axes
            ax[i].set_xscale("log")
            ax[i].set_yscale("log")

        # Create plots for rate versus age and closure temperature
        if plot_type == 3:
            # Plot ages and closure temperatures (low eU)
            ax[i][0].fill_between(
                rate_list,
                ahe_age_min,
                ahe_age_max,
                alpha=plot_alpha,
                label=f"Apatite (U-Th)/He age ± {ahe_uncertainty * 100:.0f}%",
            )
            ax[i][1].plot(rate_list, ahe_tc_list, label="Apatite (U-Th)/He")

            # Plot ages and closure temperatures (intermediate eU)
            ax[i][0].fill_between(
                rate_list,
                aft_age_min,
                aft_age_max,
                alpha=plot_alpha,
                label=f"Apatite FT age ± {aft_uncertainty * 100:.0f}%",
            )
            ax[i][1].plot(rate_list, aft_tc_list, label="Apatite FT")

            # Plot ages and closure temperatures (high eU)
            ax[i][0].fill_between(
                rate_list,
                zhe_age_min,
                zhe_age_max,
                alpha=plot_alpha,
                label=f"Zircon (U-Th)/He age ± {zhe_uncertainty * 100:.0f}%",
            )
            ax[i][1].plot(rate_list, zhe_tc_list, label="Zircon (U-Th)/He")

            # Set axis scalings
            ax[i][0].set_xscale("log")
            ax[i][0].set_yscale("log")
            ax[i][1].set_xscale("log")

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
            ax[i].set_xlim([cooling_rate_min, cooling_rate_max])
            ax[i].set_ylim([plot_tc_min, plot_tc_max])
            ax[i].set_ylabel("Closure temperature (°C)")
            if i == 2:
                ax[i].set_xlabel("Cooling rate (°C/Myr)")

        # Set axis range and add axis labels
        if plot_type == 2:
            ax[i].set_xlim([cooling_rate_min, cooling_rate_max])
            ax[i].set_ylim([plot_age_min, plot_age_max])
            ax[i].set_ylabel("Age (Ma)")
            if i == 2:
                ax[i].set_xlabel("Cooling rate (°C/Myr)")

        # Set axis ranges and add axis labels
        if plot_type == 3:
            ax[i][0].set_xlim([cooling_rate_min, cooling_rate_max])
            ax[i][0].set_ylim([plot_age_min, plot_age_max])
            ax[i][1].set_xlim([cooling_rate_min, cooling_rate_max])
            ax[i][1].set_ylim([plot_tc_min, plot_tc_max])
            ax[i][0].set_ylabel("Age (Ma)")
            ax[i][1].set_ylabel("Closure temperature (°C)")
            if i == 2:
                ax[i][0].set_xlabel("Cooling rate (°C/Myr)")
                ax[i][1].set_xlabel("Cooling rate (°C/Myr)")

        # Add subplot titles
        if plot_type == 3:
            ax[i][0].set_title(title_list[i])
            ax[i][1].set_title(title_list[i])
        else:
            ax[i].set_title(title_list[i])

        # Enable/disable gridlines
        if plot_grid:
            if plot_type == 3:
                ax[i][0].grid(visible=True)
                ax[i][1].grid(visible=True)
            else:
                ax[i].grid(visible=True)
        else:
            if plot_type == 3:
                ax[i][0].grid(visible=False)
                ax[i][1].grid(visible=False)
            else:
                ax[i].grid(visible=False)

        # Add legend
        if plot_type == 3:
            ax[i][0].legend()
            ax[i][1].legend()
        else:
            ax[i].legend()

    # Delete temporary tt file
    if clean_up_files:
        os.remove(tt_file)
        os.remove("ft_length.csv")

    # Use tight layout
    plt.tight_layout()

    # Save plot if requested
    if save_plot:
        # Define plot filename and save plot
        plot_filename = plot_filename + "_" + str(plot_dpi) + "dpi." + plot_file_format
        plt.savefig(wd_orig + "/" + plot_filename, dpi=plot_dpi)

    # Show plot if requested
    if display_plot:
        plt.show()

    # Revert to original working directory
    os.chdir(wd_orig)

    return None
