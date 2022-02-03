# T<sub>c</sub>plotter code description and instructions for use

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/HUGG/tcplotter/HEAD?urlpath=lab/tree/tcplotter.ipynb)
[![Documentation Status](https://readthedocs.org/projects/tcplotter/badge/?version=latest)](https://tcplotter.readthedocs.io/en/latest/?badge=latest) 
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/444915688.svg)](https://zenodo.org/badge/latestdoi/444915688)

T<sub>c</sub>plotter is a Python package for creating and customizing thermochronometer age and closure temperature plots presented in the article [Modelling competing effects of cooling rate, grain size and radiation damage in low temperature thermochronometers](https://gchron.copernicus.org/preprints/gchron-2021-29/) by D. Whipp, D. Kellett, I. Coutand, and R. Ketcham.
The code is designed to be easy to use to either reproduce plots from the article or customize the plots for your own use.
Below you will find some essential details about using the code and detailed documentation can be found on the [tcplotter documentation site](https://tcplotter.readthedocs.io/).

## Getting started

The easiest way to get started using T<sub>c</sub>plotter is by [using Binder](https://mybinder.org/v2/gh/HUGG/tcplotter/HEAD?urlpath=lab/tree/tcplotter.ipynb) either with the link here or by clicking on the Binder button above.
This will open an interactive web interface to the code where you can see how it works and even create your own plots.

## Installation

Currently, we recommend using the code [via Binder](https://mybinder.org/v2/gh/HUGG/tcplotter/HEAD?urlpath=lab/tree/tcplotter.ipynb) as the easiest option, as we do not yet have a Python package available for T<sub>c</sub>plotter.
If you would like to install the software for your own use, you can find detailed instructions on the [T<sub>c</sub>plotter documentation page](https://tcplotter.readthedocs.io/en/latest/installation.html#installing-the-latest-version-of-t-sub-c-sub-plotter-from-github).

## Usage

T<sub>c</sub>plotter can be used either as a function in a Python script or interpreter, or from the command line.
The four main T<sub>c</sub>plotter functions/command-line tools are:

- `time_vs_temp`
- `eu_vs_radius`
- `rate_vs_radius_eu`
- `rate_vs_age_tc`

Brief examples of possible usage for both cases can be found below.

### Usage in a Python script or interpreter

Assuming the base `tcplotter` directory is your working directory, functions available in T<sub>c</sub>plotter can be imported as follows:

```python
from tcplotter import time_vs_temp, eu_vs_radius, rate_vs_radius_eu, rate_vs_age_tc
```

Once imported, you can use functions as shown below:

```python
eu_versus_radius(save_plot=True)
```

You can find more information about the function parameters using the `help()` function:

```python
help(eu_vs_radius)
```

### Command-line usage

Command-line usage is similar to that for use in a Python script, except that the underscores in the function names have been replaced by hyphens.
For example, assuming you are in a terminal in the `tcplotter/tcplotter` directory, you can type the following to use the `eu_vs_radius()` function:

```bash
./eu-vs-radius --save-plot
```

To find more information about options available for command-line use you can include the `--help` or `-h` flags.

```bash
./eu-vs-radius -h
```

## Attribution

When including plots generated using T<sub>c</sub>plotter in publications or presentations, please cite the following article:

- Whipp, D. M., Kellett, D. A., Coutand, I, & Ketcham, R. A. (in revision). Modelling competing effects of cooling rate, grain size and radiation damage in low temperature thermochronometers. *Geochronology*, doi: [10.5194/gchron-2021-29](https://doi.org/10.5194/gchron-2021-29).

The age prediction software used for calculating apatite and zircon (U-Th)/He and apatite fission-track ages was written by Richard Ketcham at the University of Texas, USA. Results published using this software should cite the articles below:

- Ketcham, R. A., Donelick, R. A., & Carlson, W. D.: Variability of apatite fission-track annealing kinetics III: Extrapolation to geological time scales. American Mineralogist, 84, 1235-1255, doi: [10.2138/am-1999-0903](https://doi.org/10.2138/am-1999-0903), 1999.

- Ketcham, R. A., Mora, A., and Parra, M.: Deciphering exhumation and burial history with multi-sample down-well thermochronometric inverse modelling, Basin Res., 30, 48-64, [10.1111/bre.12207](https://doi.org/10.1111/bre.12207), 2018.

## License

The T<sub>c</sub>plotter software is licensed under an MIT License: [T<sub>c</sub>plotter software license](LICENSE)