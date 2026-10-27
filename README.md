# T<sub>c</sub>plotter code description and instructions for use

[![DOI](https://zenodo.org/badge/444915688.svg)](https://zenodo.org/badge/latestdoi/444915688)
[![PyPI - Version](https://img.shields.io/pypi/v/tcplotter)](https://pypi.org/project/tcplotter/)
[![Build, test, and deploy](https://github.com/HUGG/tcplotter/actions/workflows/publish-to-pypi.yml/badge.svg)](https://github.com/HUGG/tcplotter/actions/workflows/publish-to-pypi.yml)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/tcplotter/badges/version.svg)](https://anaconda.org/conda-forge/tcplotter)
[![Documentation Status](https://readthedocs.org/projects/tcplotter/badge/?version=latest)](https://tcplotter.readthedocs.io/en/latest/?badge=latest)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/HUGG/tcplotter/HEAD?urlpath=lab/tree/notebooks/tcplotter.ipynb)

T<sub>c</sub>plotter is a Python package for creating and customizing thermochronometer age and closure temperature plots presented in the article [Short communication: Modelling competing effects of cooling rate, grain size and radiation damage in low temperature thermochronometers](https://doi.org/10.5194/gchron-2021-29) by D. Whipp, D. Kellett, I. Coutand, and R. Ketcham.
The code is designed to be easy to use to either reproduce plots from the article or customize the plots for your own use.
Below you will find some essential details about using the code and detailed documentation can be found on the [tcplotter documentation site](https://tcplotter.readthedocs.io/).

## Getting started

The easiest way to get started using T<sub>c</sub>plotter is by [using Binder](https://mybinder.org/v2/gh/HUGG/tcplotter/HEAD?urlpath=lab/tree/notebooks/tcplotter.ipynb) either with the link here or by clicking on the Binder button above.
This will open an interactive web interface to the code where you can see how it works and even create your own plots.

## Installation

T<sub>c</sub>plotter can be installed in your Python environment using `pip`.

```bash
pip install tcplotter
```

You can also find installation instructions on the [T<sub>c</sub>plotter documentation page](https://tcplotter.readthedocs.io/en/latest/installation.html#installing-the-latest-version-of-t-sub-c-sub-plotter-from-github).

Installing T<sub>c</sub>plotter also installs [T<sub>c</sub>ages](https://github.com/HUGG/tcages), which provides the thermochronometer age prediction codes. Pre-built packages are available for Linux, macOS and Windows, so no compilers are needed.

**Note**: T<sub>c</sub>plotter versions up to 0.4.0 required installing separate age prediction programs from the T<sub>c</sub>core package. This is no longer necessary from version 0.5.0 onwards.

## Usage

T<sub>c</sub>plotter can be used either as a function in a Python script or interpreter, or from the command line.
The four main T<sub>c</sub>plotter functions/command-line tools are:

- `time_vs_temp()`
- `eu_vs_radius()`
- `rate_vs_radius_eu()`
- `rate_vs_age_tc()`

Brief examples of possible usage for both cases can be found below.

### Usage in a Python script or interpreter

Assuming you are working in Binder or have installed `tcplotter` locally, functions available in T<sub>c</sub>plotter can be imported as follows:

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
For example, again assuming you have installed `tcplotter` locally, you can type the following to use the `eu_vs_radius()` function:

```bash
eu-vs-radius --save-plot
```

To find more information about options available for command-line use you can include the `--help` or `-h` flags.

```bash
eu-vs-radius -h
```

## Attribution

### How to cite the code

When including plots generated using T<sub>c</sub>plotter in publications or presentations, please cite the following:

- Whipp, D. M. and Ketcham, R. A.: tcplotter: a Python package for creating and customizing thermochronometer age and closure temperature plots, Zenodo [code], <https://doi.org/10.5281/zenodo.5958939>, 2022.

### How to cite the paper

You're also welcome to cite the manuscript related to T<sub>c</sub>plotter when referencing ideas from the article:

- Whipp, D. M., Kellett, D. A., Coutand, I., and Ketcham, R. A.: Short communication: Modeling competing effects of cooling rate, grain size, and radiation damage in low-temperature thermochronometers, *Geochronology*, *4*, 143–152, <https://doi.org/10.5194/gchron-4-143-2022>, 2022. 

### How to cite related articles

The age prediction software used for calculating apatite and zircon (U-Th)/He and apatite fission-track ages (provided by [T<sub>c</sub>ages](https://github.com/HUGG/tccore)) was written by Richard Ketcham at the University of Texas, USA. Results published using this software should also cite the articles below:

- Flowers, R. M., Ketcham, R. A., Shuster, D. L., and Farley, K. A.: Apatite (U-Th)/He thermochronometry using a radiation damage accumulation and annealing model, Geochimica et Cosmochimica Acta, 73, 2347–2365, doi: [10.1016/j.gca.2009.01.015](https://doi.org/10.1016/j.gca.2009.01.015), 2009.

- Guenthner, W. R., Reiners, P. W., Ketcham, R. A., Nasdala, L., and Giester, G.: Helium diffusion in natural zircon: Radiation damage, anisotropy, and the interpretation of zircon (U-Th)/He thermochronology, American Journal of Science, 313, 145–198, doi: [10.2475/03.2013.01](https://doi.org/10.2475/03.2013.01), 2013.

- Ketcham, R. A., Donelick, R. A., & Carlson, W. D.: Variability of apatite fission-track annealing kinetics III: Extrapolation to geological time scales. American Mineralogist, 84, 1235-1255, doi: [10.2138/am-1999-0903](https://doi.org/10.2138/am-1999-0903), 1999.

- Ketcham, R. A., Mora, A., and Parra, M.: Deciphering exhumation and burial history with multi-sample down-well thermochronometric inverse modelling, Basin Res., 30, 48-64, [10.1111/bre.12207](https://doi.org/10.1111/bre.12207), 2018.

## License

T<sub>c</sub>plotter is distributed under the [GNU General Public License v3.0 or later](LICENSE).
Versions up to 0.4.0 were released under the MIT License.
