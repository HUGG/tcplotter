from .tcplotter import time_vs_temp
from .tcplotter import eu_vs_radius
from .tcplotter import rate_vs_radius_eu
from .tcplotter import rate_vs_age_tc

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("tcplotter")
except PackageNotFoundError:
    # Package is not installed (e.g. running from a source checkout)
    __version__ = "dev"
