import re
from pathlib import Path

import pytest

import tcplotter
from tcplotter.tcplotter import calc_eu

DATA_DIR = Path(__file__).resolve().parents[1] / "data"

# Ages from rate_vs_age_tc(num_points=5, plot_type=3, verbose=True), as calculated
# by the RDAAM_He and ketch_aft executables used by tcplotter 0.4.0 (2 decimals).
# Rows: low, intermediate and high eU; each row is (AHe, AFT, ZHe) for 5 cooling rates.
REFERENCE_AGES = [
    (317.15, 716.28, 1313.11),
    (67.92, 149.78, 209.10),
    (15.59, 30.87, 33.97),
    (3.45, 6.27, 5.50),
    (0.74, 1.26, 0.92),
    (595.03, 716.28, 1468.25),
    (91.93, 149.78, 279.64),
    (15.97, 30.87, 47.12),
    (3.45, 6.27, 7.80),
    (0.74, 1.26, 1.27),
    (763.16, 716.28, 0.00),
    (128.91, 149.78, 140.03),
    (20.69, 30.87, 58.31),
    (3.55, 6.27, 10.37),
    (0.74, 1.26, 1.74),
]


def test_calc_eu():
    assert calc_eu(10.0, 40.0) == pytest.approx(19.52)


def test_rate_vs_age_tc_matches_reference(capsys):
    """Predicted ages match those from the original age prediction executables."""
    tcplotter.rate_vs_age_tc(num_points=5, plot_type=3, display_plot=False, verbose=True)
    out = capsys.readouterr().out
    ages = re.findall(r"AHe: ([\d.]+) Ma.*AFT: ([\d.]+) Ma.*ZHe: ([\d.]+) Ma", out)
    assert len(ages) == len(REFERENCE_AGES)
    for predicted, expected in zip(ages, REFERENCE_AGES):
        assert [float(a) for a in predicted] == pytest.approx(expected, abs=0.011)


def test_time_vs_temp():
    tcplotter.time_vs_temp(display_plot=False)


@pytest.mark.parametrize("plot_type", [1, 2, 3])
def test_eu_vs_radius_plot_types(plot_type):
    tcplotter.eu_vs_radius(num_points=3, plot_type=plot_type, display_plot=False)


def test_eu_vs_radius_custom_history():
    tcplotter.eu_vs_radius(
        num_points=3,
        cooling_hist_type=2,
        time_hist=[0.0, 10.0, 30.0],
        temp_hist=[0.0, 80.0, 200.0],
        display_plot=False,
        tt_plot=True,
    )


def test_eu_vs_radius_history_from_file(tmp_path):
    tt_file = tmp_path / "tt.csv"
    tt_file.write_text("Time (Ma),Temperature (C)\n0.0,10.0\n5.0,60.0\n20.0,90.0\n40.0,250.0\n")
    tcplotter.eu_vs_radius(
        num_points=3, cooling_hist_type=3, input_tt_file="tt.csv", display_plot=False
    )


def test_eu_vs_radius_misfit(capsys):
    tcplotter.eu_vs_radius(
        num_points=3,
        age_data_file=str(DATA_DIR / "sample_data.csv"),
        calc_misfit=True,
        display_plot=False,
    )
    assert "Total misfit" in capsys.readouterr().out


@pytest.mark.parametrize("plot_type", [1, 2, 3])
def test_rate_vs_radius_eu(plot_type):
    tcplotter.rate_vs_radius_eu(num_points=3, plot_type=plot_type, display_plot=False)


@pytest.mark.parametrize("plot_type", [1, 2, 3])
def test_rate_vs_age_tc_plot_types(plot_type):
    tcplotter.rate_vs_age_tc(num_points=3, plot_type=plot_type, display_plot=False)


def test_no_temporary_files_left(tmp_path):
    tcplotter.eu_vs_radius(num_points=3, display_plot=False)
    tcplotter.rate_vs_radius_eu(num_points=3, display_plot=False)
    tcplotter.rate_vs_age_tc(num_points=3, display_plot=False)
    assert list(tmp_path.iterdir()) == []


def test_save_plot(tmp_path):
    tcplotter.rate_vs_age_tc(
        num_points=3, display_plot=False, save_plot=True, plot_file_format="png", plot_dpi=50
    )
    assert len(list((tmp_path / "png").glob("*.png"))) == 1


def test_clean_up_files_is_deprecated():
    with pytest.warns(DeprecationWarning, match="clean_up_files"):
        tcplotter.rate_vs_age_tc(num_points=3, display_plot=False, clean_up_files=False)

def test_version():
    assert tcplotter.__version__ not in ("", "dev")
