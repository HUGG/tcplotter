import matplotlib

# Use a non-interactive backend so tests run without a display
matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402


@pytest.fixture(autouse=True)
def run_in_tmp_path(tmp_path, monkeypatch):
    """Run each test in an empty temporary directory and close figures afterwards."""
    monkeypatch.chdir(tmp_path)
    yield tmp_path
    plt.close("all")
