"""Integration tests for ThermophysicalProperties using the sample case files."""

import os
import shutil

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest

matplotlib.use("Agg")

# Skip the entire module if PyFoam is not importable (broken install, etc.)
pytest.importorskip("PyFoam")

from FoamUtils.ThermophysicalProperties import (
    ThermophysicalProperties,
    plot_thermophysical_properties,
)

# Path to the sample thermophysicalProperties file bundled with the repo.
REPO_ROOT = os.path.join(os.path.dirname(__file__), "..")
CASE_DIR = os.path.join(REPO_ROOT, "constant")
THERMO_FILE = os.path.join(CASE_DIR, "thermophysicalProperties")
PHYSICAL_FILE = os.path.join(CASE_DIR, "physicalProperties")


@pytest.fixture(scope="module")
def thermo():
    """Create a ThermophysicalProperties object from the sample case."""
    if not os.path.exists(THERMO_FILE):
        pytest.skip("Sample thermophysicalProperties file not found.")
    return ThermophysicalProperties(THERMO_FILE)


class TestThermophysicalPropertiesInit:
    def test_loads_without_error(self, thermo):
        assert thermo is not None

    def test_repr(self, thermo):
        assert "ThermophysicalProperties" in repr(thermo)

    def test_molar_mass_positive(self, thermo):
        assert thermo.M > 0

    def test_gas_constant_positive(self, thermo):
        assert thermo.R > 0

    def test_file_not_found_raises(self):
        with pytest.raises(FileNotFoundError):
            ThermophysicalProperties("/nonexistent/path/thermophysicalProperties")


class TestThermophysicalPropertiesMethods:
    """Test the callable property functions on a representative temperature range."""

    T_RANGE = np.linspace(21.0, 30.0, 10)  # LH2 boiling-range temperatures (K)
    P = 1e5  # Pa

    def test_rho_positive(self, thermo):
        rho = thermo.rho(self.P, self.T_RANGE)
        assert np.all(rho > 0), "Density must be positive everywhere."

    def test_Cp_positive(self, thermo):
        Cp = thermo.Cp(self.T_RANGE)
        assert np.all(Cp > 0), "Heat capacity must be positive everywhere."

    def test_mu_positive(self, thermo):
        mu = thermo.mu(self.T_RANGE)
        assert np.all(mu > 0), "Dynamic viscosity must be positive everywhere."

    def test_kappa_positive(self, thermo):
        kappa = thermo.kappa(self.T_RANGE)
        assert np.all(kappa > 0), "Thermal conductivity must be positive everywhere."

    def test_beta_positive(self, thermo):
        beta = thermo.beta(self.P, self.T_RANGE)
        assert np.all(beta > 0), "Thermal expansion coefficient must be positive."

    def test_alpha_positive(self, thermo):
        alpha = thermo.alpha(self.P, self.T_RANGE)
        assert np.all(alpha > 0), "Thermal diffusivity must be positive."

    def test_unit_label_rho(self, thermo):
        label = thermo.unit_label("rho")
        assert "rho" in label or "ρ" in label

    def test_unit_label_Cp(self, thermo):
        label = thermo.unit_label("Cp")
        assert label is not None


class TestPlotThermo:
    """Test that plot_thermo can parse the properties file and produce a plot."""

    T_RANGE = np.linspace(21.0, 30.0, 10)  # LH2 boiling-range temperatures (K)
    P = 1e5  # Pa
    PROPERTIES = ["rho", "Cp", "mu", "kappa"]

    def test_plot_produces_lines_on_axes(self, thermo):
        """Verify that plot() draws lines onto the provided axes."""
        fig, axes = plt.subplots(1, len(self.PROPERTIES))
        thermo.plot(self.T_RANGE, self.P, self.PROPERTIES, axes=axes)
        for ax in axes:
            assert len(ax.lines) > 0, "Expected each axis to contain at least one line."
        plt.close(fig)

    def test_plot_thermophysical_properties_uses_thermo_file(self, tmp_path):
        """plot_thermophysical_properties reads thermophysicalProperties when present."""
        if not os.path.exists(THERMO_FILE):
            pytest.skip("Sample thermophysicalProperties file not found.")
        case_constant = tmp_path / "constant"
        case_constant.mkdir()
        shutil.copy(THERMO_FILE, case_constant / "thermophysicalProperties")
        shutil.copy(os.path.join(REPO_ROOT, "parameters"), tmp_path / "parameters")
        plot_thermophysical_properties(
            self.P, self.T_RANGE[0], self.T_RANGE[-1], self.PROPERTIES, case_dir=str(tmp_path)
        )
        plt.close("all")

    def test_plot_thermophysical_properties_falls_back_to_physical_properties(self, tmp_path):
        """plot_thermophysical_properties falls back to physicalProperties (OpenFOAM.org >= v11)."""
        if not os.path.exists(PHYSICAL_FILE):
            pytest.skip("Sample physicalProperties file not found.")
        case_constant = tmp_path / "constant"
        case_constant.mkdir()
        # Only place physicalProperties — no thermophysicalProperties
        shutil.copy(PHYSICAL_FILE, case_constant / "physicalProperties")
        shutil.copy(os.path.join(REPO_ROOT, "parameters"), tmp_path / "parameters")
        plot_thermophysical_properties(
            self.P, self.T_RANGE[0], self.T_RANGE[-1], self.PROPERTIES, case_dir=str(tmp_path)
        )
        plt.close("all")
