"""Integration tests for ThermophysicalProperties using the sample case files."""

import os

import numpy as np
import pytest

# Skip the entire module if PyFoam is not importable (broken install, etc.)
pytest.importorskip("PyFoam")

from FoamUtils.ThermophysicalProperties import ThermophysicalProperties

# Path to the sample thermophysicalProperties file bundled with the repo.
CASE_DIR = os.path.join(os.path.dirname(__file__), "..", "constant")
THERMO_FILE = os.path.join(CASE_DIR, "thermophysicalProperties")


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
