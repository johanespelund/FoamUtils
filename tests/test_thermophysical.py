"""Tests for ThermophysicalProperties pure-Python functions."""

import numpy as np
import pytest

from FoamUtils.ThermophysicalProperties import (
    rho_perfectGas,
    rho_PengRobinson,
    rPolynomial,
    sutherland_mu,
    sutherland_kappa,
)


# ── rho_perfectGas ────────────────────────────────────────────────────────────

class TestRhoPerfectGas:
    """Tests for the ideal-gas density function ρ = p / (R T)."""

    def test_scalar_300K(self):
        """Air at 300 K and 101 325 Pa with R=287 J/(kg·K)."""
        rho = rho_perfectGas(T=300.0, p=101325.0, R=287.0)
        expected = 101325.0 / (287.0 * 300.0)
        assert rho == pytest.approx(expected, rel=1e-9)

    def test_array_input(self):
        """Vectorised call returns an array of the correct shape."""
        T = np.array([200.0, 300.0, 400.0])
        rho = rho_perfectGas(T=T, p=101325.0, R=287.0)
        expected = 101325.0 / (287.0 * T)
        np.testing.assert_allclose(rho, expected, rtol=1e-9)

    def test_doubles_when_temperature_halved(self):
        """Halving T should double ρ."""
        rho1 = rho_perfectGas(T=300.0, p=1e5, R=287.0)
        rho2 = rho_perfectGas(T=150.0, p=1e5, R=287.0)
        assert rho2 == pytest.approx(2 * rho1, rel=1e-9)

    def test_proportional_to_pressure(self):
        """Doubling p should double ρ."""
        rho1 = rho_perfectGas(T=300.0, p=1e5, R=287.0)
        rho2 = rho_perfectGas(T=300.0, p=2e5, R=287.0)
        assert rho2 == pytest.approx(2 * rho1, rel=1e-9)


# ── rPolynomial ───────────────────────────────────────────────────────────────

class TestRPolynomial:
    """Tests for the rPolynomial equation of state: 1/ρ = C0 + C1 T + C2 T² - C3 p - C4 p T."""

    # Coefficients from the repository's 'parameters' file (LH2 at low T).
    COEFFS = [
        -0.12681525489315768,
         0.08138876620175439,
        -1.177435017248281e-05,
        -2.688141687658478e-07,
         3.795559337408709e-07,
    ]

    def test_returns_positive_density(self):
        rho = rPolynomial(T=22.0, p=1e5, coeffs=self.COEFFS)
        assert rho > 0

    def test_constant_coefficients(self):
        """With only C0 set, ρ = 1/C0 (independent of T and p)."""
        coeffs = [0.001, 0.0, 0.0, 0.0, 0.0]
        rho = rPolynomial(T=300.0, p=1e5, coeffs=coeffs)
        assert rho == pytest.approx(1000.0, rel=1e-9)

    def test_array_input(self):
        """Vectorised call returns an array of correct size."""
        T = np.linspace(20.0, 30.0, 5)
        rho = rPolynomial(T=T, p=1e5, coeffs=self.COEFFS)
        assert rho.shape == (5,)
        assert np.all(rho > 0)


# ── rho_PengRobinson ──────────────────────────────────────────────────────────

class TestRhoPengRobinson:
    """Tests for the Peng-Robinson equation of state."""

    # LH2 critical parameters from the repository's 'parameters' file.
    R = 8.31446261815324 / (2.02e-3)   # J/(kg·K)  for H2
    Tc = 33.145   # K
    Pc = 1.2964e6  # Pa
    omega = -0.219

    def test_scalar_returns_float(self):
        rho = rho_PengRobinson(T=25.0, p=1e5, R=self.R, Tc=self.Tc, Pc=self.Pc, omega=self.omega)
        assert isinstance(rho, float)
        assert rho > 0

    def test_array_returns_array(self):
        T = np.linspace(20.0, 30.0, 4)
        rho = rho_PengRobinson(T=T, p=1e5, R=self.R, Tc=self.Tc, Pc=self.Pc, omega=self.omega)
        assert rho.shape == (4,)
        assert np.all(rho > 0)

    def test_density_increases_with_pressure(self):
        """Higher pressure → higher density (supercritical gas-like behaviour)."""
        rho_low = rho_PengRobinson(T=40.0, p=1e5, R=self.R, Tc=self.Tc, Pc=self.Pc, omega=self.omega)
        rho_high = rho_PengRobinson(T=40.0, p=2e5, R=self.R, Tc=self.Tc, Pc=self.Pc, omega=self.omega)
        assert rho_high > rho_low

    def test_approaches_ideal_gas_at_low_pressure(self):
        """At low pressure, Peng-Robinson should agree with ideal gas within 1 %."""
        T = 100.0  # well above Tc
        p = 1000.0  # very low pressure
        rho_pr = rho_PengRobinson(T=T, p=p, R=self.R, Tc=self.Tc, Pc=self.Pc, omega=self.omega)
        rho_ideal = p / (self.R * T)
        assert rho_pr == pytest.approx(rho_ideal, rel=0.01)


# ── sutherland_mu ─────────────────────────────────────────────────────────────

class TestSutherlandMu:
    """Tests for Sutherland's viscosity law: μ = As √T / (1 + Ts/T)."""

    # Air constants
    As = 1.458e-6  # kg/(m·s·K^0.5)
    Ts = 110.4     # K

    def test_positive_viscosity(self):
        mu = sutherland_mu(self.As, self.Ts, T=300.0)
        assert mu > 0

    def test_known_value_air_300K(self):
        """μ_air ≈ 1.846e-5 Pa·s at 300 K (within 1 %)."""
        mu = sutherland_mu(self.As, self.Ts, T=300.0)
        assert mu == pytest.approx(1.846e-5, rel=0.01)

    def test_viscosity_increases_with_temperature(self):
        mu1 = sutherland_mu(self.As, self.Ts, T=200.0)
        mu2 = sutherland_mu(self.As, self.Ts, T=400.0)
        assert mu2 > mu1

    def test_array_input(self):
        T = np.array([200.0, 300.0, 400.0])
        mu = sutherland_mu(self.As, self.Ts, T)
        assert mu.shape == (3,)
        assert np.all(mu > 0)
        assert np.all(np.diff(mu) > 0)  # monotonically increasing


# ── sutherland_kappa ──────────────────────────────────────────────────────────

class TestSutherlandKappa:
    """Tests for Sutherland thermal conductivity."""

    As = 1.458e-6
    Ts = 110.4
    R = 287.0   # J/(kg·K) for air

    def test_positive_conductivity(self):
        Cp = 1005.0
        kappa = sutherland_kappa(self.As, self.Ts, T=300.0, Cp=Cp, R=self.R)
        assert kappa > 0

    def test_conductivity_increases_with_temperature(self):
        Cp = 1005.0
        k1 = sutherland_kappa(self.As, self.Ts, T=200.0, Cp=Cp, R=self.R)
        k2 = sutherland_kappa(self.As, self.Ts, T=400.0, Cp=Cp, R=self.R)
        assert k2 > k1

    def test_array_input(self):
        T = np.array([200.0, 300.0, 400.0])
        Cp = 1005.0 * np.ones_like(T)
        kappa = sutherland_kappa(self.As, self.Ts, T, Cp, self.R)
        assert kappa.shape == (3,)
        assert np.all(kappa > 0)


# ── PyFoam import-error handling ──────────────────────────────────────────────

class TestPyFoamImportErrorHandling:
    """Verify that a broken PyFoam install raises a clear, actionable error."""

    def test_broken_pyfoam_raises_helpful_importerror(self, tmp_path, monkeypatch):
        """ThermophysicalProperties.__init__ converts a broken-PyFoam ImportError."""
        import builtins
        import sys
        from FoamUtils.ThermophysicalProperties import ThermophysicalProperties

        real_import = builtins.__import__

        def mock_import(name, *args, **kwargs):
            if name.startswith("PyFoam"):
                raise ImportError("No module named 'PyFoam.ThirdParty.six.moves'")
            return real_import(name, *args, **kwargs)

        # Create a dummy (empty) file so the FileNotFoundError is not raised first.
        dummy = tmp_path / "thermophysicalProperties"
        dummy.write_text("")

        # Remove any cached PyFoam modules so our mock takes effect.
        cached = [k for k in sys.modules if k.startswith("PyFoam")]
        for k in cached:
            monkeypatch.delitem(sys.modules, k)

        monkeypatch.setattr(builtins, "__import__", mock_import)

        with pytest.raises(ImportError, match="PyFoam"):
            ThermophysicalProperties(str(dummy))
