import os

import click
import matplotlib.pyplot as plt
import numpy as np
import scienceplots
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

plt.style.use(["science", "no-latex", "grid"])


class ThermophysicalProperties:
    def __init__(self, file_path):
        """
        Create a ThermophysicalProperties class by loading and parsing
        the thermophysicalProperties file using PyFoam.

        Args:
            file_path (str): Path to the thermophysicalProperties file.
        """
        self.file_path = file_path
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")
        self.properties = ParsedParameterFile(file_path, doMacroExpansion=True)
        self.M = self.properties["mixture"]["specie"]["molWeight"] * 1e-3
        self.R = 8.31446261815324 / self.M
        self.create_functions()

    def create_functions(self):

        ### thermo ###
        if self.properties["thermoType"]["thermo"] == "hPolynomial":
            coeffs = np.flip(
                self.properties["mixture"]["thermodynamics"]["CpCoeffs<8>"]
            )
            self.Cp = np.poly1d(coeffs)
        elif self.properties["thermoType"]["thermo"] == "hConst":
            Cp = self.properties["mixture"]["thermodynamics"]["Cp"]
            self.Cp = lambda T: np.ones_like(T) * Cp

        ### transport ###
        if self.properties["thermoType"]["transport"] == "polynomial":
            coeffs = np.flip(self.properties["mixture"]["transport"]["muCoeffs<8>"])
            self.mu = np.poly1d(coeffs)
            coeffs = np.flip(self.properties["mixture"]["transport"]["kappaCoeffs<8>"])
            self.kappa = np.poly1d(coeffs)
            self.Pr = lambda T: self.mu(T) * self.Cp(T) / self.kappa(T)
        elif self.properties["thermoType"]["transport"] == "const":
            mu = self.properties["mixture"]["transport"]["mu"]
            Pr = self.properties["mixture"]["transport"]["Pr"]
            self.mu = lambda T: mu * np.ones_like(T)
            self.Pr = lambda T: Pr * np.ones_like(T)
            self.kappa = lambda T: mu * self.Cp(T) / Pr
        elif self.properties["thermoType"]["transport"] == "sutherland":
            Ts = self.properties["mixture"]["transport"]["Ts"]
            As = self.properties["mixture"]["transport"]["As"]
            Pr = self.properties["mixture"]["transport"]["Pr"]
            self.Pr = lambda T: Pr * np.ones_like(T)
            self.mu = lambda T: sutherland_mu(As, Ts, T)
            self.kappa = lambda T: sutherland_kappa(As, Ts, T, self.Cp(T), self.R)

        ### equationOfState ###
        if self.properties["thermoType"]["equationOfState"] == "rPolynomial":
            coeffs = np.flip(self.properties["mixture"]["equationOfState"]["C"])
            self.rho = lambda p, T: rPolynomial(T, p, coeffs)
        elif self.properties["thermoType"]["equationOfState"] == "perfectGas":
            self.rho = lambda p, T: rho_perfectGas(T, p, self.R)
        elif self.properties["thermoType"]["equationOfState"] == "PengRobinsonGas":
            Tc = self.properties["mixture"]["equationOfState"]["Tc"]
            Pc = self.properties["mixture"]["equationOfState"]["Pc"]
            omega = self.properties["mixture"]["equationOfState"]["omega"]
            self.rho = lambda p, T: rho_PengRobinson(T, p, self.R, Tc, Pc, omega)
        elif self.properties["thermoType"]["equationOfState"] == "icoPolynomial":
            coeffs = np.flip(self.properties["mixture"]["equationOfState"]["rhoCoeffs<8>"])
            self.rho = lambda p, T: np.poly1d(coeffs)(T)

    def beta(self, p, T):
        """
        Calculate the thermal expansion coefficient.
        """

        if self.properties["thermoType"]["equationOfState"] == "perfectGas":
            return 1 / T
        elif self.properties["thermoType"]["equationOfState"] == "PengRobinsonGas":
            # Use finite difference to calculate the derivative d(rho)/dT
            T1 = T + 1e-3
            rho1 = self.rho(p, T1)
            rho = self.rho(p, T)
            drhodT = (rho1 - rho) / (T1 - T)
            return -1 / self.rho(p, T) * drhodT

    def alpha(self, p, T):
        """
        Thermal diffusivity.
        """
        return self.kappa(T) / (self.rho(p, T) * self.Cp(T))

    def plot(self, T, p, properties, axes=None, color="C0", legend=False):
        """
        Plot the thermophysical properties.
        T (np.array): Temperature
        p (float): Pressure
        properties (list): List of properties to plot
        """
        n = len(properties)
        if axes is None:
            showplot = True
            N, M = int(np.ceil(np.sqrt(n))), int(np.ceil(n / np.ceil(np.sqrt(n))))
            fig, ax = plt.subplots(N, M, figsize=(N * 3, M * 2))
            # Set x-axis label for the last row
            for a in ax[-1]:
                a.set_xlabel("Temperature (K)")
            ax = ax.flatten()
        else:
            showplot = False
            ax = axes

        for i, prop in enumerate(properties):
            if prop == "rho":
                ax[i].plot(T, self.rho(p, T), label="Density", color=color)
            elif prop == "Cp":
                ax[i].plot(T, self.Cp(T), label="Heat capacity", color=color)
            elif prop == "mu":
                ax[i].plot(T, self.mu(T)*1e6, label="Dynamic viscosity", color=color)
            elif prop == "kappa":
                ax[i].plot(T, self.kappa(T), label="Thermal conductivity", color=color)
            elif prop == "beta":
                ax[i].plot(
                    T,
                    self.beta(p, T),
                    label="Thermal expansion coefficient",
                    color=color,
                )
            ax[i].set_ylabel(self.unit_label(prop))
            if legend:
                ax[i].legend()

        if showplot:
            fig.tight_layout()
            plt.show()

    def unit_label(self, prop):
        """
        Return the unit of a property.
        """
        if prop == "rho":
            return "$\\rho$ (kg/m³)"
        elif prop == "Cp":
            return "$C_p$ (J/kg $\cdot$ K)"
        elif prop == "mu":
            return "$\\mu$ (Pa $\cdot$ s)"
        elif prop == "kappa":
            return "$\\kappa$ (W/m $\cdot$ K)"
        elif prop == "beta":
            return "$\\beta$ (1/K)"

    def __repr__(self):
        return f"ThermophysicalProperties({self.file_path})"


def rho_perfectGas(T, p, R):
    """
    Calculate the density of a gas using the Perfect Gas law.

    Parameters:
    T (float or array-like): Temperature(s) in Kelvin.
    p (float): Pressure in Pascals.
    R (float): Specific gas constant in J/(kg·K).

    Returns:
    float or ndarray: Density in kg/m³.
    """
    return p / (R * T)


def rho_PengRobinson(T, p, R, Tc, Pc, omega):
    """
    Calculate the density using the Peng-Robinson equation of state.

    Parameters:
    T (float or array-like): Temperature(s) in Kelvin.
    p (float): Pressure in Pascals.
    R (float): Universal gas constant in J/(mol·K).
    Tc (float): Critical temperature in Kelvin.
    Pc (float): Critical pressure in Pascals.
    omega (float): Acentric factor.

    Returns:
    float or ndarray: Density in mol/m³.
    """
    T = np.atleast_1d(T)
    Tr = T / Tc

    # Calculate parameters for the Peng-Robinson EoS
    a = 0.45724 * (R * Tc) ** 2 / Pc
    b = 0.07780 * R * Tc / Pc
    kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
    alpha = (1 + kappa * (1 - np.sqrt(Tr))) ** 2
    A = alpha * a * p / (R * T) ** 2
    B = b * p / (R * T)

    Z = np.zeros_like(T)

    for i in range(len(T)):
        coeffs = [
            1,
            -(1 - B[i]),
            A[i] - 2 * B[i] - 3 * B[i] ** 2,
            -(A[i] * B[i] - B[i] ** 2 - B[i] ** 3),
        ]
        Z_roots = np.roots(coeffs)
        Z[i] = np.real(Z_roots[np.isreal(Z_roots)])[0]  # Select the real root

    rho = p / (Z * R * T)

    return rho.item() if rho.size == 1 else rho


def rPolynomial(T, p, coeffs):
    """
    Calculate the density using a polynomial equation of state.

    Parameters:
    T (float or array-like): Temperature(s) in Kelvin.
    p (float): Pressure in Pascals.
    coeffs (list): Coefficients of the polynomial equation of state.
                   1/rho = C0 + C1*T + C2*T**2 - C3*p - C4*p*T

    Returns:
    float or ndarray: Density in kg/m³.
    """
    c = coeffs
    return 1 / (c[0] + c[1] * T + c[2] * T**2 - c[3] * p - c[4] * p * T)


def sutherland_mu(As, Ts, T):
    return As * np.sqrt(T) / (1 + Ts / T)


def sutherland_kappa(As, Ts, T, Cp, R):
    Cv = Cp - R
    mu = sutherland_mu(As, Ts, T)
    # return mu*1.4*R/(0.4*0.71)
    return mu * Cv * (1.32 + 1.77 * R / Cv)


def plot_thermophysical_properties(
    p, Tstart, Tend, properties=["rho", "Cp", "mu", "kappa"], legend=False, case_dir="."
):
    """
    Plot the thermophysical properties of the OpenFOAM case.
    """
    T = np.linspace(Tstart, Tend)
    thermo = ThermophysicalProperties(f"{case_dir}/constant/thermophysicalProperties")
    thermo.plot(T, p, properties, legend=legend)


@click.command()
@click.option("--p", type=float, default=1e5, help="Pressure in Pascals.")
@click.option("--Tstart", type=float, default=20.3, help="Start temperature in Kelvin.")
@click.option("--Tend", type=float, default=100, help="End temperature in Kelvin.")
@click.option(
    "--properties", type=str, default="rho Cp mu kappa", help="Properties to plot."
)
@click.option("--legend", is_flag=True, help="Show legend.")
@click.option(
    "--case_dir", type=str, default=".", help="Path to the OpenFOAM case directory."
)
def main(p, tstart, tend, properties, legend, case_dir):
    properties = properties.split()
    plot_thermophysical_properties(p, tstart, tend, properties, legend, case_dir)


if __name__ == "__main__":
    main()
