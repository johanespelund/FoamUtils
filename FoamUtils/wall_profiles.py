#! /usr/bin/env python3 

import numpy as np
import matplotlib.pyplot as plt
from . import ThermophysicalProperties, file_utils

## Utilities for post processing data to get wall profiles for U+ and T+



def generate_polynomial(foam_string):
    """
    Generate a numpy polynomial from a string in OpenFOAM format, e.g 
    (  2.3157e+04 -1.1462e+03  4.2657e+01 -8.2923e-01  8.8072e-03 -4.8108e-05  1.0572e-07 0 )
    Here, the first element is the constant, the second is the linear term, etc.
    """
    coeffs = foam_string.replace("(", "").replace(")", "").split()
    coeffs = [float(c) for c in reversed(coeffs)]
    return np.poly1d(coeffs)


def Cp(T):
    """
    Calculate the specific heat capacity of air as a function of temperature
    """
    poly = generate_polynomial(p["NIST_Cp"])
    return poly(T)


def mu(T):
    """
    Calculate the dynamic viscosity of air as a function of temperature
    """
    poly = generate_polynomial(p["NIST_Viscosity"])
    return poly(T)

def kappa(T):
    """
    Calculate the thermal conductivity of air as a function of temperature
    """
    poly = generate_polynomial(p["NIST_ThermalConductivity"])
    return poly(T)


def Pr(T):
    """
    Calculate the Prandtl number of air as a function of temperature
    """
    return mu(T) * Cp(T) / kappa(T)

def average(T_low, T_high, property):
    """
    Calculate the average of a property over a temperature range
    T_low (float): The lower bound of the temperature range
    T_high (float): The upper bound of the temperature range
    property (function): The property to average
    """
    T = np.linspace(T_low, T_high)
    return np.trapz(property(T), T) / (T_high - T_low)


def tau_wall(u, n, T_w, thermo):
    """
    Calculate the wall shear stress given the velocity profile and wall temperature
    u (np.array): The velocity profile
    n (np.array): The normal distance profile
    T_w (float): The wall temperature
    """
    dudn = np.gradient(u, n)[0]
    mu_wall = thermo.mu(T_w)
    return mu_wall * dudn

def q_wall(T, n, thermo):
    """
    Calculate the heat flux at the wall given the temperature profile and normal distance
    T (np.array): The temperature profile
    n (np.array): The normal distance profile
    """
    dTdn = np.gradient(T, n)[0]
    kappa_wall = thermo.kappa(T[0])
    return -kappa_wall * dTdn

def friction_velocity(tau_wall, rho):
    """
    Calculate the friction velocity given the wall shear stress and density
    tau_wall (np.array): The wall shear stress
    rho (float): The density
    """
    return np.sqrt(abs(tau_wall / rho))


def friction_temperature(q_wall, tau_wall, T, rho, thermo):
    """
    Calculate the friction temperature given the heat flux, density, and specific heat capacity
    q_wall (np.array): The heat flux at the wall
    tau_wall (np.array): The wall shear stress
    rho (float): The density
    """
    return q_wall / (rho * thermo.Cp(T[0]) * friction_velocity(tau_wall, rho))

def u_plus(u, T, n, rho, thermo):
    """
    Calculate the dimensionless velocity profile
    u (np.array): The velocity profile
    T (np.array): The temperature profile
    n (np.array): The normal distance profile
    rho (float): The density
    """
    tau_w = tau_wall(u, n, T[0], thermo)
    return u / friction_velocity(tau_w, rho)

def phi_plus(u, T, n, rho, Tref, thermo):
    """
    Calculate the dimensionless temperature profile
    T (np.array): The temperature profile
    n (np.array): The normal distance profile
    rho (float): The density
    mu (np.array): The dynamic viscosity
    """
    tau_w = tau_wall(u, n, T[0], thermo)
    q_w = q_wall(T, n, thermo)
    Cp_w = thermo.Cp(T[0])
    phi_tau = friction_temperature(q_w, tau_w, T, rho, thermo)
    phi_w = T[0] - Tref
    # phi_bar = np.trapz(T - Tref, n) / (n[-1] - n[0])
    phi_bar = T - Tref
    return (phi_w - phi_bar) / (phi_tau)


def n_plus(u, T, n, rho, thermo):
    """
    Calculate the dimensionless distance profile
    u (np.array): The velocity profile
    T (np.array): The temperature profile
    n (np.array): The normal distance profile
    rho (float): The density
    """
    u_tau = friction_velocity(tau_wall(u, n, T[0], thermo), rho)
    return n * rho * u_tau / thermo.mu(T[0])


def Nusselt(q, deltaT, L, kappa):
    """
    Calculates the local Nusselt number.
    """
    return q * L / (kappa * deltaT)


def get_Ra(p):
    """
    Calculate the Rayleigh number from the parameters
    """
    g = 9.81
    T_c = p["T_left"]
    T_h = p["T_right"]
    T = (T_c + T_h) / 2
    beta = beta_PengRobinson(T, p["p_outlet"], p)
    L = p["L_x"]
    rho = rho_PengRobinson(T, p["p_outlet"], p)
    print(f"rho: {rho}, beta: {beta}, T: {T}")
    nu = mu(T) / rho
    return g * beta * (T_h - T_c) * L ** 3 / nu ** 2

def rho_PengRobinson(T, P, params):
    """
    Calculate the density of air using the Peng-Robinson equation of state
    T (np.array): The temperature profile
    p (dict): The parameters
    """
    M = params["M"]*1e-3
    R = 8.31446261815324/M
    Tc = params["Tc"]
    Pc = params["Pc"]
    omega = params["omega"]
    Tr = T / Tc

    # Calculate parameters for the Peng-Robinson EoS
    a = 0.45724 * (R * Tc)**2 / Pc
    b = 0.07780 * R * Tc / Pc
    kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
    alpha = (1 + kappa * (1 - np.sqrt(Tr)))**2
    A = alpha * a * P / (R * T)**2
    B = b * P / (R * T)

    # Cubic equation coefficients: Z^3 - (1 - B)Z^2 + (A - 2B - 3B^2)Z - (AB - B^2 - B^3) = 0
    coeffs = [1, -(1 - B), A - 2 * B - 3 * B**2, -(A * B - B**2 - B**3)]

    # Solve the cubic equation for Z
    Z_roots = np.roots(coeffs)
    Z = np.real(Z_roots[np.isreal(Z_roots)])[0]  # Select the real root

    return P / (Z * R * T)

def beta_PengRobinson(T, P, p):
    """
    Calculate the thermal expansion coefficient of air using the Peng-Robinson equation of state.
    """
    T1 = T + 1e-3
    rho1 = rho_PengRobinson(T1, P, p)
    rho = rho_PengRobinson(T, P, p)
    drhodT = (rho1 - rho) / (T1 - T)
    return -1 / rho_PengRobinson(T, P, p) * drhodT



