"""
Thermodynamic properties for ethanol-water system.
Includes VLE calculations, Antoine equation, and property estimation.
"""

import numpy as np
from config import (
    ETHANOL_ANTOINE_A,
    ETHANOL_ANTOINE_B,
    ETHANOL_ANTOINE_C,
    WATER_ANTOINE_A,
    WATER_ANTOINE_B,
    WATER_ANTOINE_C,
    ETHANOL_MW,
    WATER_MW,
    ETHANOL_DENSITY,
    WATER_DENSITY,
    ETHANOL_LATENT_HEAT,
    WATER_LATENT_HEAT,
    ETHANOL_CP,
    WATER_CP,
    ETHANOL_BP,
    WATER_BP,
)


def antoine_pressure(T, A, B, C):
    """
    Calculate vapor pressure using Antoine equation.

    Parameters:
        T: Temperature (°C)
        A, B, C: Antoine constants

    Returns:
        P: Vapor pressure (mmHg)
    """
    return 10 ** (A - B / (C + T))


def ethanol_vapor_pressure(T):
    """Vapor pressure of pure ethanol (mmHg)."""
    return antoine_pressure(T, ETHANOL_ANTOINE_A, ETHANOL_ANTOINE_B, ETHANOL_ANTOINE_C)


def water_vapor_pressure(T):
    """Vapor pressure of pure water (mmHg)."""
    return antoine_pressure(T, WATER_ANTOINE_A, WATER_ANTOINE_B, WATER_ANTOINE_C)


def mmhg_to_bar(p_mmhg):
    """Convert mmHg to bar."""
    return p_mmhg / 750.062


def bar_to_mmhg(p_bar):
    """Convert bar to mmHg."""
    return p_bar * 750.062


def vol_to_mol_fraction(x_vol):
    """
    Convert volume fraction to mole fraction for ethanol-water mixture.

    Assumes ideal mixing (volumes are additive).

    Parameters:
        x_vol: Volume fraction of ethanol

    Returns:
        x_mol: Mole fraction of ethanol
    """
    # Mass fractions from volume fractions (using pure component densities)
    # m_eth = V_eth * rho_eth, m_water = V_water * rho_water
    rho_eth = ETHANOL_DENSITY  # kg/m³
    rho_water = WATER_DENSITY  # kg/m³

    # For 1 L total: V_eth = x_vol, V_water = 1 - x_vol
    m_eth = x_vol * rho_eth
    m_water = (1 - x_vol) * rho_water

    # Moles
    n_eth = m_eth / ETHANOL_MW * 1000  # mol (MW in g/mol, m in kg)
    n_water = m_water / WATER_MW * 1000

    x_mol = n_eth / (n_eth + n_water)
    return x_mol


def mol_to_vol_fraction(x_mol):
    """
    Convert mole fraction to volume fraction for ethanol-water mixture.

    Parameters:
        x_mol: Mole fraction of ethanol

    Returns:
        x_vol: Volume fraction of ethanol
    """
    # For 1 mol total mixture
    n_eth = x_mol
    n_water = 1 - x_mol

    # Mass (g)
    m_eth = n_eth * ETHANOL_MW
    m_water = n_water * WATER_MW

    # Volume (mL, using density in kg/m³ = g/L)
    v_eth = m_eth / ETHANOL_DENSITY * 1000  # mL
    v_water = m_water / WATER_DENSITY * 1000  # mL

    x_vol = v_eth / (v_eth + v_water)
    return x_vol


def activity_coefficients_van_laar(x_eth, T=None):
    """
    Calculate activity coefficients using Van Laar equation for ethanol-water.

    Van Laar parameters for ethanol(1)-water(2) at ~1 atm:
    A12 = 1.6798, A21 = 0.9227 (from experiment)

    Parameters:
        x_eth: Mole fraction of ethanol
        T: Temperature (°C) - not used in simple Van Laar

    Returns:
        gamma_eth, gamma_water: Activity coefficients
    """
    # Van Laar parameters for ethanol-water system
    A12 = 1.6798  # ethanol parameter
    A21 = 0.9227  # water parameter

    x1 = x_eth
    x2 = 1 - x_eth

    # Avoid division by zero
    if x1 < 1e-10:
        return 1.0, np.exp(A21)
    if x2 < 1e-10:
        return np.exp(A12), 1.0

    # Van Laar equations
    denom1 = (A12 * x1 / (A12 * x1 + A21 * x2)) ** 2
    denom2 = (A21 * x2 / (A12 * x1 + A21 * x2)) ** 2

    ln_gamma1 = A12 * (A21 * x2 / (A12 * x1 + A21 * x2)) ** 2
    ln_gamma2 = A21 * (A12 * x1 / (A12 * x1 + A21 * x2)) ** 2

    gamma_eth = np.exp(ln_gamma1)
    gamma_water = np.exp(ln_gamma2)

    return gamma_eth, gamma_water


def bubble_point_temperature(x_eth, P_total_bar):
    """
    Calculate bubble point temperature for ethanol-water mixture.

    Uses modified Raoult's law with activity coefficients:
    P_total = x1 * gamma1 * P1_sat + x2 * gamma2 * P2_sat

    Parameters:
        x_eth: Mole fraction of ethanol in liquid
        P_total_bar: Total pressure (bar)

    Returns:
        T_bubble: Bubble point temperature (°C)
    """
    P_total = bar_to_mmhg(P_total_bar)

    # Handle pure components
    if x_eth <= 1e-10:
        # Pure water - solve P_water_sat(T) = P_total
        T_guess = WATER_BP
        for _ in range(50):
            P_calc = water_vapor_pressure(T_guess)
            if abs(P_calc - P_total) < 0.1:
                break
            dT = 0.01
            dP_dT = (water_vapor_pressure(T_guess + dT) - P_calc) / dT
            T_guess = T_guess - (P_calc - P_total) / dP_dT
        return T_guess

    if x_eth >= 1.0 - 1e-10:
        # Pure ethanol - solve P_eth_sat(T) = P_total
        T_guess = ETHANOL_BP
        for _ in range(50):
            P_calc = ethanol_vapor_pressure(T_guess)
            if abs(P_calc - P_total) < 0.1:
                break
            dT = 0.01
            dP_dT = (ethanol_vapor_pressure(T_guess + dT) - P_calc) / dT
            T_guess = T_guess - (P_calc - P_total) / dP_dT
        return T_guess

    # Mixture: initial guess from linear interpolation
    T_guess = x_eth * ETHANOL_BP + (1 - x_eth) * WATER_BP

    # Newton-Raphson iteration
    for _ in range(50):
        gamma_eth, gamma_water = activity_coefficients_van_laar(x_eth, T_guess)

        P_eth = ethanol_vapor_pressure(T_guess)
        P_water = water_vapor_pressure(T_guess)

        P_calc = x_eth * gamma_eth * P_eth + (1 - x_eth) * gamma_water * P_water

        # Convergence check
        if abs(P_calc - P_total) < 0.1:  # 0.1 mmHg tolerance
            break

        # Numerical derivative
        dT = 0.01
        P_eth_plus = ethanol_vapor_pressure(T_guess + dT)
        P_water_plus = water_vapor_pressure(T_guess + dT)
        P_calc_plus = (
            x_eth * gamma_eth * P_eth_plus + (1 - x_eth) * gamma_water * P_water_plus
        )

        dP_dT = (P_calc_plus - P_calc) / dT

        if abs(dP_dT) < 1e-10:
            break

        T_guess = T_guess - (P_calc - P_total) / dP_dT

        # Keep in reasonable range
        T_guess = max(70, min(110, T_guess))

    return T_guess


def vapor_composition(x_eth, T, P_total_bar):
    """
    Calculate vapor composition in equilibrium with liquid.

    Parameters:
        x_eth: Mole fraction of ethanol in liquid
        T: Temperature (°C)
        P_total_bar: Total pressure (bar)

    Returns:
        y_eth: Mole fraction of ethanol in vapor
    """
    # Handle pure components
    if x_eth <= 1e-10:
        return 0.0
    if x_eth >= 1.0 - 1e-10:
        return 1.0

    P_total = bar_to_mmhg(P_total_bar)

    gamma_eth, gamma_water = activity_coefficients_van_laar(x_eth, T)

    P_eth = ethanol_vapor_pressure(T)

    # Partial pressure of ethanol
    p_eth = x_eth * gamma_eth * P_eth

    y_eth = p_eth / P_total

    # Ensure physical bounds
    return max(0.0, min(1.0, y_eth))


def equilibrium_curve(P_bar, n_points=100):
    """
    Generate VLE equilibrium data (x-y curve) at given pressure.

    Parameters:
        P_bar: Total pressure (bar)
        n_points: Number of data points

    Returns:
        x_array: Liquid mole fractions
        y_array: Vapor mole fractions
        T_array: Bubble point temperatures
    """
    x_array = np.linspace(0, 1, n_points)
    y_array = np.zeros(n_points)
    T_array = np.zeros(n_points)

    for i, x in enumerate(x_array):
        T = bubble_point_temperature(x, P_bar)
        y = vapor_composition(x, T, P_bar)

        T_array[i] = T
        y_array[i] = min(y, 1.0)  # Cap at 1.0

    return x_array, y_array, T_array


def mixture_latent_heat(x_eth):
    """
    Calculate mixture latent heat of vaporization (J/mol).

    Simple molar average (approximation).

    Parameters:
        x_eth: Mole fraction of ethanol

    Returns:
        lambda_mix: Mixture latent heat (J/mol)
    """
    return x_eth * ETHANOL_LATENT_HEAT + (1 - x_eth) * WATER_LATENT_HEAT


def mixture_molecular_weight(x_eth):
    """
    Calculate mixture molecular weight (g/mol).

    Parameters:
        x_eth: Mole fraction of ethanol

    Returns:
        MW_mix: Mixture molecular weight (g/mol)
    """
    return x_eth * ETHANOL_MW + (1 - x_eth) * WATER_MW


def liquid_heat_capacity(x_eth):
    """
    Calculate mixture liquid heat capacity (J/(kg·K)).

    Mass-weighted average.

    Parameters:
        x_eth: Mole fraction of ethanol

    Returns:
        Cp_mix: Mixture heat capacity (J/(kg·K))
    """
    # Convert to mass fraction
    m_eth = x_eth * ETHANOL_MW
    m_water = (1 - x_eth) * WATER_MW
    w_eth = m_eth / (m_eth + m_water)
    w_water = 1 - w_eth

    return w_eth * ETHANOL_CP + w_water * WATER_CP


def liquid_density(x_vol_eth):
    """
    Calculate mixture liquid density (kg/m³).

    Simple volume-weighted average (approximation).

    Parameters:
        x_vol_eth: Volume fraction of ethanol

    Returns:
        rho_mix: Mixture density (kg/m³)
    """
    # More accurate would use excess volume, but this is reasonable approximation
    return x_vol_eth * ETHANOL_DENSITY + (1 - x_vol_eth) * WATER_DENSITY


if __name__ == "__main__":
    # Test calculations
    print("=" * 60)
    print("Thermodynamics Module Test")
    print("=" * 60)

    # Test volume to mole fraction conversion
    print("\nVolume to Mole Fraction Conversion:")
    for x_vol in [0.10, 0.50, 0.90, 0.956]:
        x_mol = vol_to_mol_fraction(x_vol)
        x_vol_back = mol_to_vol_fraction(x_mol)
        print(
            f"  {x_vol * 100:.1f}% vol -> {x_mol * 100:.2f}% mol -> {x_vol_back * 100:.2f}% vol"
        )

    # Test bubble point calculation
    print("\nBubble Point Temperatures at 1.0 bar:")
    P = 1.0  # bar
    for x_vol in [0.0, 0.10, 0.50, 0.90, 0.956, 1.0]:
        x_mol = vol_to_mol_fraction(x_vol) if x_vol > 0 and x_vol < 1 else x_vol
        T_bp = bubble_point_temperature(x_mol, P)
        y_mol = vapor_composition(x_mol, T_bp, P)
        print(
            f"  x = {x_mol * 100:.1f}% mol: T_bp = {T_bp:.2f}°C, y = {y_mol * 100:.1f}% mol"
        )

    # Test equilibrium curve
    print("\nGenerating equilibrium curve at 1.0 bar...")
    x, y, T = equilibrium_curve(1.0, n_points=11)
    print("  x (mol%)  |  y (mol%)  |  T (°C)")
    print("  " + "-" * 35)
    for i in range(len(x)):
        print(f"  {x[i] * 100:6.1f}    |  {y[i] * 100:6.1f}    |  {T[i]:6.2f}")

    print("\nLatent heats:")
    print(f"  Pure ethanol: {ETHANOL_LATENT_HEAT / 1000:.2f} kJ/mol")
    print(f"  Pure water: {WATER_LATENT_HEAT / 1000:.2f} kJ/mol")
    print(f"  90% ethanol mix: {mixture_latent_heat(0.775) / 1000:.2f} kJ/mol")
