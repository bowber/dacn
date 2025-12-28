"""
McCabe-Thiele analysis for ethanol-water distillation.
Calculates minimum reflux ratio and operating reflux ratio.
"""

import numpy as np
from scipy.optimize import brentq
from thermodynamics import (
    vol_to_mol_fraction,
    mol_to_vol_fraction,
    bubble_point_temperature,
    vapor_composition,
    equilibrium_curve,
    mixture_latent_heat,
    liquid_heat_capacity,
    mixture_molecular_weight,
)
from config import (
    FEED_CONC_VOL,
    PRODUCT_CONC_VOL,
    BOTTOM_CONC_VOL,
    FEED_TEMP,
    OPERATING_PRESSURE,
    FEED_FLOW_RATE,
    NUM_TRAYS,
    TRAY_EFFICIENCY,
    FEED_TRAY,
)


def calculate_q_factor(x_F, T_feed, P_bar):
    """
    Calculate the q-factor for feed condition.

    q = (heat needed to vaporize 1 mol of feed) / (molar latent heat)

    For subcooled liquid: q > 1
    For saturated liquid: q = 1
    For partially vaporized: 0 < q < 1
    For saturated vapor: q = 0
    For superheated vapor: q < 0

    Parameters:
        x_F: Feed mole fraction of ethanol
        T_feed: Feed temperature (°C)
        P_bar: Operating pressure (bar)

    Returns:
        q: Feed quality factor
    """
    # Bubble point of feed
    T_bp = bubble_point_temperature(x_F, P_bar)

    # Latent heat at feed composition
    lambda_F = mixture_latent_heat(x_F)  # J/mol

    # Heat capacity (need to convert to J/(mol·K))
    Cp_mass = liquid_heat_capacity(x_F)  # J/(kg·K)
    MW = mixture_molecular_weight(x_F)  # g/mol
    Cp_mol = Cp_mass * MW / 1000  # J/(mol·K)

    # Heat needed to bring feed to bubble point and vaporize
    # q = 1 + Cp * (T_bp - T_feed) / lambda
    q = 1 + Cp_mol * (T_bp - T_feed) / lambda_F

    return q


def q_line(x, x_F, q):
    """
    Calculate y on the q-line for given x.

    q-line equation: y = q/(q-1) * x - x_F/(q-1)

    Parameters:
        x: Liquid mole fraction
        x_F: Feed mole fraction
        q: Feed quality factor

    Returns:
        y: Vapor mole fraction on q-line
    """
    if abs(q - 1) < 1e-10:
        # Saturated liquid feed - vertical line at x = x_F
        return None
    return q / (q - 1) * x - x_F / (q - 1)


def rectifying_line(x, x_D, R):
    """
    Rectifying operating line.

    y = R/(R+1) * x + x_D/(R+1)

    Parameters:
        x: Liquid mole fraction
        x_D: Distillate mole fraction
        R: Reflux ratio (L/D)

    Returns:
        y: Vapor mole fraction on rectifying line
    """
    return R / (R + 1) * x + x_D / (R + 1)


def stripping_line(x, x_W, x_D, R, q, x_F):
    """
    Stripping operating line.

    Passes through (x_W, x_W) and intersection of rectifying line with q-line.

    Parameters:
        x: Liquid mole fraction
        x_W: Bottom product mole fraction
        x_D: Distillate mole fraction
        R: Reflux ratio
        q: Feed quality factor
        x_F: Feed mole fraction

    Returns:
        y: Vapor mole fraction on stripping line
    """
    # Find intersection of rectifying line and q-line
    # R/(R+1) * x + x_D/(R+1) = q/(q-1) * x - x_F/(q-1)
    # Solve for x_intersect
    if abs(q - 1) < 1e-10:
        x_int = x_F
        y_int = rectifying_line(x_int, x_D, R)
    else:
        # (q/(q-1) - R/(R+1)) * x = x_D/(R+1) + x_F/(q-1)
        slope_q = q / (q - 1)
        slope_r = R / (R + 1)
        intercept_q = -x_F / (q - 1)
        intercept_r = x_D / (R + 1)

        x_int = (intercept_r - intercept_q) / (slope_q - slope_r)
        y_int = rectifying_line(x_int, x_D, R)

    # Stripping line passes through (x_W, x_W) and (x_int, y_int)
    slope = (y_int - x_W) / (x_int - x_W) if abs(x_int - x_W) > 1e-10 else 0

    return x_W + slope * (x - x_W)


def find_equilibrium_y(x, x_eq, y_eq):
    """
    Interpolate equilibrium y value for given x.

    Parameters:
        x: Liquid mole fraction
        x_eq: Array of equilibrium x values
        y_eq: Array of equilibrium y values

    Returns:
        y: Equilibrium vapor mole fraction
    """
    return np.interp(x, x_eq, y_eq)


def find_minimum_reflux(x_F, x_D, x_W, q, x_eq, y_eq):
    """
    Find minimum reflux ratio using the Underwood method (graphical approach).

    At minimum reflux, the operating line passes through the pinch point
    (intersection of q-line and equilibrium curve).

    Parameters:
        x_F: Feed mole fraction
        x_D: Distillate mole fraction
        x_W: Bottom product mole fraction
        q: Feed quality factor
        x_eq, y_eq: Equilibrium curve data

    Returns:
        R_min: Minimum reflux ratio
    """
    # Find pinch point: intersection of q-line with equilibrium curve
    # y_eq(x) = q/(q-1) * x - x_F/(q-1)

    def objective(x):
        y_eq_val = find_equilibrium_y(x, x_eq, y_eq)
        y_q = q_line(x, x_F, q)
        if y_q is None:
            return y_eq_val - x_F  # For q=1, intersection at x=x_F
        return y_eq_val - y_q

    # Search for intersection in the feed region
    try:
        x_pinch = brentq(objective, x_W + 0.01, x_D - 0.01)
    except ValueError:
        # If no intersection found, use feed composition
        x_pinch = x_F

    y_pinch = find_equilibrium_y(x_pinch, x_eq, y_eq)

    # At minimum reflux, rectifying line passes through (x_D, x_D) and (x_pinch, y_pinch)
    # y = R/(R+1) * x + x_D/(R+1)
    # At (x_pinch, y_pinch): y_pinch = R_min/(R_min+1) * x_pinch + x_D/(R_min+1)
    # Solve for R_min:
    # y_pinch * (R_min + 1) = R_min * x_pinch + x_D
    # y_pinch * R_min + y_pinch = R_min * x_pinch + x_D
    # R_min * (y_pinch - x_pinch) = x_D - y_pinch
    # R_min = (x_D - y_pinch) / (y_pinch - x_pinch)

    if abs(y_pinch - x_pinch) < 1e-10:
        R_min = float("inf")
    else:
        R_min = (x_D - y_pinch) / (y_pinch - x_pinch)

    return R_min, x_pinch, y_pinch


def count_theoretical_stages(x_D, x_W, R, q, x_F, x_eq, y_eq, max_stages=50):
    """
    Count theoretical stages using McCabe-Thiele stepping.

    Parameters:
        x_D: Distillate mole fraction
        x_W: Bottom product mole fraction
        R: Reflux ratio
        q: Feed quality factor
        x_F: Feed mole fraction
        x_eq, y_eq: Equilibrium curve data
        max_stages: Maximum stages to count

    Returns:
        N: Number of theoretical stages
        stages_x: x values at each stage
        stages_y: y values at each stage
    """
    stages_x = [x_D]
    stages_y = [x_D]

    x = x_D
    y = x_D

    # Determine feed stage location
    # Find intersection of q-line and rectifying line
    if abs(q - 1) < 1e-10:
        x_feed_intersect = x_F
    else:
        slope_q = q / (q - 1)
        slope_r = R / (R + 1)
        intercept_q = -x_F / (q - 1)
        intercept_r = x_D / (R + 1)
        x_feed_intersect = (intercept_r - intercept_q) / (slope_q - slope_r)

    N = 0
    in_stripping = False

    while x > x_W and N < max_stages:
        # Step horizontally to equilibrium curve
        # Find x such that y_eq(x) = y
        # Binary search / interpolation
        def eq_objective(x_try):
            return find_equilibrium_y(x_try, x_eq, y_eq) - y

        try:
            x_new = brentq(eq_objective, x_W, x)
        except ValueError:
            break

        stages_x.append(x_new)
        stages_y.append(y)

        # Check if we've crossed into stripping section
        if not in_stripping and x_new < x_feed_intersect:
            in_stripping = True

        # Step vertically to operating line
        if in_stripping:
            y_new = stripping_line(x_new, x_W, x_D, R, q, x_F)
        else:
            y_new = rectifying_line(x_new, x_D, R)

        stages_x.append(x_new)
        stages_y.append(y_new)

        x = x_new
        y = y_new
        N += 1

        if x <= x_W:
            break

    return N, stages_x, stages_y


def calculate_flow_rates(F, x_F, x_D, x_W, R):
    """
    Calculate flow rates from material balance.

    Parameters:
        F: Feed flow rate (mol/h or L/h)
        x_F: Feed mole fraction
        x_D: Distillate mole fraction
        x_W: Bottom mole fraction
        R: Reflux ratio

    Returns:
        D: Distillate flow rate
        W: Bottom flow rate
        L: Reflux flow rate (in rectifying section)
        V: Vapor flow rate (in rectifying section)
    """
    # Material balance: F = D + W
    # Component balance: F * x_F = D * x_D + W * x_W
    # Solve: D = F * (x_F - x_W) / (x_D - x_W)

    D = F * (x_F - x_W) / (x_D - x_W)
    W = F - D

    # Reflux: L = R * D
    L = R * D

    # Vapor: V = L + D = (R + 1) * D
    V = (R + 1) * D

    return D, W, L, V


def run_mccabe_thiele_analysis():
    """
    Run complete McCabe-Thiele analysis and print results.
    """
    print("=" * 60)
    print("McCabe-Thiele Analysis for Ethanol-Water Distillation")
    print("=" * 60)

    # Convert volume fractions to mole fractions
    x_F = vol_to_mol_fraction(FEED_CONC_VOL)
    x_D = vol_to_mol_fraction(PRODUCT_CONC_VOL)
    x_W = vol_to_mol_fraction(BOTTOM_CONC_VOL)

    print(f"\nComposition Specifications:")
    print(f"  Feed:     {FEED_CONC_VOL * 100:.0f}% vol = {x_F * 100:.2f}% mol")
    print(f"  Distillate: {PRODUCT_CONC_VOL * 100:.0f}% vol = {x_D * 100:.2f}% mol")
    print(f"  Bottom:   {BOTTOM_CONC_VOL * 100:.0f}% vol = {x_W * 100:.2f}% mol")

    # Calculate q-factor for subcooled feed
    q = calculate_q_factor(x_F, FEED_TEMP, OPERATING_PRESSURE)
    T_bp = bubble_point_temperature(x_F, OPERATING_PRESSURE)

    print(f"\nFeed Condition:")
    print(f"  Feed temperature: {FEED_TEMP:.0f}°C")
    print(f"  Bubble point: {T_bp:.1f}°C")
    print(
        f"  Feed is {'subcooled' if FEED_TEMP < T_bp else 'at or above bubble point'}"
    )
    print(f"  q-factor: {q:.3f}")

    # Generate equilibrium curve
    x_eq, y_eq, T_eq = equilibrium_curve(OPERATING_PRESSURE, n_points=200)

    # Find minimum reflux ratio
    R_min, x_pinch, y_pinch = find_minimum_reflux(x_F, x_D, x_W, q, x_eq, y_eq)

    print(f"\nMinimum Reflux:")
    print(f"  Pinch point: x = {x_pinch * 100:.2f}%, y = {y_pinch * 100:.2f}%")
    print(f"  R_min = {R_min:.3f}")

    # Find reflux ratio that achieves separation with available trays
    # Use binary search to find R that gives N_theoretical = NUM_TRAYS * TRAY_EFFICIENCY
    target_stages = NUM_TRAYS * TRAY_EFFICIENCY

    print(
        f"\nFinding reflux ratio for {NUM_TRAYS} trays ({target_stages:.1f} theoretical stages)..."
    )

    R_low = R_min * 1.01  # Just above minimum
    R_high = R_min * 10  # High reflux

    for _ in range(50):  # Binary search
        R_test = (R_low + R_high) / 2
        N_test, _, _ = count_theoretical_stages(x_D, x_W, R_test, q, x_F, x_eq, y_eq)

        if abs(N_test - target_stages) < 0.1:
            break
        elif N_test > target_stages:
            R_low = R_test  # Need more reflux (fewer stages)
        else:
            R_high = R_test  # Can use less reflux

    R = R_test

    print(f"  Required R = {R:.3f} ({R / R_min:.2f} × R_min)")

    # Count theoretical stages with operating R
    N_theor, stages_x, stages_y = count_theoretical_stages(
        x_D, x_W, R, q, x_F, x_eq, y_eq
    )

    print(f"\nOperating Conditions:")
    print(f"  Reflux ratio R = {R:.3f}")
    print(f"  Theoretical stages: {N_theor}")
    print(f"  Actual stages needed: {N_theor / TRAY_EFFICIENCY:.1f}")
    print(f"  Available trays: {NUM_TRAYS}")

    N_actual = N_theor / TRAY_EFFICIENCY

    if N_actual <= NUM_TRAYS + 0.5:  # Allow small margin
        print(f"  ✓ Separation is FEASIBLE")
    else:
        print(f"  ✗ Separation requires more trays!")

    # Calculate flow rates
    D, W, L, V = calculate_flow_rates(FEED_FLOW_RATE, x_F, x_D, x_W, R)

    print(f"\nFlow Rates (based on {FEED_FLOW_RATE:.1f} L/h feed):")
    print(f"  Distillate D = {D:.3f} L/h")
    print(f"  Bottom W = {W:.3f} L/h")
    print(f"  Reflux L = {L:.3f} L/h")
    print(f"  Vapor V = {V:.3f} L/h")
    print(f"  Reflux ratio check: L/D = {L / D:.3f}")

    # Return key results
    return {
        "x_F": x_F,
        "x_D": x_D,
        "x_W": x_W,
        "q": q,
        "R_min": R_min,
        "R": R,
        "N_theoretical": N_theor,
        "N_actual": N_actual,
        "D": D,
        "W": W,
        "L": L,
        "V": V,
        "stages_x": stages_x,
        "stages_y": stages_y,
        "x_eq": x_eq,
        "y_eq": y_eq,
    }


if __name__ == "__main__":
    results = run_mccabe_thiele_analysis()
