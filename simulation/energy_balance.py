"""
Energy balance calculations for distillation column.
Includes reboiler and condenser heat duty calculations.

SYSTEMATIC APPROACH:
1. From McCabe-Thiele: D, L, V (product, reflux, vapor flow rates)
2. Q_feed = heat to raise feed from T_feed to bubble point
3. Q_vaporize = V × λ_vapor (heat to generate vapor V)
4. Q_condenser = Q_vaporize (optimal - just enough to condense, no subcooling)
5. Calculate optimal coolant flow from Q = U × A × ΔTLMTD

IMPORTANT: Heat transfer coefficient (U) varies with coolant flow rate!
- Q = U × A × ΔTLMTD where U depends on Reynolds number (Re)
- Lower flow rate → lower Re → lower Nusselt number (Nu) → lower U
- This is the key mechanism by which valve position affects heat removal
"""

import numpy as np
from thermodynamics import (
    vol_to_mol_fraction,
    mol_to_vol_fraction,
    mixture_latent_heat,
    mixture_molecular_weight,
    liquid_heat_capacity,
    bubble_point_temperature,
    liquid_density,
)
from config import (
    REBOILER_POWER_MAX,
    REBOILER_POWER_OPERATING,
    REBOILER_POWER_LOSS,
    CONDENSER_AREA,
    CONDENSER_U,
    COOLANT_FLOW_MAX,
    COOLANT_INLET_TEMP,
    OPERATING_PRESSURE,
    PRODUCT_CONC_VOL,
    FEED_CONC_VOL,
    FEED_TEMP,
    FEED_FLOW_RATE,
    BOTTOM_CONC_VOL,
    WATER_CP,
    WATER_DENSITY,
)


def get_effective_reboiler_power(Q_electrical):
    """
    Calculate effective thermal power from electrical input.

    Accounts for losses due to:
    - Voltage drops in wiring
    - Contact resistance
    - Heat loss to environment

    Parameters:
        Q_electrical: Electrical power input (W)

    Returns:
        Q_effective: Effective thermal power to liquid (W)
    """
    return Q_electrical * (1 - REBOILER_POWER_LOSS)


# =============================================================================
# Heat Transfer Coefficient
# =============================================================================


def calculate_heat_transfer_coefficient(coolant_flow_Lmin):
    """
    Calculate overall heat transfer coefficient based on coolant flow rate.

    The heat transfer coefficient U depends on flow regime:
    - Higher flow rate → higher Reynolds number (Re)
    - Higher Re → higher Nusselt number (Nu)
    - Higher heat transfer coefficient

    For turbulent flow in tubes: Nu ∝ Re^0.8
    Since Re ∝ velocity ∝ flow rate: U ∝ (flow rate)^0.8

    Parameters:
        coolant_flow_Lmin: Coolant flow rate (L/min)

    Returns:
        U: Overall heat transfer coefficient (W/(m²·K))
    """
    if coolant_flow_Lmin <= 0:
        return 0.0

    # Reference: U = 800 W/(m²·K) at maximum flow (7.2 L/min)
    U_ref = CONDENSER_U  # 800 W/(m²·K)
    flow_ref = COOLANT_FLOW_MAX  # 7.2 L/min

    # For turbulent flow: Nu ∝ Re^0.8 ∝ (flow)^0.8
    # U scales approximately with the coolant-side heat transfer coefficient
    flow_ratio = coolant_flow_Lmin / flow_ref

    # U ∝ (flow)^0.8 for turbulent regime
    U = U_ref * (flow_ratio**0.8)

    # Minimum U when flow is very low (natural convection limit)
    U_min = 50.0  # W/(m²·K)

    return max(U, U_min)


# =============================================================================
# Feed Heat Calculation (Q_feed)
# =============================================================================


def calculate_feed_heat(F_Lh, x_F_vol, T_feed, P_bar):
    """
    Calculate heat required to raise feed from T_feed to bubble point.

    Q_feed = F × ρ × Cp × (T_bubble - T_feed)

    Parameters:
        F_Lh: Feed flow rate (L/h)
        x_F_vol: Feed concentration (volume fraction ethanol)
        T_feed: Feed temperature (°C)
        P_bar: Operating pressure (bar)

    Returns:
        Q_feed: Heat required (W)
        T_bubble: Bubble point temperature (°C)
        details: Dictionary with calculation details
    """
    # Convert to mole fraction
    x_F_mol = vol_to_mol_fraction(x_F_vol)

    # Bubble point of feed
    T_bubble = bubble_point_temperature(x_F_mol, P_bar)

    # Feed properties
    rho_feed = liquid_density(x_F_vol)  # kg/m³
    Cp_feed = liquid_heat_capacity(x_F_mol)  # J/(kg·K)

    # Convert flow rate: L/h → kg/s
    F_kgs = F_Lh / 3600 * rho_feed / 1000  # L/h → L/s → m³/s → kg/s

    # Heat required
    delta_T = T_bubble - T_feed
    Q_feed = F_kgs * Cp_feed * delta_T  # W

    details = {
        "F_Lh": F_Lh,
        "F_kgs": F_kgs,
        "x_F_vol": x_F_vol,
        "x_F_mol": x_F_mol,
        "T_feed": T_feed,
        "T_bubble": T_bubble,
        "delta_T": delta_T,
        "rho_feed": rho_feed,
        "Cp_feed": Cp_feed,
    }

    return Q_feed, T_bubble, details


# =============================================================================
# Vapor Heat Calculation (Q_vaporize)
# =============================================================================


def calculate_vapor_heat(V_Lh, x_V_vol):
    """
    Calculate heat required to vaporize vapor stream V.

    Q_vaporize = V × ρ × λ / MW

    In rectifying section, vapor composition ≈ distillate composition (x_D).

    Parameters:
        V_Lh: Vapor flow rate (L/h as liquid equivalent)
        x_V_vol: Vapor composition (volume fraction ethanol)

    Returns:
        Q_vaporize: Heat required (W)
        details: Dictionary with calculation details
    """
    # Convert to mole fraction
    x_V_mol = vol_to_mol_fraction(x_V_vol)

    # Properties at vapor composition
    rho_V = liquid_density(x_V_vol)  # kg/m³
    MW_V = mixture_molecular_weight(x_V_mol)  # g/mol
    lambda_V = mixture_latent_heat(x_V_mol)  # J/mol

    # Convert flow rate: L/h → mol/s
    V_m3s = V_Lh / 3600 / 1000  # L/h → m³/s
    V_kgs = V_m3s * rho_V  # kg/s
    V_mols = V_kgs * 1000 / MW_V  # mol/s

    # Heat required
    Q_vaporize = V_mols * lambda_V  # W

    details = {
        "V_Lh": V_Lh,
        "V_mols": V_mols,
        "x_V_vol": x_V_vol,
        "x_V_mol": x_V_mol,
        "rho_V": rho_V,
        "MW_V": MW_V,
        "lambda_V": lambda_V,
    }

    return Q_vaporize, details


# =============================================================================
# Condenser Heat Transfer Calculation
# =============================================================================


def calculate_LMTD(T_vapor, T_coolant_in, T_coolant_out):
    """
    Calculate Log Mean Temperature Difference for condenser.

    For condenser (vapor condensing at constant T_vapor):
    - Hot side: T_vapor (constant, phase change)
    - Cold side: T_coolant_in → T_coolant_out

    LMTD = (ΔT1 - ΔT2) / ln(ΔT1/ΔT2)
    where ΔT1 = T_vapor - T_coolant_in, ΔT2 = T_vapor - T_coolant_out

    Parameters:
        T_vapor: Condensing vapor temperature (°C)
        T_coolant_in: Coolant inlet temperature (°C)
        T_coolant_out: Coolant outlet temperature (°C)

    Returns:
        LMTD: Log mean temperature difference (K or °C)
    """
    delta_T1 = T_vapor - T_coolant_in
    delta_T2 = T_vapor - T_coolant_out

    if delta_T2 <= 0:
        # Coolant outlet exceeds vapor temp - physically impossible
        return 0.0

    if abs(delta_T1 - delta_T2) < 0.01:
        # Nearly equal - use arithmetic mean
        return (delta_T1 + delta_T2) / 2

    LMTD = (delta_T1 - delta_T2) / np.log(delta_T1 / delta_T2)
    return LMTD


def calculate_optimal_coolant_flow(
    Q_required, T_vapor, T_coolant_in, approach_temp=3.0, max_iterations=50
):
    """
    Calculate optimal coolant flow rate for given heat duty.

    Uses iterative solution since U depends on flow rate:
    Q = U(F) × A × LMTD(F)

    Parameters:
        Q_required: Required heat duty (W)
        T_vapor: Condensing vapor temperature (°C)
        T_coolant_in: Coolant inlet temperature (°C)
        approach_temp: Minimum approach temperature (°C)
        max_iterations: Maximum iteration count

    Returns:
        F_optimal: Optimal coolant flow rate (L/min)
        details: Dictionary with calculation details
    """
    if Q_required <= 0:
        return 0.0, {"converged": True, "Q_actual": 0, "iterations": 0}

    # Minimum coolant outlet temp (approach temperature)
    T_coolant_out_max = T_vapor - approach_temp

    # Initial guess: based on coolant heat capacity
    # Q = m_dot × Cp × ΔT_coolant
    delta_T_coolant_guess = (T_coolant_out_max - T_coolant_in) / 2
    F_guess = Q_required / (
        WATER_CP * delta_T_coolant_guess * WATER_DENSITY / 1000 / 60
    )
    F_guess = max(0.5, min(COOLANT_FLOW_MAX, F_guess))

    F = F_guess
    converged = False

    for i in range(max_iterations):
        # Calculate U at current flow
        U = calculate_heat_transfer_coefficient(F)

        # Calculate coolant outlet temp from heat balance
        # Q = m_dot × Cp × (T_out - T_in)
        m_dot = F / 60 * WATER_DENSITY / 1000  # L/min → kg/s
        if m_dot > 0:
            T_coolant_out = T_coolant_in + Q_required / (m_dot * WATER_CP)
        else:
            T_coolant_out = T_coolant_out_max

        # Check approach temperature constraint
        if T_coolant_out > T_coolant_out_max:
            T_coolant_out = T_coolant_out_max

        # Calculate LMTD
        LMTD = calculate_LMTD(T_vapor, T_coolant_in, T_coolant_out)

        # Calculate Q from heat transfer equation
        Q_ht = U * CONDENSER_AREA * LMTD

        # Calculate required flow from heat balance
        delta_T_coolant = T_coolant_out - T_coolant_in
        if delta_T_coolant > 0:
            F_new = Q_required / (
                WATER_CP * delta_T_coolant * WATER_DENSITY / 1000 / 60
            )
        else:
            F_new = COOLANT_FLOW_MAX

        # Check convergence
        if abs(F_new - F) < 0.01 and abs(Q_ht - Q_required) / Q_required < 0.01:
            converged = True
            break

        # Update with relaxation
        F = 0.5 * F + 0.5 * F_new
        F = max(0.1, min(COOLANT_FLOW_MAX * 1.5, F))
    else:
        i = max_iterations - 1

    # Final calculation
    U_final = calculate_heat_transfer_coefficient(F)
    m_dot_final = F / 60 * WATER_DENSITY / 1000
    T_out_final = (
        T_coolant_in + Q_required / (m_dot_final * WATER_CP)
        if m_dot_final > 0
        else T_coolant_in
    )
    LMTD_final = calculate_LMTD(T_vapor, T_coolant_in, T_out_final)
    Q_actual = U_final * CONDENSER_AREA * LMTD_final

    details = {
        "converged": converged,
        "iterations": i + 1,
        "F_Lmin": F,
        "U": U_final,
        "LMTD": LMTD_final,
        "T_coolant_out": T_out_final,
        "Q_required": Q_required,
        "Q_actual": Q_actual,
        "valve_pct": F / COOLANT_FLOW_MAX * 100,
    }

    return F, details


# =============================================================================
# Condenser Capacity (for simulation)
# =============================================================================


def condenser_capacity(coolant_flow_Lmin, T_coolant_in, T_vapor):
    """
    Calculate maximum condenser heat removal capacity.

    Uses both heat transfer (with flow-dependent U) and coolant heat balance.

    Parameters:
        coolant_flow_Lmin: Coolant flow rate (L/min)
        T_coolant_in: Coolant inlet temperature (°C)
        T_vapor: Vapor temperature (°C)

    Returns:
        Q_max: Maximum heat removal capacity (W)
        T_coolant_out: Coolant outlet temperature (°C)
    """
    if coolant_flow_Lmin <= 0:
        return 0.0, T_coolant_in

    # Convert flow rate to kg/s
    coolant_flow_kgs = coolant_flow_Lmin / 60 * WATER_DENSITY / 1000  # L/min -> kg/s

    # Calculate flow-dependent heat transfer coefficient
    U = calculate_heat_transfer_coefficient(coolant_flow_Lmin)

    # Method 1: Coolant heat capacity limit
    # Q = m_dot * Cp * (T_out - T_in)
    # Maximum T_out should be less than T_vapor (approach temperature)
    approach_temp = 5.0  # °C minimum approach
    T_coolant_out_max = T_vapor - approach_temp

    # Maximum Q from coolant capacity
    Q_coolant_max = coolant_flow_kgs * WATER_CP * (T_coolant_out_max - T_coolant_in)

    # Method 2: Heat transfer limit with flow-dependent U
    # Q = U * A * LMTD
    # For condenser with phase change, use simple delta T approximation
    delta_T_avg = T_vapor - (T_coolant_in + T_coolant_out_max) / 2
    Q_ht = U * CONDENSER_AREA * delta_T_avg

    # Actual Q is minimum of both limits
    Q_max = min(max(0, Q_coolant_max), max(0, Q_ht))

    # Calculate actual outlet temperature
    if coolant_flow_kgs > 0 and Q_max > 0:
        T_coolant_out = T_coolant_in + Q_max / (coolant_flow_kgs * WATER_CP)
    else:
        T_coolant_out = T_coolant_in

    return Q_max, T_coolant_out


def condenser_heat_duty(valve_opening_pct, T_vapor):
    """
    Calculate condenser heat duty for given valve opening.

    Parameters:
        valve_opening_pct: Valve opening (0-100%)
        T_vapor: Vapor temperature at condenser inlet (°C)

    Returns:
        Q_c: Condenser heat duty (W)
        coolant_flow: Actual coolant flow rate (L/min)
    """
    # Coolant flow is proportional to valve opening (linear valve)
    coolant_flow = COOLANT_FLOW_MAX * valve_opening_pct / 100.0  # L/min

    # Get condenser capacity at this flow (with flow-dependent U)
    Q_c, _ = condenser_capacity(coolant_flow, COOLANT_INLET_TEMP, T_vapor)

    return Q_c, coolant_flow


# =============================================================================
# Reboiler Calculations
# =============================================================================


def reboiler_vapor_rate(Q_reboiler, x_bottom_mol):
    """
    Calculate vapor generation rate from reboiler heat input.

    Parameters:
        Q_reboiler: Reboiler heat duty (W)
        x_bottom_mol: Bottom product mole fraction ethanol

    Returns:
        V_mol: Vapor rate (mol/s)
        V_Lh: Vapor rate (L/h as liquid equivalent)
    """
    # Latent heat of vaporization at bottom composition
    lambda_mix = mixture_latent_heat(x_bottom_mol)  # J/mol

    # Vapor rate in mol/s
    V_mol = Q_reboiler / lambda_mix  # mol/s

    # Convert to L/h (as liquid equivalent for comparison)
    MW = mixture_molecular_weight(x_bottom_mol)  # g/mol
    x_vol = 0.01  # Approximate bottom as ~1% vol
    rho = liquid_density(x_vol)  # kg/m³

    # mol/s -> g/s -> kg/s -> L/s -> L/h
    mass_rate = V_mol * MW / 1000  # kg/s
    vol_rate = mass_rate / rho * 1000  # L/s
    V_Lh = vol_rate * 3600  # L/h

    return V_mol, V_Lh


def condensation_rate(Q_condenser, x_distillate_mol):
    """
    Calculate condensation rate from condenser heat removal.

    Parameters:
        Q_condenser: Condenser heat duty (W)
        x_distillate_mol: Distillate mole fraction ethanol

    Returns:
        V_condensed_mol: Condensation rate (mol/s)
        V_condensed_Lh: Condensation rate (L/h)
    """
    # Latent heat at distillate composition
    lambda_mix = mixture_latent_heat(x_distillate_mol)  # J/mol

    # Condensation rate
    V_condensed_mol = Q_condenser / lambda_mix  # mol/s

    # Convert to L/h
    MW = mixture_molecular_weight(x_distillate_mol)  # g/mol
    x_vol = 0.90  # Distillate ~90% vol
    rho = liquid_density(x_vol)  # kg/m³

    mass_rate = V_condensed_mol * MW / 1000  # kg/s
    vol_rate = mass_rate / rho * 1000  # L/s
    V_condensed_Lh = vol_rate * 3600  # L/h

    return V_condensed_mol, V_condensed_Lh


def calculate_product_rate(V_condensed_Lh, R):
    """
    Calculate product (distillate) rate from condensation and reflux ratio.

    V_condensed = L + D = R*D + D = (R+1)*D
    D = V_condensed / (R+1)

    Parameters:
        V_condensed_Lh: Total condensation rate (L/h)
        R: Reflux ratio

    Returns:
        D: Product (distillate) rate (L/h)
        L: Reflux rate (L/h)
    """
    D = V_condensed_Lh / (R + 1)
    L = R * D
    return D, L


# =============================================================================
# Complete Energy Balance Analysis
# =============================================================================


def systematic_energy_balance():
    """
    Perform systematic energy balance calculation.

    Steps:
    1. Get McCabe-Thiele results: D, L, V, R
    2. Calculate Q_feed (heat to raise feed to bubble point)
    3. Calculate Q_vaporize from vapor rate V
    4. Set Q_condenser = Q_vaporize (optimal operation)
    5. Calculate optimal coolant flow rate

    Returns:
        results: Dictionary with all calculated values
    """
    print("=" * 70)
    print("SYSTEMATIC ENERGY BALANCE ANALYSIS")
    print("=" * 70)

    # -------------------------------------------------------------------------
    # Step 1: Get McCabe-Thiele results
    # -------------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("Step 1: McCabe-Thiele Material Balance")
    print("-" * 70)

    # Import here to avoid circular dependency
    from mccabe_thiele import run_mccabe_thiele_analysis

    mt_results = run_mccabe_thiele_analysis()

    D_Lh = mt_results["D"]  # L/h
    W_Lh = mt_results["W"]  # L/h
    L_Lh = mt_results["L"]  # L/h
    V_Lh = mt_results["V"]  # L/h
    R = mt_results["R"]

    print(f"\nMcCabe-Thiele Results:")
    print(f"  Feed F = {FEED_FLOW_RATE:.2f} L/h")
    print(f"  Distillate D = {D_Lh:.3f} L/h")
    print(f"  Bottom W = {W_Lh:.3f} L/h")
    print(f"  Reflux L = {L_Lh:.3f} L/h")
    print(f"  Vapor V = {V_Lh:.3f} L/h")
    print(f"  Reflux ratio R = {R:.2f}")

    # -------------------------------------------------------------------------
    # Step 2: Calculate Q_feed
    # -------------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("Step 2: Feed Heating (Q_feed)")
    print("-" * 70)

    Q_feed, T_bubble, feed_details = calculate_feed_heat(
        FEED_FLOW_RATE, FEED_CONC_VOL, FEED_TEMP, OPERATING_PRESSURE
    )

    print(f"\nFeed Conditions:")
    print(
        f"  Flow rate: {FEED_FLOW_RATE:.2f} L/h = {feed_details['F_kgs'] * 1000:.2f} g/s"
    )
    print(
        f"  Concentration: {FEED_CONC_VOL * 100:.0f}% vol = {feed_details['x_F_mol'] * 100:.2f}% mol"
    )
    print(f"  Temperature: {FEED_TEMP:.0f}°C")
    print(f"  Bubble point: {T_bubble:.1f}°C")
    print(f"  ΔT to bubble: {feed_details['delta_T']:.1f}°C")
    print(f"  Heat capacity: {feed_details['Cp_feed']:.0f} J/(kg·K)")
    print(f"\n  Q_feed = {Q_feed:.1f} W")

    # -------------------------------------------------------------------------
    # Step 3: Calculate Q_vaporize
    # -------------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("Step 3: Vaporization Heat (Q_vaporize)")
    print("-" * 70)

    # Vapor composition in rectifying section ≈ distillate composition
    Q_vaporize, vapor_details = calculate_vapor_heat(V_Lh, PRODUCT_CONC_VOL)

    print(f"\nVapor Stream (rectifying section):")
    print(f"  Flow rate: {V_Lh:.3f} L/h = {vapor_details['V_mols'] * 3600:.2f} mol/h")
    print(
        f"  Composition: {PRODUCT_CONC_VOL * 100:.0f}% vol = {vapor_details['x_V_mol'] * 100:.2f}% mol"
    )
    print(f"  Latent heat: {vapor_details['lambda_V'] / 1000:.2f} kJ/mol")
    print(f"\n  Q_vaporize = {Q_vaporize:.1f} W")

    # -------------------------------------------------------------------------
    # Step 4: Total Reboiler Heat
    # -------------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("Step 4: Total Reboiler Heat Duty")
    print("-" * 70)

    Q_reboiler_total = Q_feed + Q_vaporize

    print(f"\nEnergy Balance at Reboiler:")
    print(f"  Q_feed (heating feed to bubble point): {Q_feed:.1f} W")
    print(f"  Q_vaporize (generating vapor V): {Q_vaporize:.1f} W")
    print(f"  -----------------------------------------")
    print(f"  Q_reboiler_total = {Q_reboiler_total:.1f} W")
    print(f"\n  Note: This is the MINIMUM reboiler power needed.")
    print(f"  Current setting: {REBOILER_POWER_OPERATING:.0f} W")

    if Q_reboiler_total > REBOILER_POWER_OPERATING:
        print(f"  WARNING: Required power exceeds operating setting!")
        print(f"  Need to increase reboiler power or reduce feed rate.")

    # -------------------------------------------------------------------------
    # Step 5: Condenser Heat (optimal - no subcooling)
    # -------------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("Step 5: Condenser Heat Duty (Optimal)")
    print("-" * 70)

    # For optimal operation: Q_condenser = Q_vaporize
    # (just enough to condense vapor, no subcooling)
    Q_condenser_optimal = Q_vaporize

    print(f"\nOptimal Condenser Operation:")
    print(f"  Q_condenser = Q_vaporize = {Q_condenser_optimal:.1f} W")
    print(f"  (Just enough to condense vapor, no subcooling)")

    # -------------------------------------------------------------------------
    # Step 6: Calculate Optimal Coolant Flow
    # -------------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("Step 6: Optimal Coolant Flow Rate")
    print("-" * 70)

    # Vapor temperature at condenser (dew point of distillate)
    x_D_mol = vol_to_mol_fraction(PRODUCT_CONC_VOL)
    T_vapor = bubble_point_temperature(x_D_mol, OPERATING_PRESSURE)

    F_optimal, coolant_details = calculate_optimal_coolant_flow(
        Q_condenser_optimal, T_vapor, COOLANT_INLET_TEMP
    )

    print(f"\nCondenser Conditions:")
    print(f"  Vapor temperature: {T_vapor:.1f}°C")
    print(f"  Coolant inlet: {COOLANT_INLET_TEMP:.0f}°C")
    print(f"  Condenser area: {CONDENSER_AREA:.3f} m²")

    print(f"\nHeat Transfer Calculation:")
    print(f"  Q = U × A × LMTD")
    print(
        f"  {Q_condenser_optimal:.1f} = {coolant_details['U']:.0f} × {CONDENSER_AREA:.3f} × {coolant_details['LMTD']:.1f}"
    )

    print(f"\nOptimal Coolant Flow:")
    print(f"  Flow rate: {F_optimal:.2f} L/min")
    print(f"  Valve opening: {coolant_details['valve_pct']:.1f}%")
    print(f"  Coolant outlet temp: {coolant_details['T_coolant_out']:.1f}°C")
    print(f"  Heat transfer coefficient U: {coolant_details['U']:.0f} W/(m²·K)")
    print(f"  LMTD: {coolant_details['LMTD']:.1f}°C")

    # -------------------------------------------------------------------------
    # Step 7: Comparison with Baseline (100% valve)
    # -------------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("Step 7: Comparison with Baseline (100% valve)")
    print("-" * 70)

    Q_baseline, _ = condenser_heat_duty(100.0, T_vapor)
    F_baseline = COOLANT_FLOW_MAX

    print(f"\nBaseline Operation (valve 100%):")
    print(f"  Coolant flow: {F_baseline:.1f} L/min")
    print(f"  Condenser capacity: {Q_baseline:.0f} W")
    print(f"  Excess cooling: {Q_baseline - Q_vaporize:.0f} W (wasted)")

    print(f"\nOptimal Operation (valve {coolant_details['valve_pct']:.0f}%):")
    print(f"  Coolant flow: {F_optimal:.2f} L/min")
    print(f"  Condenser duty: {Q_condenser_optimal:.0f} W")

    coolant_savings = (F_baseline - F_optimal) / F_baseline * 100
    print(f"\nSavings:")
    print(
        f"  Coolant flow reduction: {F_baseline - F_optimal:.2f} L/min ({coolant_savings:.0f}%)"
    )

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    results = {
        # McCabe-Thiele
        "D_Lh": D_Lh,
        "W_Lh": W_Lh,
        "L_Lh": L_Lh,
        "V_Lh": V_Lh,
        "R": R,
        # Feed heating
        "Q_feed": Q_feed,
        "T_feed": FEED_TEMP,
        "T_bubble": T_bubble,
        # Vaporization
        "Q_vaporize": Q_vaporize,
        "T_vapor": T_vapor,
        # Total reboiler
        "Q_reboiler_total": Q_reboiler_total,
        # Optimal condenser
        "Q_condenser_optimal": Q_condenser_optimal,
        "F_optimal": F_optimal,
        "valve_optimal_pct": coolant_details["valve_pct"],
        "U_optimal": coolant_details["U"],
        "LMTD_optimal": coolant_details["LMTD"],
        "T_coolant_out": coolant_details["T_coolant_out"],
        # Baseline
        "Q_baseline": Q_baseline,
        "F_baseline": F_baseline,
        # Savings
        "coolant_savings_pct": coolant_savings,
    }

    print(
        f"\n  Feed: {FEED_FLOW_RATE:.1f} L/h @ {FEED_CONC_VOL * 100:.0f}% vol, {FEED_TEMP}°C"
    )
    print(f"  Product: {D_Lh:.3f} L/h @ {PRODUCT_CONC_VOL * 100:.0f}% vol")
    print(f"  Reflux ratio: {R:.2f}")
    print(f"\n  Q_reboiler = {Q_reboiler_total:.0f} W (Q_feed + Q_vaporize)")
    print(f"  Q_condenser = {Q_condenser_optimal:.0f} W (optimal)")
    print(f"\n  Optimal valve: {coolant_details['valve_pct']:.0f}%")
    print(f"  Coolant savings: {coolant_savings:.0f}%")

    return results


def energy_balance_analysis():
    """
    Perform complete energy balance analysis and print results.
    (Legacy function for compatibility)
    """
    print("=" * 60)
    print("Energy Balance Analysis")
    print("=" * 60)

    # Get compositions
    x_D_mol = vol_to_mol_fraction(PRODUCT_CONC_VOL)
    x_W_mol = vol_to_mol_fraction(0.01)  # 1% vol bottom

    # Vapor temperature (bubble point at distillate composition)
    T_vapor = bubble_point_temperature(x_D_mol, OPERATING_PRESSURE)

    print(f"\nOperating Conditions:")
    print(f"  Pressure: {OPERATING_PRESSURE:.2f} bar")
    print(f"  Vapor temperature: {T_vapor:.1f}°C")
    print(f"  Coolant inlet: {COOLANT_INLET_TEMP:.0f}°C")

    # Reboiler analysis
    print(f"\nReboiler Analysis:")
    print(f"  Maximum power: {REBOILER_POWER_MAX:.0f} W")
    print(f"  Operating power: {REBOILER_POWER_OPERATING:.0f} W")

    V_mol_max, V_Lh_max = reboiler_vapor_rate(REBOILER_POWER_MAX, x_W_mol)
    V_mol_op, V_Lh_op = reboiler_vapor_rate(REBOILER_POWER_OPERATING, x_W_mol)

    print(f"  Max vapor rate: {V_mol_max * 3600:.2f} mol/h ({V_Lh_max:.2f} L/h)")
    print(f"  Operating vapor rate: {V_mol_op * 3600:.2f} mol/h ({V_Lh_op:.2f} L/h)")

    # Condenser analysis at different valve openings
    print(f"\nCondenser Analysis (with flow-dependent U):")
    print(f"  Area: {CONDENSER_AREA:.3f} m²")
    print(f"  U_max (at 100% flow): {CONDENSER_U:.0f} W/(m²·K)")
    print(f"  Max coolant flow: {COOLANT_FLOW_MAX:.1f} L/min")

    print(f"\n  Valve%  | Coolant | U        | Q_cond  | Condensation")
    print(f"          | (L/min) | W/(m²·K) | (W)     | (L/h)")
    print(f"  " + "-" * 55)

    for valve in [20, 40, 60, 80, 100]:
        flow = COOLANT_FLOW_MAX * valve / 100.0
        U = calculate_heat_transfer_coefficient(flow)
        Q_c, _ = condenser_heat_duty(valve, T_vapor)
        _, V_cond_Lh = condensation_rate(Q_c, x_D_mol)
        print(
            f"  {valve:5.0f}   | {flow:5.1f}   | {U:7.0f}  | {Q_c:7.0f} | {V_cond_Lh:6.2f}"
        )

    # Product rate analysis with typical reflux ratio
    print(f"\nProduct Rate vs Valve Opening (R = 16.8):")
    R = 16.8  # From McCabe-Thiele analysis

    print(f"\n  Valve%  | Condensed | Product D | Reflux L")
    print(f"          | (L/h)     | (L/h)     | (L/h)")
    print(f"  " + "-" * 45)

    results = []
    for valve in [20, 40, 60, 80, 100]:
        Q_c, flow = condenser_heat_duty(valve, T_vapor)
        _, V_cond_Lh = condensation_rate(Q_c, x_D_mol)
        D, L = calculate_product_rate(V_cond_Lh, R)
        print(f"  {valve:5.0f}   | {V_cond_Lh:7.2f}   | {D:7.3f}   | {L:6.2f}")
        results.append(
            {"valve": valve, "V_cond": V_cond_Lh, "D": D, "L": L, "Q_c": Q_c}
        )

    # Key insight
    print(f"\n" + "=" * 60)
    print("KEY INSIGHT:")
    print("=" * 60)
    D_100 = results[-1]["D"]
    D_20 = results[0]["D"]
    Q_100 = results[-1]["Q_c"]
    Q_20 = results[0]["Q_c"]

    print(f"  At 100% valve: Q = {Q_100:.0f} W, D = {D_100:.3f} L/h")
    print(f"  At 20% valve:  Q = {Q_20:.0f} W, D = {D_20:.3f} L/h")
    print(f"\n  Heat transfer coefficient U decreases with flow rate!")
    print(f"  Q = U × A × ΔTLMTD, where U ∝ (flow)^0.8")
    print(f"\n  For atmospheric column:")
    print(f"  - Pressure stays at 1.0 bar (open to atmosphere)")
    print(f"  - If Q_condenser < Q_reboiler: excess vapor escapes (product loss!)")
    print(f"  - Optimal: Q_condenser ≈ Q_reboiler = {REBOILER_POWER_OPERATING:.0f} W")

    return results


def actual_operation_analysis():
    """
    Analyze actual operation with configured reboiler power.

    The reboiler power (4500W electrical) minus 10% losses gives effective
    thermal power. This analysis shows optimal valve position for actual
    equipment configuration.
    """
    print("=" * 70)
    print("ACTUAL OPERATION ANALYSIS")
    print("(With configured reboiler power and losses)")
    print("=" * 70)

    # Calculate effective thermal power
    Q_electrical = REBOILER_POWER_OPERATING
    Q_effective = get_effective_reboiler_power(Q_electrical)

    # Get condenser heat duty at different valve positions
    x_D_mol = vol_to_mol_fraction(PRODUCT_CONC_VOL)
    T_vapor = bubble_point_temperature(x_D_mol, OPERATING_PRESSURE)

    print(f"\nReboiler Power:")
    print(f"  Electrical input: {Q_electrical:.0f} W")
    print(f"  Power loss: {REBOILER_POWER_LOSS * 100:.0f}% (voltage drop, heat loss)")
    print(f"  Effective thermal: {Q_effective:.0f} W")

    print(f"\nOperating Conditions:")
    print(f"  Vapor temperature: {T_vapor:.1f}°C")
    print(f"  Coolant inlet: {COOLANT_INLET_TEMP:.0f}°C")
    print(f"  Condenser area: {CONDENSER_AREA:.3f} m²")

    print(f"\nCondenser Performance vs Valve Position:")
    print(f"  Valve | Flow     | U          | Q_cond  | Match")
    print(f"  " + "-" * 55)

    optimal_valve = None
    min_diff = float("inf")

    for valve in range(10, 101, 5):
        Q_c, flow = condenser_heat_duty(valve, T_vapor)
        U = calculate_heat_transfer_coefficient(flow)
        diff = abs(Q_c - Q_effective)

        if diff < min_diff:
            min_diff = diff
            optimal_valve = valve

        marker = ""
        if diff < 200:
            marker = " <-- optimal"

        print(
            f"  {valve:4}% | {flow:.2f} L/min | {U:4.0f} W/(m²K) | {Q_c:5.0f} W{marker}"
        )

    # Calculate results at optimal valve
    Q_optimal, F_optimal = condenser_heat_duty(optimal_valve, T_vapor)
    Q_baseline, F_baseline = condenser_heat_duty(100, T_vapor)

    print(f"\n" + "-" * 70)
    print("COMPARISON")
    print("-" * 70)

    print(f"\nBaseline (100% valve):")
    print(f"  Coolant flow: {F_baseline:.1f} L/min")
    print(f"  Q_condenser: {Q_baseline:.0f} W")
    print(f"  Excess: {Q_baseline - Q_effective:.0f} W (wasted)")

    print(f"\nOptimal ({optimal_valve}% valve):")
    print(f"  Coolant flow: {F_optimal:.1f} L/min")
    print(f"  Q_condenser: {Q_optimal:.0f} W")
    print(f"  Match: {Q_optimal - Q_effective:+.0f} W")

    savings_pct = (F_baseline - F_optimal) / F_baseline * 100
    print(f"\nSavings:")
    print(f"  Coolant reduction: {F_baseline - F_optimal:.1f} L/min")
    print(f"  Percentage: {savings_pct:.0f}%")

    return {
        "Q_electrical": Q_electrical,
        "Q_effective": Q_effective,
        "optimal_valve_pct": optimal_valve,
        "Q_optimal": Q_optimal,
        "F_optimal": F_optimal,
        "Q_baseline": Q_baseline,
        "F_baseline": F_baseline,
        "coolant_savings_pct": savings_pct,
    }


if __name__ == "__main__":
    # Run the systematic analysis (theoretical)
    print("\n" + "=" * 70)
    print("PART 1: THEORETICAL DESIGN (McCabe-Thiele)")
    print("=" * 70)
    results_theory = systematic_energy_balance()

    # Run the actual operation analysis
    print("\n\n")
    print("=" * 70)
    print("PART 2: ACTUAL OPERATION (With 4500W Reboiler)")
    print("=" * 70)
    results_actual = actual_operation_analysis()
