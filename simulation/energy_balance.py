"""
Energy balance calculations for distillation column.
Includes reboiler and condenser heat duty calculations.

IMPORTANT: Heat transfer coefficient (U) varies with coolant flow rate!
- Q = U × A × ΔTLMTD where U depends on Reynolds number (Re)
- Lower flow rate → lower Re → lower Nusselt number (Nu) → lower U
- This is the key mechanism by which valve position affects heat removal
"""

import numpy as np
from thermodynamics import (
    vol_to_mol_fraction,
    mixture_latent_heat,
    mixture_molecular_weight,
    liquid_heat_capacity,
    bubble_point_temperature,
    liquid_density,
)
from config import (
    REBOILER_POWER_MAX,
    REBOILER_POWER_OPERATING,
    CONDENSER_AREA,
    CONDENSER_U,
    COOLANT_FLOW_MAX,
    COOLANT_INLET_TEMP,
    OPERATING_PRESSURE,
    PRODUCT_CONC_VOL,
    WATER_CP,
    WATER_DENSITY,
)


def calculate_heat_transfer_coefficient(coolant_flow_Lmin):
    """
    Calculate overall heat transfer coefficient based on coolant flow rate.

    The heat transfer coefficient U depends on flow regime:
    - Higher flow rate → higher Reynolds number (Re)
    - Higher Re → higher Nusselt number (Nu)
    - Higher Nu → higher heat transfer coefficient

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


def energy_balance_analysis():
    """
    Perform complete energy balance analysis and print results.
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


if __name__ == "__main__":
    results = energy_balance_analysis()
