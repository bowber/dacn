"""
Energy balance calculations for distillation column.
Includes reboiler and condenser heat duty calculations.
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

    Uses both heat transfer and coolant heat balance approaches.

    Parameters:
        coolant_flow_Lmin: Coolant flow rate (L/min)
        T_coolant_in: Coolant inlet temperature (°C)
        T_vapor: Vapor temperature (°C)

    Returns:
        Q_max: Maximum heat removal capacity (W)
        T_coolant_out: Coolant outlet temperature (°C)
    """
    # Convert flow rate to kg/s
    coolant_flow_kgs = coolant_flow_Lmin / 60 * WATER_DENSITY / 1000  # L/min -> kg/s

    # Method 1: Heat transfer limit
    # Q = U * A * LMTD
    # For condenser with phase change, use simple delta T approximation
    # LMTD ≈ T_vapor - (T_in + T_out)/2
    # Assume T_out limited to avoid too close approach

    # Method 2: Coolant heat capacity limit
    # Q = m_dot * Cp * (T_out - T_in)
    # Maximum T_out should be less than T_vapor (approach temperature)
    approach_temp = 5.0  # °C minimum approach
    T_coolant_out_max = T_vapor - approach_temp

    # Maximum Q from coolant capacity
    Q_coolant_max = coolant_flow_kgs * WATER_CP * (T_coolant_out_max - T_coolant_in)

    # Heat transfer Q with average delta T
    delta_T_avg = T_vapor - (T_coolant_in + T_coolant_out_max) / 2
    Q_ht = CONDENSER_U * CONDENSER_AREA * delta_T_avg

    # Actual Q is minimum of both limits
    Q_max = min(Q_coolant_max, Q_ht)

    # Calculate actual outlet temperature
    if coolant_flow_kgs > 0:
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

    # Get condenser capacity at this flow
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
    print(f"\nCondenser Analysis:")
    print(f"  Area: {CONDENSER_AREA:.3f} m²")
    print(f"  U: {CONDENSER_U:.0f} W/(m²·K)")
    print(f"  Max coolant flow: {COOLANT_FLOW_MAX:.1f} L/min")

    print(f"\n  Valve%  | Coolant | Q_cond  | Condensation")
    print(f"          | (L/min) | (W)     | (L/h)")
    print(f"  " + "-" * 45)

    for valve in [20, 40, 60, 80, 100]:
        Q_c, flow = condenser_heat_duty(valve, T_vapor)
        _, V_cond_Lh = condensation_rate(Q_c, x_D_mol)
        print(f"  {valve:5.0f}   | {flow:5.1f}   | {Q_c:7.0f} | {V_cond_Lh:6.2f}")

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
    D_60 = results[2]["D"]
    improvement = (D_100 - D_60) / D_60 * 100 if D_60 > 0 else 0

    print(f"  At 100% valve: D = {D_100:.3f} L/h")
    print(f"  At 60% valve:  D = {D_60:.3f} L/h")
    print(
        f"  Product rate is {abs(improvement):.1f}% {'higher' if improvement > 0 else 'lower'} at 100% valve"
    )
    print(f"\n  However, 100% valve causes OVERCOOLING:")
    print(f"  - More condensation than needed for the reflux ratio")
    print(f"  - Pressure drops below setpoint")
    print(f"  - System must increase reboiler power to compensate")
    print(f"\n  With PID control at optimal valve opening:")
    print(f"  - Maintain pressure at setpoint with less cooling")
    print(f"  - Same product rate with less energy waste")

    return results


if __name__ == "__main__":
    results = energy_balance_analysis()
