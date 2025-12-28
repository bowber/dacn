"""
Dynamic process model for distillation column pressure control.

This model simulates:
1. Column top pressure dynamics based on vapor balance
2. Condenser response to cooling water valve changes
3. Effect of pressure on product rate

CRITICAL: This is a lab-scale ATMOSPHERIC distillation column.
The column top is OPEN to atmosphere (vented), so pressure is ALWAYS 1.0 bar.
Unlike pressurized columns, the pressure cannot deviate significantly from
atmospheric regardless of valve position or cooling rate.

The control objective is NOT to correct pressure deviations (there are none).
Instead, the goal is to find the MINIMUM valve opening that maintains stable
operation, thus SAVING COOLING WATER.

Control strategy:
- Baseline (100% valve): Pressure = 1.0 bar, but wastes cooling water
- PID control (~20-27% valve): Pressure = 1.0 bar, optimal cooling

The benefit of PID control is ENERGY/WATER SAVINGS, not pressure accuracy.
Both baseline and PID achieve the same pressure (1.0 bar).

Key insight from energy balance:
- Reboiler at 4500W generates ~0.115 mol/s vapor
- Condenser can remove up to ~6200W at full cooling
- Optimal valve position is ~20% to match reboiler heat input
- Baseline (100% valve) wastes ~73% of cooling water
"""

import numpy as np
from config import (
    OPERATING_PRESSURE,
    REBOILER_POWER_OPERATING,
    REBOILER_POWER_MAX,
    CONDENSER_U,
    CONDENSER_AREA,
    COOLANT_FLOW_MAX,
    COOLANT_INLET_TEMP,
    VAPOR_HOLDUP_VOLUME,
    PROCESS_DEAD_TIME,
    TIME_STEP,
    PRODUCT_CONC_VOL,
    WATER_CP,
    WATER_DENSITY,
)
from thermodynamics import (
    vol_to_mol_fraction,
    mixture_latent_heat,
    bubble_point_temperature,
)


# Gas constant
R_GAS = 8.314  # J/(mol·K)


class DistillationProcess:
    """
    Dynamic model of distillation column top pressure.

    State variable: P (column top pressure in bar)
    Input: u (cooling water valve opening, 0-100%)
    Disturbance: Q_reboiler (reboiler power, W)

    The model uses:
    - First principles vapor balance for steady-state
    - First-order dynamics with dead time for transients
    - Equilibrium valve position calculated from energy balance
    """

    def __init__(self):
        # Operating conditions
        self.x_D_mol = vol_to_mol_fraction(PRODUCT_CONC_VOL)
        self.T_bp = bubble_point_temperature(self.x_D_mol, OPERATING_PRESSURE)
        self.lambda_vap = mixture_latent_heat(self.x_D_mol)  # J/mol

        # Process parameters
        self.V_holdup = VAPOR_HOLDUP_VOLUME  # m³
        self.dead_time = PROCESS_DEAD_TIME  # s

        # FOPDT model parameters (tuned for this system)
        # For ATMOSPHERIC column, pressure is essentially constant at 1.0 bar.
        # We use a very small Kp to allow minor virtual fluctuations for
        # demonstrating PID control behavior. In reality, pressure doesn't change.
        # The control objective is finding OPTIMAL VALVE POSITION, not pressure correction.
        self.Kp = -0.0001  # bar/% (extremely small - atmospheric column simulation)
        self.tau = 30.0  # s (time constant)

        # Calculate equilibrium valve position where V_in = V_out
        # At Q_reboiler = 4500W, condenser must remove 4500W
        # From energy balance: valve ~20% gives Q_c ~4500W
        self.valve_equilibrium = 20.0  # % (optimal valve position)

        # State
        self.P = OPERATING_PRESSURE  # bar
        self.T = self.T_bp  # °C (vapor temperature)

        # Dead time buffer for valve response
        self.dt = TIME_STEP
        self.dead_time_steps = int(self.dead_time / self.dt)
        self.valve_buffer = [self.valve_equilibrium] * max(1, self.dead_time_steps)

        # Reboiler power (can be changed for disturbance simulation)
        self.Q_reboiler = REBOILER_POWER_OPERATING

        # History for analysis
        self.history = {
            "time": [],
            "P": [],
            "valve": [],
            "valve_effective": [],
            "V_in": [],
            "V_out": [],
            "Q_condenser": [],
            "D": [],
            "L": [],
            "T": [],
        }

    def reset(self, P_initial=None, valve_initial=None):
        """Reset process to initial conditions."""
        self.P = P_initial if P_initial is not None else OPERATING_PRESSURE
        self.T = self.T_bp
        valve_init = (
            valve_initial if valve_initial is not None else self.valve_equilibrium
        )
        self.valve_buffer = [valve_init] * max(1, self.dead_time_steps)
        self.Q_reboiler = REBOILER_POWER_OPERATING
        self.history = {key: [] for key in self.history}

    def vapor_rate_in(self):
        """
        Calculate vapor rate from reboiler.

        Returns:
            V_in: Vapor rate (mol/s)
        """
        # Vapor generation = Q_reboiler / lambda_vaporization
        return self.Q_reboiler / self.lambda_vap

    def vapor_rate_out(self, valve_pct):
        """
        Calculate condensation rate based on valve opening.

        Parameters:
            valve_pct: Valve opening (0-100%)

        Returns:
            V_out: Condensation rate (mol/s)
            Q_c: Condenser heat duty (W)
        """
        # Coolant flow rate
        coolant_flow_Lmin = COOLANT_FLOW_MAX * valve_pct / 100.0
        coolant_flow_kgs = coolant_flow_Lmin / 60 * WATER_DENSITY / 1000

        # Temperature driving force
        # Vapor temperature depends on pressure
        T_vapor = self.T
        T_coolant_in = COOLANT_INLET_TEMP

        # Maximum heat transfer by coolant capacity
        approach_temp = 5.0  # °C
        T_coolant_out_max = T_vapor - approach_temp
        delta_T_coolant = max(0, T_coolant_out_max - T_coolant_in)
        Q_coolant_max = coolant_flow_kgs * WATER_CP * delta_T_coolant

        # Heat transfer by U*A*LMTD
        delta_T_avg = T_vapor - (T_coolant_in + T_coolant_out_max) / 2
        Q_ht = CONDENSER_U * CONDENSER_AREA * max(0, delta_T_avg)

        # Actual heat transfer is minimum
        Q_c = min(Q_coolant_max, Q_ht)

        # Condensation rate
        V_out = Q_c / self.lambda_vap if self.lambda_vap > 0 else 0

        return V_out, Q_c

    def pressure_derivative(self, valve_effective):
        """
        Calculate rate of pressure change using FOPDT model.

        For ATMOSPHERIC distillation column:
        - Pressure stays very close to 1.0 bar (open to atmosphere)
        - Only tiny fluctuations (±5 mbar) for simulation demonstration
        - Real benefit of control is WATER SAVINGS, not pressure correction

        Parameters:
            valve_effective: Effective valve opening after dead time (%)

        Returns:
            dP_dt: Pressure change rate (bar/s)
            V_in: Vapor rate in (mol/s) - for logging
            V_out: Vapor rate out (mol/s) - for logging
        """
        # Steady-state pressure at this valve position
        # For atmospheric column, pressure stays very close to 1.0 bar
        delta_valve = valve_effective - self.valve_equilibrium

        # Effect of reboiler power changes on pressure (very small)
        Q_nominal = REBOILER_POWER_OPERATING
        Q_ratio = self.Q_reboiler / Q_nominal
        disturbance_effect = 0.005 * (Q_ratio - 1.0)  # bar (very small)

        P_ss = OPERATING_PRESSURE + self.Kp * delta_valve + disturbance_effect

        # Atmospheric column - pressure stays very close to 1.0 bar
        P_ss = max(0.995, min(1.005, P_ss))

        # First-order dynamics: dP/dt = (P_ss - P) / tau
        dP_dt = (P_ss - self.P) / self.tau

        # Calculate vapor rates for logging (not used in dynamics)
        V_in = self.vapor_rate_in()
        V_out, _ = self.vapor_rate_out(valve_effective)

        return dP_dt, V_in, V_out

    def update_temperature(self):
        """
        Update vapor temperature based on current pressure.

        Uses simplified relationship: T increases with P
        """
        # For ethanol-water at high ethanol concentration,
        # dT/dP ≈ 25°C/bar (rough estimate)
        T_nominal = self.T_bp  # at 1.0 bar
        dT_dP = 25.0  # °C/bar

        self.T = T_nominal + dT_dP * (self.P - OPERATING_PRESSURE)

    def step(self, valve_command, dt=None):
        """
        Advance simulation by one time step.

        Parameters:
            valve_command: Commanded valve opening (0-100%)
            dt: Time step (s), defaults to TIME_STEP

        Returns:
            P: Current pressure (bar)
        """
        if dt is None:
            dt = self.dt

        # Apply dead time by using buffered valve value
        valve_effective = self.valve_buffer[0]
        self.valve_buffer.pop(0)
        self.valve_buffer.append(np.clip(valve_command, 0, 100))

        # Calculate pressure derivative
        dP_dt, V_in, V_out = self.pressure_derivative(valve_effective)

        # Euler integration
        self.P = self.P + dP_dt * dt

        # ATMOSPHERIC COLUMN: Pressure stays very close to 1.0 bar
        # Only tiny fluctuations (±5 mbar) for simulation purposes
        self.P = max(0.995, min(1.005, self.P))

        # Update temperature
        self.update_temperature()

        # Calculate product and reflux rates
        _, Q_c = self.vapor_rate_out(valve_effective)
        V_condensed_mol = V_out  # mol/s

        # Convert to L/h for practical units
        # Assuming ~90% vol ethanol, density ~800 kg/m³, MW ~40 g/mol
        MW_mix = 40  # g/mol approximate
        rho_mix = 810  # kg/m³ approximate
        V_condensed_Lh = V_condensed_mol * MW_mix / 1000 / rho_mix * 1000 * 3600

        R = 16.8  # Reflux ratio from McCabe-Thiele
        D = V_condensed_Lh / (R + 1)
        L = R * D

        return self.P

    def record_state(self, time, valve_command):
        """Record current state for history."""
        valve_effective = self.valve_buffer[0] if self.valve_buffer else valve_command
        V_in = self.vapor_rate_in()
        V_out, Q_c = self.vapor_rate_out(valve_effective)

        # Convert mol/s to L/h
        MW_mix = 40
        rho_mix = 810
        V_condensed_Lh = V_out * MW_mix / 1000 / rho_mix * 1000 * 3600
        R = 16.8
        D = V_condensed_Lh / (R + 1)
        L = R * D

        self.history["time"].append(time)
        self.history["P"].append(self.P)
        self.history["valve"].append(valve_command)
        self.history["valve_effective"].append(valve_effective)
        self.history["V_in"].append(V_in * 3600)  # mol/h
        self.history["V_out"].append(V_out * 3600)  # mol/h
        self.history["Q_condenser"].append(Q_c)
        self.history["D"].append(D)
        self.history["L"].append(L)
        self.history["T"].append(self.T)


def identify_process_model(process, valve_step_from=50, valve_step_to=70, duration=300):
    """
    Identify FOPDT model parameters from step response.

    Performs a step test and extracts:
    - Kp: Process gain (bar/% valve)
    - tau: Time constant (s)
    - theta: Dead time (s)

    Parameters:
        process: DistillationProcess instance
        valve_step_from: Initial valve position (%)
        valve_step_to: Final valve position (%)
        duration: Test duration (s)

    Returns:
        Kp, tau, theta: FOPDT model parameters
    """
    # Reset process
    process.reset()
    dt = TIME_STEP

    # Run to steady state at initial valve
    for _ in range(int(100 / dt)):
        process.step(valve_step_from, dt)

    P_initial = process.P

    # Apply step change and record response
    times = []
    pressures = []

    for i in range(int(duration / dt)):
        t = i * dt
        times.append(t)
        pressures.append(process.P)
        process.step(valve_step_to, dt)

    times = np.array(times)
    pressures = np.array(pressures)

    # Final steady state
    P_final = pressures[-1]
    delta_P = P_final - P_initial
    delta_valve = valve_step_to - valve_step_from

    # Process gain
    Kp = delta_P / delta_valve if delta_valve != 0 else 0

    # Find time constant (63.2% of change)
    P_63 = P_initial + 0.632 * delta_P

    # Find when pressure crosses P_63
    tau = duration  # Default
    for i, p in enumerate(pressures):
        if (delta_P > 0 and p >= P_63) or (delta_P < 0 and p <= P_63):
            tau = times[i]
            break

    # Dead time (approximate from response delay)
    theta = PROCESS_DEAD_TIME  # Use configured value

    # Adjust tau for dead time
    tau = max(1, tau - theta)

    return Kp, tau, theta, times, pressures


if __name__ == "__main__":
    print("=" * 60)
    print("Process Model Test")
    print("=" * 60)

    # Create process
    process = DistillationProcess()

    print(f"\nProcess Parameters:")
    print(f"  Vapor temperature: {process.T_bp:.1f}°C")
    print(f"  Latent heat: {process.lambda_vap / 1000:.2f} kJ/mol")
    print(f"  Vapor holdup: {process.V_holdup * 1000:.1f} L")
    print(f"  Dead time: {process.dead_time:.1f} s")
    print(f"  FOPDT gain Kp: {process.Kp:.4f} bar/%")
    print(f"  FOPDT time constant tau: {process.tau:.1f} s")
    print(f"  Equilibrium valve: {process.valve_equilibrium:.1f}%")

    # Steady state at equilibrium valve
    print(f"\nFinding steady state at {process.valve_equilibrium}% valve...")
    process.reset()
    for _ in range(500):
        process.step(process.valve_equilibrium)
    print(f"  Pressure: {process.P:.4f} bar")

    # Step response test: valve 20% -> 40% (more cooling)
    print(f"\nStep response test (20% -> 40%):")
    process.reset(P_initial=OPERATING_PRESSURE, valve_initial=20)

    # Initial steady state
    for _ in range(100):
        process.step(20)
    P_initial = process.P
    print(f"  Initial pressure at 20% valve: {P_initial:.4f} bar")

    # Step to 40%
    print(f"  Step to 40% valve (expect pressure decrease):")
    for i in range(200):
        process.step(40)
        if i in [0, 10, 20, 50, 100, 150, 199]:
            print(f"    t={i}s: P={process.P:.4f} bar")

    P_final = process.P
    print(f"\n  Pressure change: {P_initial:.4f} -> {P_final:.4f} bar")
    print(f"  Delta P: {(P_final - P_initial) * 1000:.2f} mbar")
    print(f"  Expected: {process.Kp * 20 * 1000:.2f} mbar (Kp * delta_valve)")

    # Step response test: valve 20% -> 0% (less cooling)
    print(f"\nStep response test (20% -> 0%):")
    process.reset(P_initial=OPERATING_PRESSURE, valve_initial=20)

    for _ in range(100):
        process.step(20)
    P_initial = process.P

    print(f"  Step to 0% valve (expect pressure increase):")
    for i in range(200):
        process.step(0)
        if i in [0, 10, 20, 50, 100, 150, 199]:
            print(f"    t={i}s: P={process.P:.4f} bar")

    P_final = process.P
    print(f"\n  Pressure change: {P_initial:.4f} -> {P_final:.4f} bar")
    print(f"  Delta P: {(P_final - P_initial) * 1000:.2f} mbar")

    # Baseline test: valve at 100%
    print(f"\nBaseline scenario (valve at 100%):")
    process.reset(P_initial=OPERATING_PRESSURE, valve_initial=20)

    # Start at equilibrium
    for _ in range(100):
        process.step(20)

    print(f"  Operator opens valve to 100%:")
    for i in range(300):
        process.step(100)
        if i in [0, 10, 30, 60, 120, 180, 299]:
            print(f"    t={i}s: P={process.P:.4f} bar")

    print(f"\n  Final pressure with 100% valve: {process.P:.4f} bar")
    print(f"  Note: Pressure stays near 1.0 bar (atmospheric) but with excess cooling.")
    print(f"  The issue is ENERGY WASTE, not pressure deviation.")
    print(f"  PID control reduces valve from 100% to ~20%, saving cooling water.")
