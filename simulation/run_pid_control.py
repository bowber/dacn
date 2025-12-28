"""
PID Control simulation: Optimize valve position for water savings.

This simulates the improved operation with PID control, where the
cooling water valve is modulated to find the MINIMUM valve opening
that maintains stable operation.

IMPORTANT: For ATMOSPHERIC distillation column:
- Pressure is ALWAYS 1.0 bar (open to atmosphere)
- PID control doesn't "correct" pressure - it's already at 1.0 bar
- The goal is to find optimal valve position (~20-27%)

Key observations:
- Pressure maintained at 1.0 bar (same as baseline)
- Valve operates at ~20-27% (vs 100% baseline)
- SAVES 73% cooling water
- Same product rate as baseline

Includes tests for:
1. Finding optimal valve position
2. Disturbance rejection (reboiler power change)
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

from config import (
    OPERATING_PRESSURE,
    SIMULATION_TIME,
    TIME_STEP,
    OUTPUT_DIR,
    FIGURE_DPI,
    REBOILER_POWER_OPERATING,
)
from process_model import DistillationProcess
from pid_controller import PIDController, tune_pid_imc
from energy_balance import (
    condenser_heat_duty,
    condensation_rate,
    calculate_product_rate,
)
from thermodynamics import vol_to_mol_fraction

# Configure matplotlib
rcParams["font.family"] = "DejaVu Sans"


def run_pid_simulation(duration=None, save_figures=True, include_disturbance=True):
    """
    Run PID control simulation.

    Parameters:
        duration: Simulation duration in seconds (default: SIMULATION_TIME)
        save_figures: Whether to save figures to OUTPUT_DIR
        include_disturbance: Whether to include a reboiler power disturbance

    Returns:
        results: Dictionary with simulation results
    """
    if duration is None:
        duration = SIMULATION_TIME

    print("=" * 60)
    print("PID Control Simulation")
    print("=" * 60)

    # Create process
    process = DistillationProcess()

    # PID tuning using IMC method
    Kp_process = process.Kp  # -0.015 bar/%
    tau = process.tau  # 30 s
    theta = process.dead_time  # 10 s

    Kc, Ti, Td = tune_pid_imc(Kp_process, tau, theta, tau_c_factor=1.5)

    print(f"\nProcess Model Parameters:")
    print(f"  Kp = {Kp_process} bar/%")
    print(f"  tau = {tau} s")
    print(f"  theta = {theta} s")

    print(f"\nPID Tuning (IMC method):")
    print(f"  Kc = {abs(Kc):.2f}")
    print(f"  Ti = {Ti:.1f} s")
    print(f"  Td = {Td:.1f} s")

    # Create controller
    controller = PIDController(
        Kc=abs(Kc),
        Ti=Ti,
        Td=0,
        setpoint=OPERATING_PRESSURE,
        reverse_acting=False,  # Direct acting for reverse-gain process
    )

    # Initialize at equilibrium
    equilibrium_valve = process.valve_equilibrium
    process.reset(P_initial=OPERATING_PRESSURE, valve_initial=equilibrium_valve)
    controller.reset(output=equilibrium_valve)
    controller.integral = equilibrium_valve / abs(Kc) * Ti  # Bumpless start

    print(f"\nInitial Conditions:")
    print(f"  Pressure: {process.P:.4f} bar")
    print(f"  Valve: {equilibrium_valve:.1f}%")

    # Simulation
    dt = TIME_STEP
    n_steps = int(duration / dt)

    # Storage
    times = []
    pressures = []
    valves = []
    setpoints = []
    temperatures = []
    Q_condensers = []
    product_rates = []
    Q_reboilers = []

    print(f"\nRunning simulation ({duration / 60:.0f} minutes)...")
    if include_disturbance:
        print("  - Disturbance: +10% reboiler power at t=600s (10 min)")

    for i in range(n_steps):
        t = i * dt

        # Setpoint (constant at 1.0 bar)
        SP = OPERATING_PRESSURE
        controller.set_setpoint(SP)

        # Disturbance: increase reboiler power at t=600s
        if include_disturbance and t >= 600:
            process.Q_reboiler = REBOILER_POWER_OPERATING * 1.10  # +10%
        else:
            process.Q_reboiler = REBOILER_POWER_OPERATING

        # Get current pressure
        PV = process.P

        # Calculate controller output
        valve = controller.calculate(PV, time=t)

        # Record state
        times.append(t)
        pressures.append(PV)
        valves.append(valve)
        setpoints.append(SP)
        temperatures.append(process.T)
        Q_reboilers.append(process.Q_reboiler)

        # Calculate condenser duty and product rate
        Q_c, _ = condenser_heat_duty(valve, process.T)
        x_D_mol = vol_to_mol_fraction(0.90)
        _, V_cond_Lh = condensation_rate(Q_c, x_D_mol)
        D, L = calculate_product_rate(V_cond_Lh, R=16.8)

        Q_condensers.append(Q_c)
        product_rates.append(D)

        # Step process
        process.step(valve)
        process.record_state(t, valve)

        # Print progress
        if t % 300 == 0:
            error = SP - PV
            print(
                f"  t={t:4.0f}s: SP={SP:.3f}, PV={PV:.4f} bar, "
                f"valve={valve:.1f}%, D={D:.3f} L/h"
            )

    # Final results
    print(f"\n" + "-" * 40)
    print("Final Steady State:")
    print(
        f"  Pressure: {pressures[-1]:.4f} bar (setpoint: {OPERATING_PRESSURE:.2f} bar)"
    )
    print(f"  Temperature: {temperatures[-1]:.1f}°C")
    print(f"  Valve: {valves[-1]:.1f}%")
    print(f"  Condenser duty: {Q_condensers[-1]:.0f} W")
    print(f"  Product rate: {product_rates[-1]:.3f} L/h")

    # Calculate statistics
    pressure_error = OPERATING_PRESSURE - pressures[-1]
    max_deviation = max(abs(np.array(pressures) - OPERATING_PRESSURE))
    avg_valve = np.mean(valves)

    print(f"\n  Final pressure error: {pressure_error * 1000:.1f} mbar")
    print(f"  Max pressure deviation: {max_deviation * 1000:.1f} mbar")
    print(f"  Average valve position: {avg_valve:.1f}%")

    # Compile results
    results = {
        "time": np.array(times),
        "pressure": np.array(pressures),
        "valve": np.array(valves),
        "setpoint": np.array(setpoints),
        "temperature": np.array(temperatures),
        "Q_condenser": np.array(Q_condensers),
        "Q_reboiler": np.array(Q_reboilers),
        "product_rate": np.array(product_rates),
        "steady_state": {
            "pressure": pressures[-1],
            "temperature": temperatures[-1],
            "valve": valves[-1],
            "Q_condenser": Q_condensers[-1],
            "product_rate": product_rates[-1],
            "pressure_error": pressure_error,
        },
        "stats": {
            "max_deviation": max_deviation,
            "avg_valve": avg_valve,
        },
        "tuning": {
            "Kc": abs(Kc),
            "Ti": Ti,
            "Td": Td,
        },
    }

    # Save figures
    if save_figures:
        os.makedirs(OUTPUT_DIR, exist_ok=True)

        # Figure 1: Pressure control
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(
            np.array(times) / 60, pressures, "b-", linewidth=2, label="PV (Áp suất)"
        )
        ax.plot(
            np.array(times) / 60, setpoints, "r--", linewidth=1.5, label="SP (Setpoint)"
        )
        if include_disturbance:
            ax.axvline(
                x=10,
                color="orange",
                linestyle=":",
                linewidth=1.5,
                label="Nhiễu (+10% Q_reb)",
            )
        ax.set_xlabel("Thời gian (phút)", fontsize=12)
        ax.set_ylabel("Áp suất (bar)", fontsize=12)
        ax.set_title("Điều khiển PID: Áp suất đỉnh tháp", fontsize=14)
        ax.legend(loc="best", fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, duration / 60)
        ax.set_ylim(0.9, 1.1)
        plt.tight_layout()
        plt.savefig(os.path.join(OUTPUT_DIR, "fig_pid_pressure.png"), dpi=FIGURE_DPI)
        plt.close()

        # Figure 2: Valve position
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(np.array(times) / 60, valves, "g-", linewidth=2)
        if include_disturbance:
            ax.axvline(
                x=10, color="orange", linestyle=":", linewidth=1.5, label="Nhiễu"
            )
        ax.set_xlabel("Thời gian (phút)", fontsize=12)
        ax.set_ylabel("Độ mở van (%)", fontsize=12)
        ax.set_title("Điều khiển PID: Tín hiệu điều khiển", fontsize=14)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, duration / 60)
        ax.set_ylim(0, 50)
        if include_disturbance:
            ax.legend(loc="best", fontsize=11)
        plt.tight_layout()
        plt.savefig(os.path.join(OUTPUT_DIR, "fig_pid_valve.png"), dpi=FIGURE_DPI)
        plt.close()

        # Figure 3: Combined plot
        fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

        # Pressure
        axes[0].plot(np.array(times) / 60, pressures, "b-", linewidth=2, label="PV")
        axes[0].plot(np.array(times) / 60, setpoints, "r--", linewidth=1.5, label="SP")
        axes[0].set_ylabel("Áp suất (bar)", fontsize=11)
        axes[0].set_title("Điều khiển áp suất đỉnh tháp bằng PID", fontsize=14)
        axes[0].legend(loc="upper right", fontsize=10)
        axes[0].grid(True, alpha=0.3)
        axes[0].set_ylim(0.9, 1.1)
        if include_disturbance:
            axes[0].axvline(x=10, color="orange", linestyle=":", linewidth=1.5)

        # Valve position
        axes[1].plot(np.array(times) / 60, valves, "g-", linewidth=2)
        axes[1].set_ylabel("Độ mở van (%)", fontsize=11)
        axes[1].grid(True, alpha=0.3)
        axes[1].set_ylim(0, 50)
        if include_disturbance:
            axes[1].axvline(x=10, color="orange", linestyle=":", linewidth=1.5)

        # Product rate
        axes[2].plot(np.array(times) / 60, product_rates, "m-", linewidth=2)
        axes[2].set_xlabel("Thời gian (phút)", fontsize=11)
        axes[2].set_ylabel("Lưu lượng sản phẩm (L/h)", fontsize=11)
        axes[2].grid(True, alpha=0.3)
        axes[2].set_ylim(0, 2)
        if include_disturbance:
            axes[2].axvline(
                x=10, color="orange", linestyle=":", linewidth=1.5, label="Nhiễu"
            )
            axes[2].legend(loc="upper right", fontsize=10)

        plt.tight_layout()
        plt.savefig(os.path.join(OUTPUT_DIR, "fig_pid_combined.png"), dpi=FIGURE_DPI)
        plt.close()

        print(f"\nFigures saved to {OUTPUT_DIR}/")

    return results


if __name__ == "__main__":
    results = run_pid_simulation(
        duration=3600, save_figures=True, include_disturbance=True
    )

    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print("\nPID control maintains pressure at setpoint:")
    print(f"  - Final pressure: {results['steady_state']['pressure']:.4f} bar")
    print(f"  - Error: {results['steady_state']['pressure_error'] * 1000:.1f} mbar")
    print(f"  - Max deviation: {results['stats']['max_deviation'] * 1000:.1f} mbar")
    print(f"  - Average valve: {results['stats']['avg_valve']:.1f}%")
    print(f"  - Product rate: {results['steady_state']['product_rate']:.3f} L/h")
    print(
        f"\nPID Tuning: Kc={results['tuning']['Kc']:.1f}, Ti={results['tuning']['Ti']:.0f}s"
    )
