"""
Baseline simulation: Valve fixed at 100% (maximum cooling).

This simulates the traditional operation where cooling water valve
is fully open at all times, causing EXCESSIVE COOLING and wasting water.

IMPORTANT: For ATMOSPHERIC distillation column:
- Pressure is ALWAYS 1.0 bar (open to atmosphere)
- Baseline and PID both achieve same pressure (1.0 bar)
- The difference is WATER USAGE, not pressure

Key observations:
- Pressure stays at 1.0 bar (atmospheric)
- Valve at 100% wastes cooling water
- Same product rate as PID control

This baseline provides comparison data for the PID control case.
The benefit of PID is WATER SAVINGS (73% reduction), not pressure correction.
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
)
from process_model import DistillationProcess
from energy_balance import (
    condenser_heat_duty,
    condensation_rate,
    calculate_product_rate,
)
from thermodynamics import vol_to_mol_fraction, bubble_point_temperature

# Configure matplotlib for Vietnamese text
rcParams["font.family"] = "DejaVu Sans"


def run_baseline_simulation(duration=None, save_figures=True):
    """
    Run baseline simulation with valve at 100%.

    Parameters:
        duration: Simulation duration in seconds (default: SIMULATION_TIME)
        save_figures: Whether to save figures to OUTPUT_DIR

    Returns:
        results: Dictionary with simulation results
    """
    if duration is None:
        duration = SIMULATION_TIME

    print("=" * 60)
    print("Baseline Simulation: Valve at 100%")
    print("=" * 60)

    # Create process
    process = DistillationProcess()

    # Initialize at equilibrium, then step to 100% valve
    process.reset(P_initial=OPERATING_PRESSURE, valve_initial=20)

    # Run to steady state at equilibrium first (30 seconds)
    print("\nInitializing at equilibrium (valve=20%)...")
    for _ in range(30):
        process.step(20)

    # Now switch to 100% valve (baseline operation)
    print(f"Switching to baseline operation (valve=100%)...")

    dt = TIME_STEP
    n_steps = int(duration / dt)

    # Storage
    times = []
    pressures = []
    valves = []
    temperatures = []
    Q_condensers = []
    product_rates = []

    # Fixed valve at 100%
    valve_baseline = 100.0

    for i in range(n_steps):
        t = i * dt

        # Record state
        process.record_state(t, valve_baseline)

        times.append(t)
        pressures.append(process.P)
        valves.append(valve_baseline)
        temperatures.append(process.T)

        # Calculate condenser duty and product rate
        Q_c, _ = condenser_heat_duty(valve_baseline, process.T)
        x_D_mol = vol_to_mol_fraction(0.90)
        _, V_cond_Lh = condensation_rate(Q_c, x_D_mol)
        D, L = calculate_product_rate(V_cond_Lh, R=16.8)

        Q_condensers.append(Q_c)
        product_rates.append(D)

        # Step process
        process.step(valve_baseline)

        # Print progress
        if t % 300 == 0:
            print(
                f"  t={t:4.0f}s: P={process.P:.4f} bar, T={process.T:.1f}°C, D={D:.3f} L/h"
            )

    # Final results
    print(f"\n" + "-" * 40)
    print("Final Steady State:")
    print(
        f"  Pressure: {pressures[-1]:.4f} bar (setpoint: {OPERATING_PRESSURE:.2f} bar)"
    )
    print(f"  Temperature: {temperatures[-1]:.1f}°C")
    print(f"  Valve: {valves[-1]:.0f}%")
    print(f"  Condenser duty: {Q_condensers[-1]:.0f} W")
    print(f"  Product rate: {product_rates[-1]:.3f} L/h")

    # Calculate deviation from setpoint
    pressure_error = OPERATING_PRESSURE - pressures[-1]
    print(f"\n  Pressure at 1.0 bar (atmospheric column)")
    print(f"  Note: Valve at 100% wastes cooling water")

    # Compile results
    results = {
        "time": np.array(times),
        "pressure": np.array(pressures),
        "valve": np.array(valves),
        "temperature": np.array(temperatures),
        "Q_condenser": np.array(Q_condensers),
        "product_rate": np.array(product_rates),
        "steady_state": {
            "pressure": pressures[-1],
            "temperature": temperatures[-1],
            "valve": valves[-1],
            "Q_condenser": Q_condensers[-1],
            "product_rate": product_rates[-1],
            "pressure_error": pressure_error,
        },
    }

    # Save figures
    if save_figures:
        os.makedirs(OUTPUT_DIR, exist_ok=True)

        # Figure 1: Pressure response
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(np.array(times) / 60, pressures, "b-", linewidth=2, label="Áp suất")
        ax.axhline(
            y=OPERATING_PRESSURE,
            color="r",
            linestyle="--",
            linewidth=1.5,
            label="Setpoint",
        )
        ax.set_xlabel("Thời gian (phút)", fontsize=12)
        ax.set_ylabel("Áp suất (bar)", fontsize=12)
        ax.set_title("Baseline: Áp suất đỉnh tháp (van 100%)", fontsize=14)
        ax.legend(loc="best", fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, duration / 60)
        ax.set_ylim(0.9, 1.1)
        plt.tight_layout()
        plt.savefig(
            os.path.join(OUTPUT_DIR, "fig_baseline_pressure.png"), dpi=FIGURE_DPI
        )
        plt.close()

        # Figure 2: Combined plot
        fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

        # Pressure
        axes[0].plot(np.array(times) / 60, pressures, "b-", linewidth=2)
        axes[0].axhline(
            y=OPERATING_PRESSURE, color="r", linestyle="--", linewidth=1.5, label="SP"
        )
        axes[0].set_ylabel("Áp suất (bar)", fontsize=12)
        axes[0].set_title("Baseline: Van nước làm mát mở 100%", fontsize=14)
        axes[0].legend(loc="upper right", fontsize=10)
        axes[0].grid(True, alpha=0.3)
        axes[0].set_ylim(0.9, 1.1)

        # Product rate
        axes[1].plot(np.array(times) / 60, product_rates, "g-", linewidth=2)
        axes[1].set_xlabel("Thời gian (phút)", fontsize=12)
        axes[1].set_ylabel("Lưu lượng sản phẩm (L/h)", fontsize=12)
        axes[1].grid(True, alpha=0.3)
        axes[1].set_ylim(0, 2)

        plt.tight_layout()
        plt.savefig(
            os.path.join(OUTPUT_DIR, "fig_baseline_combined.png"), dpi=FIGURE_DPI
        )
        plt.close()

        print(f"\nFigures saved to {OUTPUT_DIR}/")

    return results


if __name__ == "__main__":
    results = run_baseline_simulation(duration=3600, save_figures=True)

    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print("\nBaseline operation (100% valve):")
    print(f"  - Pressure: {results['steady_state']['pressure']:.3f} bar (atmospheric)")
    print(f"  - Valve: 100% (wastes cooling water)")
    print(f"  - Product rate: {results['steady_state']['product_rate']:.3f} L/h")
    print("\nWith PID control, valve can be reduced to ~20-27%, saving 73% water.")
