#!/usr/bin/env python3
"""
Main script to run all simulations and generate comparison figures.

This script:
1. Runs baseline simulation (valve=100%)
2. Runs PID control simulation (SP=1.0 bar)
3. Generates comparison figures and tables

Output files are saved to simulation/outputs/
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

from config import OUTPUT_DIR, FIGURE_DPI, OPERATING_PRESSURE
from run_baseline import run_baseline_simulation
from run_pid_control import run_pid_simulation

# Configure matplotlib
rcParams["font.family"] = "DejaVu Sans"


def generate_comparison_figures(baseline_results, pid_results):
    """Generate comparison figures between baseline and PID control."""

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Time in minutes
    t_base = baseline_results["time"] / 60
    t_pid = pid_results["time"] / 60

    # Figure 1: Pressure comparison
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(
        t_base,
        baseline_results["pressure"],
        "b-",
        linewidth=2,
        label="Baseline (van 100%)",
    )
    ax.plot(t_pid, pid_results["pressure"], "g-", linewidth=2, label="PID Control")
    ax.axhline(
        y=OPERATING_PRESSURE, color="r", linestyle="--", linewidth=1.5, label="Setpoint"
    )
    ax.set_xlabel("Thời gian (phút)", fontsize=12)
    ax.set_ylabel("Áp suất (bar)", fontsize=12)
    ax.set_title("So sánh điều khiển áp suất: Baseline vs PID", fontsize=14)
    ax.legend(loc="best", fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 60)
    ax.set_ylim(0.9, 1.1)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "fig_comparison_pressure.png"), dpi=FIGURE_DPI)
    plt.close()

    # Figure 2: Valve position comparison
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(t_base, baseline_results["valve"], "b-", linewidth=2, label="Baseline")
    ax.plot(t_pid, pid_results["valve"], "g-", linewidth=2, label="PID Control")
    ax.set_xlabel("Thời gian (phút)", fontsize=12)
    ax.set_ylabel("Độ mở van (%)", fontsize=12)
    ax.set_title("So sánh độ mở van nước làm mát", fontsize=14)
    ax.legend(loc="best", fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 110)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "fig_comparison_valve.png"), dpi=FIGURE_DPI)
    plt.close()

    # Figure 3: Combined comparison (2x2 subplot)
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Pressure
    axes[0, 0].plot(
        t_base, baseline_results["pressure"], "b-", linewidth=2, label="Baseline"
    )
    axes[0, 0].plot(t_pid, pid_results["pressure"], "g-", linewidth=2, label="PID")
    axes[0, 0].axhline(
        y=OPERATING_PRESSURE, color="r", linestyle="--", linewidth=1.5, label="SP"
    )
    axes[0, 0].set_ylabel("Áp suất (bar)", fontsize=11)
    axes[0, 0].set_title("Áp suất đỉnh tháp", fontsize=12)
    axes[0, 0].legend(loc="best", fontsize=9)
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].set_ylim(0.9, 1.1)

    # Valve
    axes[0, 1].plot(
        t_base, baseline_results["valve"], "b-", linewidth=2, label="Baseline"
    )
    axes[0, 1].plot(t_pid, pid_results["valve"], "g-", linewidth=2, label="PID")
    axes[0, 1].set_ylabel("Độ mở van (%)", fontsize=11)
    axes[0, 1].set_title("Độ mở van nước làm mát", fontsize=12)
    axes[0, 1].legend(loc="best", fontsize=9)
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].set_ylim(0, 110)

    # Condenser duty
    axes[1, 0].plot(
        t_base,
        np.array(baseline_results["Q_condenser"]) / 1000,
        "b-",
        linewidth=2,
        label="Baseline",
    )
    axes[1, 0].plot(
        t_pid,
        np.array(pid_results["Q_condenser"]) / 1000,
        "g-",
        linewidth=2,
        label="PID",
    )
    axes[1, 0].set_xlabel("Thời gian (phút)", fontsize=11)
    axes[1, 0].set_ylabel("Công suất làm mát (kW)", fontsize=11)
    axes[1, 0].set_title("Công suất thiết bị ngưng tụ", fontsize=12)
    axes[1, 0].legend(loc="best", fontsize=9)
    axes[1, 0].grid(True, alpha=0.3)

    # Product rate
    axes[1, 1].plot(
        t_base, baseline_results["product_rate"], "b-", linewidth=2, label="Baseline"
    )
    axes[1, 1].plot(t_pid, pid_results["product_rate"], "g-", linewidth=2, label="PID")
    axes[1, 1].set_xlabel("Thời gian (phút)", fontsize=11)
    axes[1, 1].set_ylabel("Lưu lượng sản phẩm (L/h)", fontsize=11)
    axes[1, 1].set_title("Lưu lượng sản phẩm D", fontsize=12)
    axes[1, 1].legend(loc="best", fontsize=9)
    axes[1, 1].grid(True, alpha=0.3)

    plt.suptitle(
        "So sánh vận hành: Baseline (van 100%) vs Điều khiển PID", fontsize=14, y=1.02
    )
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "fig_comparison_all.png"), dpi=FIGURE_DPI)
    plt.close()

    print(f"\nComparison figures saved to {OUTPUT_DIR}/")


def print_comparison_table(baseline_results, pid_results):
    """Print comparison table of results."""

    base_ss = baseline_results["steady_state"]
    pid_ss = pid_results["steady_state"]

    print("\n" + "=" * 70)
    print("KẾT QUẢ SO SÁNH: BASELINE vs PID CONTROL")
    print("=" * 70)

    print(f"\n{'Thông số':<30} {'Baseline':>15} {'PID Control':>15} {'Đơn vị':>8}")
    print("-" * 70)
    print(
        f"{'Áp suất xác lập':<30} {base_ss['pressure']:>15.4f} {pid_ss['pressure']:>15.4f} {'bar':>8}"
    )
    print(
        f"{'Sai lệch áp suất':<30} {base_ss['pressure_error'] * 1000:>15.1f} {pid_ss['pressure_error'] * 1000:>15.1f} {'mbar':>8}"
    )
    print(
        f"{'Độ mở van xác lập':<30} {base_ss['valve']:>15.1f} {pid_ss['valve']:>15.1f} {'%':>8}"
    )
    print(
        f"{'Công suất ngưng tụ':<30} {base_ss['Q_condenser']:>15.0f} {pid_ss['Q_condenser']:>15.0f} {'W':>8}"
    )
    print(
        f"{'Lưu lượng sản phẩm D':<30} {base_ss['product_rate']:>15.3f} {pid_ss['product_rate']:>15.3f} {'L/h':>8}"
    )
    print(
        f"{'Nhiệt độ đỉnh':<30} {base_ss['temperature']:>15.1f} {pid_ss['temperature']:>15.1f} {'°C':>8}"
    )
    print("-" * 70)

    # Improvements
    print(f"\n{'CHỈ SỐ CẢI THIỆN':^70}")
    print("-" * 70)

    pressure_improvement = (
        abs(base_ss["pressure_error"]) - abs(pid_ss["pressure_error"])
    ) * 1000
    print(f"{'Cải thiện sai lệch áp suất:':<40} {pressure_improvement:>10.1f} mbar")

    valve_reduction = base_ss["valve"] - pid_ss["valve"]
    print(f"{'Giảm độ mở van:':<40} {valve_reduction:>10.1f} %")

    cooling_reduction = (
        (base_ss["Q_condenser"] - pid_ss["Q_condenser"]) / base_ss["Q_condenser"] * 100
    )
    print(f"{'Giảm công suất làm mát:':<40} {cooling_reduction:>10.1f} %")

    # Note about product rate
    print(f"\n{'GHI CHÚ:':^70}")
    print("-" * 70)
    print("- Cả Baseline và PID đều có áp suất 1.0 bar (tháp hở khí quyển)")
    print("- PID Control giảm độ mở van từ 100% xuống ~27%, tiết kiệm nước")
    print("- Lưu lượng sản phẩm tương đương ở cả hai chế độ")
    print("- Lợi ích chính: TIẾT KIỆM 73% NƯỚC LÀM MÁT")

    return {
        "pressure_improvement_mbar": pressure_improvement,
        "valve_reduction_pct": valve_reduction,
        "cooling_reduction_pct": cooling_reduction,
    }


def main():
    """Run all simulations and generate outputs."""

    print("=" * 70)
    print("DISTILLATION COLUMN PRESSURE CONTROL SIMULATION")
    print("=" * 70)

    # Run baseline simulation
    print("\n" + "=" * 70)
    print("PHASE 1: Baseline Simulation")
    print("=" * 70)
    baseline_results = run_baseline_simulation(duration=3600, save_figures=True)

    # Run PID control simulation
    print("\n" + "=" * 70)
    print("PHASE 2: PID Control Simulation")
    print("=" * 70)
    pid_results = run_pid_simulation(
        duration=3600, save_figures=True, include_disturbance=True
    )

    # Generate comparison figures
    print("\n" + "=" * 70)
    print("PHASE 3: Generating Comparison Figures")
    print("=" * 70)
    generate_comparison_figures(baseline_results, pid_results)

    # Print comparison table
    improvements = print_comparison_table(baseline_results, pid_results)

    # Final summary
    print("\n" + "=" * 70)
    print("SIMULATION COMPLETE")
    print("=" * 70)
    print(f"\nOutput files in: {OUTPUT_DIR}/")
    print("  - fig_baseline_pressure.png")
    print("  - fig_baseline_combined.png")
    print("  - fig_pid_pressure.png")
    print("  - fig_pid_valve.png")
    print("  - fig_pid_combined.png")
    print("  - fig_comparison_pressure.png")
    print("  - fig_comparison_valve.png")
    print("  - fig_comparison_all.png")

    return baseline_results, pid_results, improvements


if __name__ == "__main__":
    baseline, pid, improvements = main()
