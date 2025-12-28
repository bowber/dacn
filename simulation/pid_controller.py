"""
PID Controller with anti-windup for distillation column pressure control.

Features:
- Discrete-time implementation suitable for simulation
- Anti-windup using back-calculation method
- Output clamping to valve limits (0-100%)
- Manual/Auto mode switching
- Bumpless transfer

The controller adjusts cooling water valve (MV) to maintain
column top pressure (PV) at setpoint (SP).

Process characteristics (from identification):
- Kp = -0.015 bar/% (negative: more cooling -> lower pressure)
- tau = 30 s
- theta = 10 s

PID tuning using IMC rules for FOPDT:
- Kc = tau / (Kp * (tau_c + theta))
- Ti = tau
- Td = 0 (or small, for noise filtering)

Where tau_c is the desired closed-loop time constant.
"""

import numpy as np
from config import VALVE_MIN, VALVE_MAX, TIME_STEP, PRESSURE_SETPOINT


class PIDController:
    """
    Discrete PID controller with anti-windup.

    Implements the velocity form of PID with back-calculation anti-windup.
    """

    def __init__(
        self,
        Kc=1.0,
        Ti=30.0,
        Td=0.0,
        dt=TIME_STEP,
        output_min=VALVE_MIN,
        output_max=VALVE_MAX,
        setpoint=PRESSURE_SETPOINT,
        reverse_acting=True,
    ):
        """
        Initialize PID controller.

        Parameters:
            Kc: Controller gain
            Ti: Integral time (s), set to inf to disable
            Td: Derivative time (s), set to 0 to disable
            dt: Sample time (s)
            output_min: Minimum output value
            output_max: Maximum output value
            setpoint: Initial setpoint
            reverse_acting: If True, output increases when PV < SP
                           For this system: valve opens more when pressure drops
                           But actually we want: pressure high -> open valve more
                           So reverse_acting=False for direct action
        """
        # Controller parameters
        self.Kc = Kc
        self.Ti = Ti
        self.Td = Td
        self.dt = dt

        # Output limits
        self.output_min = output_min
        self.output_max = output_max

        # Setpoint
        self.SP = setpoint

        # Action direction
        # reverse_acting=True: error = SP - PV (standard)
        # For cooling valve: when P < SP, we want LESS cooling (close valve)
        # So output should DECREASE when error is negative
        # This is direct acting (reverse_acting=False)
        self.reverse_acting = reverse_acting

        # State variables
        self.integral = 0.0
        self.prev_error = 0.0
        self.prev_pv = None
        self.output = (output_min + output_max) / 2  # Start at mid-range

        # Anti-windup tracking
        self.anti_windup_gain = 1.0 / Ti if Ti > 0 else 0.0

        # Mode
        self.auto_mode = True
        self.manual_output = self.output

        # History
        self.history = {
            "time": [],
            "SP": [],
            "PV": [],
            "error": [],
            "P_term": [],
            "I_term": [],
            "D_term": [],
            "output": [],
            "output_clamped": [],
        }

    def reset(self, output=None, integral=None):
        """Reset controller state."""
        if output is not None:
            self.output = output
            self.manual_output = output
        if integral is not None:
            self.integral = integral
        else:
            self.integral = 0.0
        self.prev_error = 0.0
        self.prev_pv = None
        self.history = {key: [] for key in self.history}

    def set_tuning(self, Kc=None, Ti=None, Td=None):
        """Update controller tuning parameters."""
        if Kc is not None:
            self.Kc = Kc
        if Ti is not None:
            self.Ti = Ti
            self.anti_windup_gain = 1.0 / Ti if Ti > 0 else 0.0
        if Td is not None:
            self.Td = Td

    def set_setpoint(self, SP):
        """Change setpoint."""
        self.SP = SP

    def set_mode(self, auto=True, bumpless=True):
        """
        Switch between auto and manual mode.

        Parameters:
            auto: True for automatic, False for manual
            bumpless: If True, perform bumpless transfer
        """
        if auto and not self.auto_mode:
            # Switching to auto - bumpless transfer
            if bumpless:
                # Initialize integral to current output
                self.integral = self.output / self.Kc if self.Kc != 0 else 0
        self.auto_mode = auto

    def calculate(self, PV, time=None):
        """
        Calculate controller output.

        Parameters:
            PV: Process variable (measured pressure)
            time: Current time (for logging)

        Returns:
            output: Controller output (valve opening %)
        """
        if not self.auto_mode:
            return self.manual_output

        # Calculate error
        error = self.SP - PV
        if not self.reverse_acting:
            error = -error

        # Proportional term
        P_term = self.Kc * error

        # Integral term
        if self.Ti > 0:
            self.integral += error * self.dt
            I_term = self.Kc * self.integral / self.Ti
        else:
            I_term = 0.0

        # Derivative term (on PV to avoid derivative kick on SP change)
        if self.Td > 0 and self.prev_pv is not None:
            dPV = PV - self.prev_pv
            D_term = -self.Kc * self.Td * dPV / self.dt
        else:
            D_term = 0.0

        # Calculate unclamped output
        output_unclamped = P_term + I_term + D_term

        # Clamp output
        output_clamped = np.clip(output_unclamped, self.output_min, self.output_max)

        # Anti-windup: back-calculation
        if self.Ti > 0:
            # Reduce integral if output is saturated
            saturation_error = output_clamped - output_unclamped
            self.integral += self.anti_windup_gain * saturation_error * self.dt

        # Store state
        self.prev_error = error
        self.prev_pv = PV
        self.output = output_clamped

        # Record history
        if time is not None:
            self.history["time"].append(time)
            self.history["SP"].append(self.SP)
            self.history["PV"].append(PV)
            self.history["error"].append(error if self.reverse_acting else -error)
            self.history["P_term"].append(P_term)
            self.history["I_term"].append(I_term)
            self.history["D_term"].append(D_term)
            self.history["output"].append(output_unclamped)
            self.history["output_clamped"].append(output_clamped)

        return output_clamped


def tune_pid_imc(Kp_process, tau, theta, tau_c_factor=1.0):
    """
    Tune PID using IMC (Internal Model Control) rules for FOPDT process.

    For FOPDT: G(s) = Kp * exp(-theta*s) / (tau*s + 1)

    IMC-PID tuning:
    - Kc = tau / (Kp * (tau_c + theta))
    - Ti = tau
    - Td = 0.5 * theta (optional, for more aggressive)

    Parameters:
        Kp_process: Process gain (bar/%)
        tau: Process time constant (s)
        theta: Process dead time (s)
        tau_c_factor: Closed-loop time constant factor (default 1.0)
                     Larger = more conservative, Smaller = more aggressive

    Returns:
        Kc, Ti, Td: PID tuning parameters
    """
    # Desired closed-loop time constant
    tau_c = max(0.1 * tau, tau_c_factor * theta)

    # IMC tuning
    Kc = tau / (Kp_process * (tau_c + theta))
    Ti = tau
    Td = 0  # PI control usually sufficient

    return Kc, Ti, Td


def tune_pid_zn(Kp_process, tau, theta):
    """
    Tune PID using Ziegler-Nichols rules for FOPDT process.

    Parameters:
        Kp_process: Process gain (bar/%)
        tau: Process time constant (s)
        theta: Process dead time (s)

    Returns:
        Kc, Ti, Td: PID tuning parameters for PI controller
    """
    # Z-N rules for PI controller
    Kc = 0.9 * tau / (Kp_process * theta)
    Ti = 3.33 * theta
    Td = 0

    return Kc, Ti, Td


if __name__ == "__main__":
    print("=" * 60)
    print("PID Controller Test")
    print("=" * 60)

    # Process parameters (from process_model.py)
    Kp_process = -0.015  # bar/%
    tau = 30.0  # s
    theta = 10.0  # s

    print(f"\nProcess Parameters:")
    print(f"  Kp = {Kp_process} bar/%")
    print(f"  tau = {tau} s")
    print(f"  theta = {theta} s")

    # IMC tuning
    Kc_imc, Ti_imc, Td_imc = tune_pid_imc(Kp_process, tau, theta, tau_c_factor=1.5)
    print(f"\nIMC Tuning (tau_c = 1.5*theta):")
    print(f"  Kc = {Kc_imc:.2f}")
    print(f"  Ti = {Ti_imc:.1f} s")
    print(f"  Td = {Td_imc:.1f} s")

    # Z-N tuning
    Kc_zn, Ti_zn, Td_zn = tune_pid_zn(Kp_process, tau, theta)
    print(f"\nZiegler-Nichols Tuning:")
    print(f"  Kc = {Kc_zn:.2f}")
    print(f"  Ti = {Ti_zn:.1f} s")
    print(f"  Td = {Td_zn:.1f} s")

    # Test controller with simulated process
    print(f"\n" + "=" * 60)
    print("Closed-Loop Simulation Test")
    print("=" * 60)

    from process_model import DistillationProcess

    # Create process and controller
    process = DistillationProcess()

    # Use IMC tuning (more conservative)
    # For this system:
    # - Process: valve increase -> pressure decrease (Kp < 0)
    # - Control: when P < SP (error > 0), we need less cooling (decrease valve)
    # - This requires negative controller gain to make output decrease when error positive
    # OR use direct acting (reverse_acting=False) with positive gain
    #
    # With reverse_acting=False: output = -Kc * error
    # error = SP - PV, if error > 0 (P low), output decreases -> correct!
    controller = PIDController(
        Kc=abs(Kc_imc),  # Use positive gain
        Ti=Ti_imc,
        Td=0,
        setpoint=1.0,
        reverse_acting=False,  # Direct acting for this reverse-gain process
    )

    # Initialize at equilibrium
    process.reset(P_initial=1.0, valve_initial=20)
    controller.reset(output=20)  # Start at equilibrium valve position
    controller.integral = (
        20.0 / abs(Kc_imc) * Ti_imc
    )  # Pre-load integral for bumpless start

    # Simulate with setpoint change
    dt = TIME_STEP
    duration = 600  # 10 minutes

    print(f"\nSimulation: Setpoint step 1.0 -> 0.95 bar at t=60s")
    print(f"Controller: Kc={Kc_imc:.2f}, Ti={Ti_imc:.1f}s")

    times = []
    pressures = []
    valves = []
    setpoints = []

    for i in range(int(duration / dt)):
        t = i * dt

        # Setpoint step at t=60s
        SP = 1.0 if t < 60 else 0.95
        controller.set_setpoint(SP)

        # Get current pressure
        PV = process.P

        # Calculate valve position
        valve = controller.calculate(PV, time=t)

        # Update process
        process.step(valve)
        process.record_state(t, valve)

        # Store for plotting
        times.append(t)
        pressures.append(PV)
        valves.append(valve)
        setpoints.append(SP)

        # Print progress
        if t in [0, 30, 60, 90, 120, 180, 300, 600]:
            print(f"  t={t:3.0f}s: SP={SP:.3f}, PV={PV:.4f} bar, valve={valve:.1f}%")

    # Results
    print(f"\nResults:")
    final_error = abs(setpoints[-1] - pressures[-1])
    print(f"  Final error: {final_error * 1000:.2f} mbar")
    print(f"  Final valve: {valves[-1]:.1f}%")

    # Check if settled
    settled = final_error < 0.005
    print(f"  Settled: {'Yes' if settled else 'No'}")
