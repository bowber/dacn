"""
Configuration parameters for distillation column simulation.
All values based on equipment specifications from Chapter 3.
"""

# =============================================================================
# Operating Conditions
# =============================================================================
OPERATING_PRESSURE = 1.0  # bar (atmospheric)
FEED_TEMP = 70.0  # °C (subcooled feed)
COOLANT_INLET_TEMP = 28.0  # °C
HEAT_LOSS = 0.0  # W (assume insulated)

# =============================================================================
# Feed and Product Specifications
# =============================================================================
FEED_CONC_VOL = 0.10  # 10% vol ethanol
PRODUCT_CONC_VOL = 0.90  # 90% vol ethanol
BOTTOM_CONC_VOL = 0.01  # ~1% vol ethanol (assumed)
FEED_FLOW_RATE = 4.8  # L/h

# =============================================================================
# Column Specifications
# =============================================================================
COLUMN_DIAMETER = 0.150  # m (150 mm)
COLUMN_HEIGHT = 0.800  # m (800 mm)
NUM_TRAYS = 6  # sieve trays
TRAY_EFFICIENCY = 0.80  # 80% Murphree efficiency (good operation)
FEED_TRAY = 4  # from bottom (2 stripping, 4 rectifying)
TRAY_SPACING = 0.140  # m (140 mm)

# =============================================================================
# Reboiler Specifications
# =============================================================================
REBOILER_POWER_MAX = 6000.0  # W (2 x 3000W heaters)
REBOILER_POWER_OPERATING = 6000.0  # W (SSR control - ON/OFF only, no analog)

# Power loss factor (voltage drops, contact resistance, heat loss to environment)
# Effective power = REBOILER_POWER * (1 - REBOILER_POWER_LOSS)
REBOILER_POWER_LOSS = 0.10  # 10% loss
# Effective power = 6000 * 0.9 = 5400W

# =============================================================================
# Condenser Specifications (Coil type)
# =============================================================================
CONDENSER_TUBE_DIAMETER = 0.013  # m (13 mm)
CONDENSER_COIL_TURNS = 22  # turns
CONDENSER_COIL_DIAMETER = 0.100  # m (100 mm)
CONDENSER_AREA = 0.28  # m² (calculated)
COOLANT_FLOW_MAX = 7.2  # L/min

# Estimated overall heat transfer coefficient for coil condenser
# Typical range: 500-1500 W/(m²·K) for steam/water systems
CONDENSER_U = 800.0  # W/(m²·K)

# =============================================================================
# Reflux Drum
# =============================================================================
REFLUX_DRUM_DIAMETER = 0.020  # m (20 mm)
REFLUX_DRUM_HEIGHT = 0.245  # m (245 mm)

# =============================================================================
# Physical Properties - Ethanol
# =============================================================================
ETHANOL_MW = 46.07  # g/mol
ETHANOL_DENSITY = 789.0  # kg/m³ at 20°C
ETHANOL_BP = 78.37  # °C at 1 atm
ETHANOL_LATENT_HEAT = 38.56e3  # J/mol (at boiling point)
ETHANOL_CP = 2.44e3  # J/(kg·K) liquid heat capacity

# Antoine constants for ethanol (log10(P[mmHg]) = A - B/(C + T[°C]))
ETHANOL_ANTOINE_A = 8.20417
ETHANOL_ANTOINE_B = 1642.89
ETHANOL_ANTOINE_C = 230.300

# =============================================================================
# Physical Properties - Water
# =============================================================================
WATER_MW = 18.015  # g/mol
WATER_DENSITY = 998.0  # kg/m³ at 20°C
WATER_BP = 100.0  # °C at 1 atm
WATER_LATENT_HEAT = 40.66e3  # J/mol (at boiling point)
WATER_CP = 4.18e3  # J/(kg·K) liquid heat capacity

# Antoine constants for water (log10(P[mmHg]) = A - B/(C + T[°C]))
WATER_ANTOINE_A = 8.07131
WATER_ANTOINE_B = 1730.63
WATER_ANTOINE_C = 233.426

# =============================================================================
# Azeotrope Data
# =============================================================================
AZEOTROPE_COMP_MOL = 0.8943  # mol fraction ethanol (≈ 95.6 mol%)
AZEOTROPE_COMP_WT = 0.955  # weight fraction (≈ 95.5 wt%)
AZEOTROPE_TEMP = 78.15  # °C at 1 atm

# =============================================================================
# Simulation Parameters
# =============================================================================
SIMULATION_TIME = 3600.0  # s (1 hour)
TIME_STEP = 1.0  # s

# Process dynamics (estimated for lab-scale system)
PROCESS_DEAD_TIME = 10.0  # s (transport delay)
VAPOR_HOLDUP_VOLUME = 0.005  # m³ (estimated column vapor space)

# =============================================================================
# Control Parameters
# =============================================================================
PRESSURE_SETPOINT = 1.0  # bar
VALVE_MIN = 0.0  # % (fully closed)
VALVE_MAX = 100.0  # % (fully open)

# =============================================================================
# Output Settings
# =============================================================================
OUTPUT_DIR = "outputs"
FIGURE_DPI = 150
FIGURE_FORMAT = "png"
