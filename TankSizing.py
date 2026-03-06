import math
import CoolProp.CoolProp as CP
from BasicSizing import BasicSizing

# RUN BASIC SIZING 
mode = "Hotfire"
sizing = BasicSizing(mode)

# INPUT PARAMETERS
# Conversion Factors
lbf_to_N = 4.44822162
lbm_to_kg = 0.453592
psi_to_pa = 6894.76
in_to_m = 0.0254

# From BasicSizing
isp_expected = sizing.isp           # Predicted Isp in seconds
OF = sizing.OF                      # O/F Ratio
thrust = sizing.thrust   # Thrust values [N], 400 lbf
m_dot_total = sizing.m_dot_total

# General vals
tank_id_in = 3.75            # Tank Inner diameter
burn_time = 4                # Burn time in seconds

# Prop Temps
# Mode selection
if mode == "Hotfire": # Hotfire Input Values
    ox_temp = 298            # NOs temp [K], 0deg C
    P_ox = 850               # Target tank pressure
if mode == "Waterflow": # Water Input values
    water_rho = 1000            # Water density

# Densities (kg/m^3)
if mode == "Hotfire":
    ox_rho = CP.PropsSI ("D", "T", ox_temp, "P", P_ox, "NitrousOxide")    # N2O density [kg/m^3]
    fuel_rho = 789    # E98 density [kg/m^3]
elif mode == "Waterflow":
    ox_rho = 1000     # Water density [kg/m^3]
    fuel_rho = 1000   # Water density [kg/m^3]

# Ullage Factors (1.1 means 10% extra space)
ullage_ox = 1.15
ullage_fuel = 1.10


# CALCULATIONS
# Total Propellant Mass
total_impulse = thrust * burn_time
mass_total_req = m_dot_total * burn_time

# Individual Mass
# mass_ox / mass_fuel = of_ratio -> mass_ox = of_ratio * mass_fuel
# mass_total = mass_fuel * (of_ratio + 1)
mass_fuel = mass_total_req / (OF + 1)
mass_ox = mass_total_req - mass_fuel

# Required Volumes (Cubic Meters)
vol_ox_net = mass_ox / ox_rho
vol_fuel_net = mass_fuel / fuel_rho

# Apply Ullage
vol_ox_total = vol_ox_net * ullage_ox
vol_fuel_total = vol_fuel_net * ullage_fuel

# Tank Lengths
# Volume = Area * Length -> Length = Volume / (pi * r^2)
tank_area_m2 = math.pi * ((tank_id_in * 0.0254) / 2)**2
len_ox_m = vol_ox_total / tank_area_m2
len_fuel_m = vol_fuel_total / tank_area_m2


# OUTPUT
print(f"--- Sizing for a Burn Time of {burn_time} s ---")
print(f"Total Propellant Mass: {mass_total_req:.3f} kg")
print(f"Oxidizer ({mass_ox:.2f} kg) Length: {len_ox_m * 100:.2f} cm ({len_ox_m * 39.37:.2f} inches)")
print(f"Fuel ({mass_fuel:.2f} kg) Length: {len_fuel_m * 100:.2f} cm ({len_fuel_m * 39.37:.2f} inches)")
print(f"The total impulse for the engine is: {total_impulse} Ns")