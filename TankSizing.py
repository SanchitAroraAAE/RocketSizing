import numpy as np
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

# General Cals
# Tank vals
id_tank = 3.75 * in_to_m          # Tank Inner diameter
od_tank = 4 * in_to_m             # Tank Outer diameter
t_tank = (od_tank-id_tank)/2      # Tank Thickness
p_tank = 776 * psi_to_pa          # Target Tank Pressure

# Bolt Vals
num_fastener = 8                  # Number of fasteners on each side of the tank (from HalfCat, would ideally find this val based on fos)
id_fastener  = 0.2614 * in_to_m   # Fastener minor diamter (HalfCat uses 5/16-24 for the tanks)
od_fastener = 0.3125 * in_to_m    # Fastener major diameter / bolt hole dia.
edge_dist = 0.5 * in_to_m         # Distance from center of bolt to edge of casing

burn_time = 4                     # Burn time in seconds

# Material Properties
# 6061 Aluminum
yield_tensile_6061 = 40*1000 * psi_to_pa   # Yield Tensile Strength 40 ksi
ult_tensile_6061 = 45 * 1000 * psi_to_pa   # Ultimate Tensile Strength 45 ksi
yield_bea_6061 = 56 * 1000 * psi_to_pa     # Yield Bearing Strength 56 ksi
ult_bea_6061 = 88 * 1000 * psi_to_pa       # Ultimate Bearing Strength 88 ksi
shear_6061 = 30 * 1000 * psi_to_pa         # Shear Strength 30 ksi

# Alloy Steel (MCM)
yield_tensile_steel = 120 * 1000 * psi_to_pa   # Yield Tensile Strength 120 ksi (from HalfCat. MCM says 170 ksi)
shear_steel = yield_tensile_steel * 0.6        # Shear Strength ~72ksi (from HalfCat)


# Prop Temps
# Mode selection
if mode == "Hotfire": # Hotfire Input Values
    ox_temp = 295            # NOs temp [K], 0deg C
if mode == "Waterflow": # Water Input values
    water_rho = 1000         # Water density

# Densities (kg/m^3)
if mode == "Hotfire":
    ox_rho = CP.PropsSI ("D", "T", ox_temp, "P", p_tank, "NitrousOxide")    # N2O density [kg/m^3]
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
tank_area_m2 = np.pi * (id_tank / 2)**2
len_ox = vol_ox_total / tank_area_m2
len_fuel = vol_fuel_total / tank_area_m2


# COMPONENT SIZING

# Thickness Sizing

# Bolt Sizing
# Force on each bolt
F_bolt = (np.pi/4) * (id_tank**2) * p_tank / num_fastener

# Shear failiure- Occurs when the fasteners holding the closures to the casing break in shear due to force
# applied perpendicularly from pressure
bolt_shear = F_bolt / ( (np.pi/4) * (id_fastener**2))
bolt_fos = shear_steel / bolt_shear

# Tank Sizing
# Hoop Stress- The tensile stress developed in tank wall in the tangential direction as a function of 
# internal pressure pressing outwards (trying to make the tank a larger circle)
hoop_stress = p_tank * od_tank / (2*t_tank)
hoop_fos = yield_tensile_6061 / hoop_stress

# Axial Stress- The tensile stress deceloped in thge tank wall parallel to the axis of the cylinder as 
# a function of internal pressure stretching the tank (acting at the ends, trying to make the tank longer)
axial_stress = p_tank * od_tank / (4*t_tank)
axial_fos = yield_tensile_6061 / axial_stress

# Tear out- Ocurs when the fasteners tear through the end of the casing via shear failiure of the casing material (alum)
min_dist = edge_dist - od_fastener / 2
bolt_tear_out = F_bolt / (min_dist * 2 * t_tank)
tear_out_fos =  shear_6061 / bolt_tear_out

# Casing Tensile Failiure- Occurs when portion of the casing ebtween the fastener holes is stretched beyond breaking
tensile_stress = (np.pi/4) * (id_tank**2) * p_tank / ((np.pi * (od_tank - t_tank) - num_fastener * od_fastener) * t_tank)
tensile_fos = yield_tensile_6061 / tensile_stress

# Bearing Failiure- Occurs when the foce of the fasteners pushing against the edges of their holes causes the casing 
# material to fail in compression (like squeezing a sandwich till the fillings fall out)
bearing_stress = F_bolt / (od_fastener * t_tank)
bearing_fos = yield_bea_6061 / bearing_stress

# Ensure all FOS are atleast 2
fos_req = 2
safe = False
fos_list = [hoop_fos, axial_fos, bolt_fos, tear_out_fos, tensile_fos, bearing_fos]
if all(fos>fos_req for fos in fos_list):
    safe=True
else:
    print(fos_list)

# OUTPUT
print(f"--- Sizing for a Burn Time of {burn_time} s ---")
print(f"Total Propellant Mass: {mass_total_req:.3f} kg")
print(f"Oxidizer ({mass_ox:.2f} kg) Length: {len_ox * 100:.2f} cm ({len_ox * 39.37:.2f} inches)")
print(f"Fuel ({mass_fuel:.2f} kg) Length: {len_fuel * 100:.2f} cm ({len_fuel * 39.37:.2f} inches)")
print(f"Does everything have a FOS of atleast 2: {safe}")
print(f"The total impulse for the engine is: {total_impulse} Ns")