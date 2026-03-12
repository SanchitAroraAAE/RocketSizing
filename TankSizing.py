import numpy as np
import pandas as pd
from BasicSizing import BasicSizing


# INITIALIZATION
mode = "Hotfire"
sizing = BasicSizing(mode)

# Conversion Factors
in_to_m = 0.0254
lbf_to_N = 4.44822162
lbm_to_kg = 0.453592
psi_to_pa = 6894.76


# INPUT PARAMETERS
isp_expected = sizing.isp
OF = sizing.OF
thrust = sizing.thrust        # [N]
m_dot_total = sizing.m_dot_total
burn_time = 4                 # [s]

# Tank Geometry
id_tank = 3.75 * in_to_m
od_tank = 4.0 * in_to_m
t_tank = (od_tank - id_tank) / 2
p_tank = 776 * psi_to_pa      # Target Tank Pressure

# Fastener Geometry (Based on 5/16-24 specs)
num_fastener = 8
id_fastener = 0.2614 * in_to_m # Minor diameter
od_fastener = 0.3125 * in_to_m # Major diameter / hole dia
edge_dist = 0.5 * in_to_m      # Center of bolt to edge of casing


# MATERIAL PROPERTIES
# 6061-T6 Aluminum
yield_tensile_6061 = 40000 * psi_to_pa
ult_tensile_6061 = 45000 * psi_to_pa
yield_bea_6061 = 56000 * psi_to_pa
ult_bea_6061 = 88000 * psi_to_pa
shear_6061 = 30000 * psi_to_pa

# Alloy Steel Fasteners
yield_tensile_steel = 120000 * psi_to_pa
shear_steel = yield_tensile_steel * 0.6


# FLUID PROPERTIES & VOLUMES
# N2O Density Interpolation
n2o_data = pd.read_excel(r"N20 Densities.xlsx")
ox_temp_f = 70 if mode == "Hotfire" else 68 # Default room temp for water

ox_temps = pd.to_numeric(n2o_data.iloc[:, 0], errors="coerce").to_numpy()
ox_pressures = pd.to_numeric(n2o_data.iloc[:, 1], errors="coerce").to_numpy()
ox_rhos = pd.to_numeric(n2o_data.iloc[:, 2], errors="coerce").to_numpy()

if mode == "Hotfire":
    ox_rho = np.interp(ox_temp_f, ox_temps, ox_rhos)
    fuel_rho = 789  # E98
else:
    ox_rho = 1000   # Water
    fuel_rho = 1000

# Mass Calculations
total_impulse = thrust * burn_time
mass_total_req = m_dot_total * burn_time
mass_fuel = mass_total_req / (OF + 1)
mass_ox = mass_total_req - mass_fuel

# Volume & Length Calculations (with Ullage)
ullage_ox, ullage_fuel = 1.15, 1.10
vol_ox_total = (mass_ox / ox_rho) * ullage_ox
vol_fuel_total = (mass_fuel / fuel_rho) * ullage_fuel

tank_area = np.pi * (id_tank / 2)**2
len_ox = vol_ox_total / tank_area
len_fuel = vol_fuel_total / tank_area


# STRUCTURAL ANALYSIS (FOS)

# Load per bolt
F_total_pressure = (np.pi / 4) * (id_tank**2) * p_tank
F_bolt = F_total_pressure / num_fastener

# Bolt Shear
bolt_shear_stress = F_bolt / ((np.pi / 4) * (id_fastener**2))
bolt_fos = shear_steel / bolt_shear_stress

# Tank Wall Stresses (Thin-Walled Theory)
hoop_stress = (p_tank * od_tank) / (2 * t_tank)
hoop_fos = yield_tensile_6061 / hoop_stress

axial_stress = (p_tank * od_tank) / (4 * t_tank)
axial_fos = yield_tensile_6061 / axial_stress

# Casing Failure Modes
# Tear out: Shear failure of casing between bolt and edge
min_dist = edge_dist - (od_fastener / 2)
bolt_tear_out_stress = F_bolt / (min_dist * 2 * t_tank)
tear_out_fos = shear_6061 / bolt_tear_out_stress

# Tensile: Net section failure between bolt holes
net_circumference = (np.pi * (od_tank - t_tank)) - (num_fastener * od_fastener)
tensile_stress_net = F_total_pressure / (net_circumference * t_tank)
tensile_fos = yield_tensile_6061 / tensile_stress_net

# Bearing: Crushing of casing material at bolt contact
bearing_stress = F_bolt / (od_fastener * t_tank)
bearing_fos = yield_bea_6061 / bearing_stress


# VALIDATION & OUTPUT
fos_req = 2.0
fos_dict = {
    "Hoop": hoop_fos, 
    "Axial": axial_fos, 
    "Bolt Shear": bolt_fos, 
    "Tear Out": tear_out_fos, 
    "Net Tensile": tensile_fos, 
    "Bearing": bearing_fos
}

is_safe = all(fos >= fos_req for fos in fos_dict.values())

print(f"--- Sizing for {burn_time}s Burn ---")
print(f"Total Impulse: {total_impulse:.1f} Ns")
print(f"Propellant Mass: {mass_total_req:.3f} kg")
print(f"Oxidizer Tank Length: {len_ox*100:.2f} cm | {len_ox*39.37:.2f} in")
print(f"Fuel Tank Length:     {len_fuel*100:.2f} cm | {len_fuel*39.37:.2f} in")
print("-" * 30)
print(f"All FOS > {fos_req}: {is_safe}")

if not is_safe:
    print("WARNING: Low Factors of Safety detected:")
    for mode, value in fos_dict.items():
        if value < fos_req:
            print(f"  > {mode}: {value:.2f}")