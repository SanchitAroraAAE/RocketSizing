import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from BasicSizing import BasicSizing

# --- RUN BASIC SIZING ---
mode = "Hotfire"
sizing = BasicSizing(mode)

# --- CONSTANTS & CONVERSIONS ---
psi_to_pa = 6894.76
in_to_m = 0.0254

m_dot_fuel = sizing.m_dot_fuel    
m_dot_ox = sizing.m_dot_ox        
OF = sizing.OF
Pc = sizing.Pc
d_c = sizing.d_c                  

# --- USER INPUTS ---
discharge_coef = 0.65     # For radial holes
Cd_annulus = 0.60         # For fuel annulus
skip_distance = 1        
shaft_ratio = 1/8         # Reduced from 1/5 to help thicken the gap

# TARGET ANNULAR THICKNESS (Back-solving for pressure from here)
target_annular_in = 0.015 
annular_thk_phys = target_annular_in * in_to_m

target_LMR_min, target_LMR_max = 1.0, 5.0
target_TMR_min, target_TMR_max = 0.5, 5 # Relaxed slightly for initial solve

# --- CALCULATIONS ---
film_percent = 0.05                   
m_dot_fuel_pint = m_dot_fuel * (1 - film_percent)

shaft_dia = d_c * shaft_ratio
shaft_rad = shaft_dia / 2

# N2O Properties (Simplified interpolation)
n2o = pd.read_excel(r"N20 Densities.xlsx")
ox_temp = 70
ox_pressure = np.interp(ox_temp, pd.to_numeric(n2o.iloc[:,0], errors="coerce"), pd.to_numeric(n2o.iloc[:,1], errors="coerce")) * 1000
ox_rho = np.interp(ox_temp, pd.to_numeric(n2o.iloc[:,0], errors="coerce"), pd.to_numeric(n2o.iloc[:,2], errors="coerce"))
fuel_rho = 789

# --- BACK-SOLVE FOR FUEL DELTA P ---
# 1. Calculate Physical Area of the gap
A_fuel_phys = np.pi * ((shaft_rad + annular_thk_phys)**2 - shaft_rad**2)

# 2. Calculate Effective Area (A_eff = Cd * A_phys)
A_fuel_eff = A_fuel_phys * Cd_annulus

# 3. Back-solve for required Delta P: delta_P = (m_dot / (Cd * A_phys))**2 / (2 * rho)
# Note: Since A_eff = Cd * A_phys, this simplifies to (m_dot / A_eff)**2 / (2 * rho)
delta_P_fuel_required = (m_dot_fuel_pint / A_fuel_eff)**2 / (2 * fuel_rho)

# 4. Calculate Fuel Velocity (Momentum Velocity)
vel_fuel = m_dot_fuel_pint / (A_fuel_eff * fuel_rho)

# --- OXIDIZER PRESSURE DROP ---
ox_feed_loss = 135.7 * psi_to_pa
inlet_P_ox = ox_pressure - ox_feed_loss
delta_P_ox = inlet_P_ox - Pc

# --- DRILL BITS ---
drills_list = pd.read_excel(r"Drill_Bits.xlsx")
drills = np.sort(pd.to_numeric(drills_list["Decimal Value (mm)"], errors="coerce").dropna() * 0.001)

results = []

for num_holes in range(10, 120, 2):
    area_ox_req = m_dot_ox / (discharge_coef * np.sqrt(2 * ox_rho * delta_P_ox))
    hole_diameter_req = 2 * np.sqrt(area_ox_req / (np.pi * num_holes))

    idx = np.argmin(np.abs(hole_diameter_req - drills))
    act_dia_ox = drills[idx]
    act_A_ox_phys = num_holes * np.pi * (act_dia_ox / 2)**2
    
    # Effective Ox Velocity
    vel_ox = m_dot_ox / (act_A_ox_phys * discharge_coef * ox_rho)
    
    # Momentum Ratios
    TMR = (m_dot_ox * vel_ox) / (m_dot_fuel_pint * vel_fuel)
    BF = (num_holes * act_dia_ox) / (np.pi * shaft_dia)
    LMR = TMR / BF
    
    act_delta_P_ox = (m_dot_ox / (discharge_coef * act_A_ox_phys))**2 / (2 * ox_rho)

    if (target_TMR_min <= TMR <= target_TMR_max) and (target_LMR_min <= LMR <= target_LMR_max):
        spray_angle = np.degrees(2 * 0.7 * np.arctan(2 * LMR))
        
        results.append({
            "num_holes": num_holes,
            "hole_dia_mm": act_dia_ox * 1000,
            "annular_thk_in": annular_thk_phys / in_to_m,
            "LMR": LMR,
            "TMR": TMR,
            "BF": BF,
            "spray_deg": spray_angle,
            "vel_ox": vel_ox,
            "vel_fuel": vel_fuel,
            "delta_P_fuel_psi": delta_P_fuel_required / psi_to_pa,
            "delta_P_ox_psi": act_delta_P_ox / psi_to_pa,
            "req_tank_P_fuel_psi": (Pc + delta_P_fuel_required + 65.6*psi_to_pa) / psi_to_pa
        })

# --- OUTPUT ---
if results:
    results_df = pd.DataFrame(results)
    print(f"Sizing for Target Gap: {target_annular_in} in")
    print(f"Required Fuel Delta P: {delta_P_fuel_required/psi_to_pa:.2f} psi")
    print(results_df.head(5).round(4))
else:
    print("No valid configurations. Try changing shaft_ratio or target_annular_in.")