import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
import FeedPressureDrop
from BasicSizing import BasicSizing


# REFERENCES
# Injector Design: https://purdue-space-program.atlassian.net/wiki/spaces/PL/pages/180486437/Injector+Design+and+Analysis
# Phoenix Injector: https://purdue-space-program.atlassian.net/wiki/spaces/PL/pages/1248264194/Phoenix+Injector
# Film Cooling: https://events.iist.ac.in/phd/thesis/SC09D002%20FT.pdf


# NITIALIZATION & INPUTS
mode = "Hotfire"
sizing = BasicSizing(mode)

# Conversion Factors
LBM_TO_KG = 0.453592
PSI_TO_PA = 6894.76
IN_TO_M = 0.0254

# Sizing Parameters
m_dot_total = sizing.m_dot_total  # [kg/s]
m_dot_fuel = sizing.m_dot_fuel    # [kg/s]
m_dot_ox = sizing.m_dot_ox        # [kg/s]
OF = sizing.OF
Pc = sizing.Pc
d_c = sizing.d_c                  # Chamber diameter [m]

# Thermodynamic Conditions
if mode == "Hotfire":
    ox_temp = 70                  # [F] (26.6 C)
    fuel_temp = 298               # [K]
else: # Waterflow
    ox_temp = 293                 # [K]
    fuel_temp = 293               # [K]

# Design Constraints & Targets
discharge_coef = 0.65
Cd_annulus = 0.6
skip_distance = 1                 # Ratio of skip length to pintle diameter
shaft_ratio = 1/5                 # Ratio for shaft diameter calculation
film_percent = 0.05               # 5% film cooling fraction

target_LMR_min, target_LMR_max = 1.0, 3.0
target_TMR_min, target_TMR_max = 0.9, 2.0

# Pressure Loss Constants
piston_loss = 15 * PSI_TO_PA
ox_feed_loss = 135.7 * PSI_TO_PA
fuel_feed_loss = 65.6 * PSI_TO_PA


# PRE-LOOP CALCULATIONS

# Mass Flow & Geometry
m_dot_fuel_pint = m_dot_fuel * (1 - film_percent)
shaft_dia = d_c * shaft_ratio
shaft_rad = shaft_dia / 2
skip_len = skip_distance * shaft_dia

# Oxidizer Properties (Interpolation)
n2o = pd.read_excel(r"N20 Densities.xlsx")
ox_temps = pd.to_numeric(n2o.iloc[:, 0], errors="coerce").to_numpy()
ox_pressures = pd.to_numeric(n2o.iloc[:, 1], errors="coerce").to_numpy()
ox_rhos = pd.to_numeric(n2o.iloc[:, 2], errors="coerce").to_numpy()

ox_pressure = np.interp(ox_temp, ox_temps, ox_pressures) * 1000 # [Pa]
ox_rho_interp = np.interp(ox_temp, ox_temps, ox_rhos)

if mode == "Hotfire":
    ox_rho = ox_rho_interp
    fuel_rho = 789
else: # Waterflow
    ox_rho = 1000
    fuel_rho = 1000

# Available Drills
drills_list = pd.read_excel(r"Drill_Bits.xlsx")
drills = np.sort(pd.to_numeric(drills_list["Decimal Value (mm)"], errors="coerce").dropna() * 0.001)

# System Pressure Drops
if mode == "Hotfire":
    inlet_P_ox = ox_pressure - ox_feed_loss
    inlet_P_fuel = ox_pressure - fuel_feed_loss
    delta_P_ox = inlet_P_ox - Pc
    delta_P_fuel = inlet_P_fuel - Pc
else:
    min_drop = 40 * PSI_TO_PA
    delta_P_ox = max(Pc * 0.8, min_drop)
    delta_P_fuel = delta_P_ox # Assuming similar for waterflow


# PTIMIZATION LOOP
results = []

for num_holes in range(10, 120, 2):
    # Determine Hole Size based on pressure drop
    area_ox_req = m_dot_ox / (discharge_coef * np.sqrt(2 * ox_rho * delta_P_ox))
    ideal_hole_dia = 2 * np.sqrt(area_ox_req / (np.pi * num_holes))

    # Match to real drill bit
    idx = np.argmin(np.abs(ideal_hole_dia - drills))
    act_dia_ox = drills[idx]
    act_A_ox = num_holes * np.pi * (act_dia_ox / 2)**2
    
    # Kinematics
    vel_ox = m_dot_ox / (act_A_ox * ox_rho)

    # Annulus Geometry (Eq 1.9 from PSP / NASA SP-8089)
    annular_thk = (np.pi * ox_rho * act_dia_ox) / (4 * fuel_rho * (OF**2))
    A_fuel_eff = np.pi * ((shaft_rad + annular_thk)**2 - shaft_rad**2)
    vel_fuel = m_dot_fuel_pint / (A_fuel_eff * fuel_rho)
    
    # Momentum Ratios & Mixing
    TMR = (m_dot_ox * vel_ox) / (m_dot_fuel_pint * vel_fuel)
    BF = (num_holes * act_dia_ox) / (np.pi * shaft_dia) # Blockage Factor
    LMR = TMR / BF
    
    num_rows = 2 if BF > 1 else 1

    # Final Verification
    act_delta_P = (m_dot_ox / (discharge_coef * act_A_ox))**2 / (2 * ox_rho)
    
    # Filtering and Result Collection
    if (target_TMR_min <= TMR <= target_TMR_max) and (target_LMR_min <= LMR <= target_LMR_max):
        spray_angle = np.degrees(2 * 0.7 * np.arctan(2 * LMR))
        
        results.append({
            "num_holes": num_holes,
            "num_rows": num_rows,
            "hole_diam_in": act_dia_ox / IN_TO_M,
            "hole_dia_mm": act_dia_ox * 1000,
            "annular_thk": annular_thk / IN_TO_M,
            "LMR": LMR,
            "TMR": TMR,
            "blockage_factor": BF,
            "spray_angle_deg": spray_angle,
            "vel_ox": vel_ox,
            "vel_fuel": vel_fuel,
            "area_ox_in": act_A_ox / IN_TO_M**2,
            "area_fuel_in": A_fuel_eff / IN_TO_M**2,
            "actual_delta_P_psi": act_delta_P / PSI_TO_PA,
            "delta_P_error_percent": ((act_delta_P / delta_P_ox) - 1) * 100,
        })


# OTPUTS & PLOTTING
# Output Results
if results:
    # Sort by proximity to center of TMR target range
    target_TMR_mid = (target_TMR_min + target_TMR_max) / 2
    results.sort(key=lambda x: abs(x["LMR"] - target_TMR_mid))
    
    # Create DataFrame and save to Excel
    results_df = pd.DataFrame(results)
    #results_df.to_excel("optimized_injector_configs.xlsx", index=False)
    print(f"Found {len(results)} valid configurations. Top 3:")
    top3 = results_df.head(10).round(5).reset_index(drop=True)
    top3.index += 1
    print(top3)
else:
    print("No valid configurations found. Try relaxing constraints.")

# Plot Relationship between Number of Holes and LMR
plot_enabled = True # Set to 1 to enable plotting
if len(results) > 1 and plot_enabled:
    # Create figure with two subplots
    plt.figure(figsize=(12, 5))

    # Plot TMR vs Hole Count
    plt.subplot(1, 2, 1)
    plt.plot(results_df["num_holes"], results_df["TMR"], "bo-")
    plt.xlabel("Number of Holes")
    plt.ylabel("TMR")
    plt.title("TMR vs Hole Count")
    plt.grid(True)
    
    # Plot LMR vs Hole Count
    plt.subplot(1, 2, 2)
    plt.plot(results_df["num_holes"], results_df["LMR"], "ro-")
    plt.xlabel("Number of Holes")
    plt.ylabel("LMR")
    plt.title("LMR vs Hole Count")
    plt.grid(True)
    
    # Add horizontal lines showing target LMR range
    plt.axhline(y=target_LMR_min, color="gray", linestyle="--")
    plt.axhline(y=target_LMR_max, color="gray", linestyle="--")
    
    plt.tight_layout()
    
    # Save and show plot
    # plt.savefig("hole_count_vs_momentum_ratios.png", dpi=300)
    plt.show()
else:
    if plot_enabled == 1:
        print("Not enough data points to generate meaningful plots")





