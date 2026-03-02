import numpy as np
import pandas as pd
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt

#https://purdue-space-program.atlassian.net/wiki/spaces/PL/pages/180486437/Injector+Design+and+Analysis
#https://purdue-space-program.atlassian.net/wiki/spaces/PL/pages/1248264194/Phoenix+Injector
#https://events.iist.ac.in/phd/thesis/SC09D002%20FT.pdf Film cooling study (in addition to basics from NASA SP-125)

# INPUTS
# Mode selection
mode = "Hotfire"

if mode == "Hotfire": # Hotfire Input Values
    m_dot_total = 1.419   # kg/s
    OF = 3
    Pc = 2068427.184         # Chamber pressure [Pa], 300 psia
    ox_temp = 308            # NOs temp [K], 0deg C
    fuel_temp = 293          # E98 temp [k], 20deg C
if mode == "Water": # Water Input values
    m_dot_total = 1.158       # kg/s
    OF = 1
    Pc = 101352.93201599999   # Chamber pressre [Pa], 14.7psia
    ox_temp = 293             # Water temp [K] (water replacement for NOs)
    fuel_temp = 293           # Water temp [K] (water replacement for E98)

# User Inputs
d_c = 0.08255            # Chamber diameter [m], 3.25"
discharge_coef = 0.65
skip_distance = 1        # Ratio of skip length (distance from annular to radial flow) to pintle diameter.
shaft_ratio = 1/5        # Ratio used with BZ1 and BZB

target_LMR_min = 1.0     # Minimum LMR (Flow characteristics of a pintle injector element https://www.sciencedirect.com/science/article/pii/S0094576518309883#fd2)
target_LMR_max = 3.0     # Maximum LMR (Range between 1.5 and 3.0 recommended for best atomization and Wide, uniform spray pattern)

target_TMR_min = 0.9     # keep around this range to have efficient shear mixing and optimize C*
target_TMR_max = 1.5

# Conversion Factors and Constants
lbm_to_kg = 0.453592
psi_to_pa = 6894.76
in_to_m = 0.0254

# Mass Flow Calcs
film_percent = 0.05                   # 5% film cooling (from phoenix)
m_dot_fuel = m_dot_total / (1+OF)     # Fuel mass flow [kg/s]
m_dot_ox = m_dot_fuel * OF            # oxidizer mass flow [kg/s]
m_dot_fuel_pint = m_dot_fuel * (1-film_percent)

# Pintle Geo. Calcs
shaft_dia = d_c * shaft_ratio
shaft_rad = shaft_dia /2
skip_len = skip_distance * shaft_dia
print(f"Skip length: {skip_len}m")

# Stiffness/Pressure Drops
if mode =="Hotfire":
    delta_P= Pc * 0.2         # 20% Is standard value in industry
elif mode == "Water":
    min_drop = 40 * psi_to_pa # 40psi min
    delta_P = max(Pc * 0.8, min_drop)

inlet_P = Pc + delta_P    # Required injector inlet pressure [Pa]

# Fluid Properties
ox_rho = CP.PropsSI ("D", "T", ox_temp, "P", inlet_P, "NitrousOxide")
fuel_rho = CP.PropsSI("D", "T", fuel_temp, "P", inlet_P, "Ethanol")
print(ox_rho)
print(fuel_rho)

# Available Drill Bit Sizes
drills_list = pd.read_csv(r"Drill_Bits.csv")
drills = np.sort(pd.to_numeric(drills_list["Decimal Value (mm)"], errors="coerce").dropna() * 0.001) # Convert mm to m
results = []

#Optimization Loop (from flowchart on PSP confluence)
# --- MODIFIED OPTIMIZATION LOOP (FUEL CENTERED) ---
for num_holes in range(10, 120, 2): 

    # 1. Calculate theoretical FUEL hole diameter (Holes are now FUEL)
    area_fuel_holes = m_dot_fuel_pint / (discharge_coef * np.sqrt(2 * fuel_rho * delta_P))
    hole_diameter = 2 * np.sqrt(area_fuel_holes / (np.pi * num_holes))

    # Find nearest drill size
    idx = np.argmin(np.abs(hole_diameter - drills))
    act_dia_fuel = drills[idx]
    act_A_fuel = num_holes * np.pi * (act_dia_fuel / 2)**2
    
    # Calc Velocity (Fuel)
    vel_fuel = m_dot_fuel_pint / (act_A_fuel * fuel_rho)
    
    # 2. Calc Annulus (Annulus is now OXIDIZER)
    # Using the same ratio logic: (Area_ox / Area_fuel) ~ (rho_fuel / rho_ox) * (1/OF)
    # Re-deriving annular thickness for OX in the annulus:
    annular_thk = (np.pi * fuel_rho * act_dia_fuel * OF**2) / (4 * ox_rho) # Flipped OF and rho
    A_ox = np.pi * ((shaft_rad + annular_thk)**2 - shaft_rad**2)
    
    vel_ox = m_dot_ox / (A_ox * ox_rho)
    
    # 3. Momentum Ratios (TMR = Radial Momentum / Axial Momentum)
    # In Fuel-Centered: Radial = Fuel, Axial = Ox
    TMR = (m_dot_fuel_pint * vel_fuel) / (m_dot_ox * vel_ox) 
    
    # Blockage Factor (BF) remains based on holes vs circumference
    BF = (num_holes * act_dia_fuel) / (np.pi * shaft_dia)
    LMR = TMR / BF 
    
    num_rows = 1 if BF <= 1 else 2

    # Check against targets
    if (target_TMR_min <= TMR <= target_TMR_max) and (target_LMR_min <= LMR <= target_LMR_max):
        spray_angle = np.degrees(2 * 0.7 * np.arctan(2 * LMR))
        
        results.append({
            "num_holes": num_holes,
            "num_rows": num_rows,
            "hole_dia_mm": act_dia_fuel * 1000,
            "annular_thk_in": annular_thk / in_to_m,
            "LMR": LMR,
            "TMR": TMR,
            "blockage_factor": BF,
            "spray_angle_deg": spray_angle,
            "vel_ox": vel_ox,
            "vel_fuel": vel_fuel,
            "area_ox_in": A_ox / in_to_m**2,
            "area_fuel_in": act_A_fuel / in_to_m**2,
        })

# Output Results
if results:
    target_TMR_mid = (target_TMR_min + target_TMR_max) / 2
    results.sort(key=lambda x: abs(x["LMR"] - target_TMR_mid))
    
    # Create DataFrame and save to Excel
    results_df = pd.DataFrame(results)
    #results_df.to_excel("optimized_injector_configs.xlsx", index=False)
    print(f"Found {len(results)} valid configurations. Top 3:")
    top3 = results_df.head(4).round(3).reset_index(drop=True)
    top3.index += 1
    print(top3)
else:
    print("No valid configurations found. Try relaxing constraints.")

# Plot Relationship between Number of Holes and LMR
plot = 1 # Set to 1 to enable plotting
if results and len(results) > 1 and plot:
    # Create figure with two subplots
    plt.figure(figsize=(12, 5))
    
    # TMR vs Hole Count
    plt.subplot(1, 2, 1)
    plt.plot(results_df["num_holes"], results_df["TMR"], "bo-")
    plt.xlabel("Number of Holes")
    plt.ylabel("TMR")
    plt.title("TMR vs Hole Count")
    plt.grid(True)
    
    # LMR vs Hole Count
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
    #plt.savefig("hole_count_vs_momentum_ratios.png", dpi=300)
    plt.show()
else:
    if plot == 1:
        print("Not enough data points to generate meaningful plots")





