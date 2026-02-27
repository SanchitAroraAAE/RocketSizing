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
    m_dot_total = 0.944   # kg/s
    OF = 3
    Pc = 2068427.184          # Chamber pressure [Pa], 300 psia
    ox_temp = 305.37          # NOs temp [K]
    fuel_temp = 3298          # E98 temp [k]
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
ox_rho = CP.PropsSI ("D", "T", ox_temp, "P", inlet_P, "oxygen")
fuel_rho = CP.PropsSI("D", "T", fuel_temp, "P", inlet_P, "Ethanol")

# Available Drill Bit Sizes
drills_list = pd.read_csv(r"Drill_Bits.csv")
drills = np.sort(pd.to_numeric(drills_list["Decimal Value (mm)"], errors="coerce").dropna() * 0.001) # Convert mm to m
results = []

#Optimization Loop (from flowchart on PSP confluence)
for num_holes in range(10, 120, 2): # Needs to have atleast 10 holes. Increment by 2 for efficiency

    # Calculate theoretical hole diameter
    area_ox = m_dot_ox / (discharge_coef * np.sqrt(2 * ox_rho * delta_P))     # Standard Orifice Equation
    hole_diameter = 2 * np.sqrt(area_ox / (np.pi * num_holes))                # Area of a circle times the number of holes needs to be total ox area

    # Find nearest drill size
    idx = np.argmin(np.abs(hole_diameter - drills))
    act_dia_ox = drills[idx]
    act_A_ox = num_holes * np.pi * (act_dia_ox / 2)**2
    
    
    # Calc velocities
    vel_ox = m_dot_ox / (act_A_ox*ox_rho)                                    # Find exit velocity of oxdizer
    
    annular_thk = (np.pi*ox_rho*act_dia_ox) / (4*fuel_rho*(OF**2))           # Eqt 1.9 from PSP Injector Design and Analysis page (which I think got it from NASA SP-8089)
    A_fuel = np.pi * ((shaft_rad + annular_thk)**2 - shaft_rad**2)           
    vel_fuel = m_dot_fuel_pint / (A_fuel * fuel_rho)                         
    
    # Momentum Ratios
    TMR = (m_dot_ox * vel_ox) / (m_dot_fuel_pint * vel_fuel)                 # Eqt. 1.7 from PSP page
    BF = (num_holes * act_dia_ox) / (np.pi * shaft_dia)         # Eqt. 1.11 from PSP page
    LMR = TMR / BF                                            # Eqt. 1.13 from PSP page
    
    # Check how many rows are needed
    num_rows = 1
    if (BF>1):
        num_rows = 2

    # Only store configurations within target LMR range
    if (target_TMR_min <= TMR <= target_TMR_max) and (target_LMR_min <= LMR <= target_LMR_max):
        spray_angle = np.degrees(2 * 0.7 * np.arctan(2 * LMR))               # Eqt. 1.14 from PSP page (from Flow characteristics of a pintle injector element Eq(3))
        
        results.append({
            "num_holes": num_holes,
            "num_rows": num_rows,
            "hole_diam_in": act_dia_ox / in_to_m,    # Convert back to inches
            "hole_dia_mm": act_dia_ox * 1000,        # Diameter in mm (convert m to mm)
            "annular_thk": annular_thk / in_to_m,
            "LMR": LMR,
            "TMR": TMR,
            "blockage_factor": BF,
            "spray_angle_deg": spray_angle,
            "vel_ox": vel_ox,
            "vel_fuel": vel_fuel,
            "area_ox_in": act_A_ox/ in_to_m**2,      # Convert area to inches^2
            "area_fuel_in": A_fuel/ in_to_m**2,      # Convert area to inches^2
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
plot = 0 # Set to 1 to enable plotting
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
    plt.savefig("hole_count_vs_momentum_ratios.png", dpi=300)
    plt.show()
else:
    if plot == 1:
        print("Not enough data points to generate meaningful plots")





