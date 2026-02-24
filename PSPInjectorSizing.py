#    Greer McDevitt    #
# Last Updated 7/19/25 #

# https://cearun.grc.nasa.gov/cgi-bin/CEARUN/donecea3.cgi
# https://www.eucass.eu/doi/EUCASS2017-474.pdf Source used for L* and Pc value

#INPUTS
import numpy as np
import pandas as pd
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt

mode = "Hotfire"

if mode == "Hotfire": # Hotfire Input Values
    mass_flow_total = 3.1437919  # lbm/s
    OF_ratio = 1.58
    chamber_pressure = 200       # psi
    ox_temp = 90                 # K (LOX)
    fuel_temp = 298
if mode == "Water":  # Water Input values
    mass_flow_total = 1.06       # lbm/s
    OF_ratio = 1
    chamber_pressure = 14.7      # psi
    ox_temp = 293                # K (water proxy for LOX)
    fuel_temp = 293              # K (water proxy for ethanol)
    
# User Inputs
chamber_diameter = 3.975     # in
discharge_coef = 0.65
skip_distance = 1     # ratio of skip length(distance from annular to radial flow) to pintle diameter.
shaft_ratio = 1/5     # ratio used with BZ1 and BZB

target_LMR_min = 1.0  # Minimum LMR (Flow characteristics of a pintle injector element https://www.sciencedirect.com/science/article/pii/S0094576518309883#fd2)
target_LMR_max = 3.0  # Maximum LMR (Range between 1.5 and 3.0 recommended for best atomization and Wide, uniform spray pattern)

target_TMR_min = 0.9 # keep around this range to have Efficient shear mixing and optimize C*
target_TMR_max = 1.5

# Conversion Factors and Constants
lbm_to_kg = 0.453592
psi_to_pa = 6894.76
in_to_m = 0.0254

ox_temp = 90 #293          # Liquid Oxygen Temp K
fuel_temp = 298 #293       # Ethanol Temp K

# Unit Conversions
mass_flow_total *= lbm_to_kg
chamber_pressure *= psi_to_pa
chamber_diameter *= in_to_m

# Pressure Drops
if mode == "Hotfire":
    ox_pressure_drop = chamber_pressure * 0.20
    fuel_pressure_drop = chamber_pressure * 0.20
elif mode == "Water":  # Water test
    min_drop = 40 * psi_to_pa  # 40 psi minimum
    ox_pressure_drop = max(chamber_pressure * 0.80, min_drop)
    fuel_pressure_drop = max(chamber_pressure * 0.80, min_drop)

ox_inlet_pressure = chamber_pressure + ox_pressure_drop
fuel_inlet_pressure = chamber_pressure + fuel_pressure_drop

# Fluid Properties
ox_density = CP.PropsSI('D', 'T', ox_temp, 'P', ox_inlet_pressure, 'Oxygen')               #[kg/m^3] density of LOX
fuel_density = CP.PropsSI('D', 'T', fuel_temp, 'P', fuel_inlet_pressure, 'Ethanol')          #[kg/m^3] density of FUEL

# Pintle Geometry
shaft_diameter = chamber_diameter * shaft_ratio  #in
shaft_radius = shaft_diameter / 2                #in
skip_length = shaft_diameter/in_to_m * skip_distance
print(skip_length)
# Mass Flow Calculations
fc_fuel_ratio = 0.05  # 5% film cooling
fuel_flow_total = mass_flow_total / (OF_ratio + 1)
ox_flow_rate = fuel_flow_total * OF_ratio
fuel_flow_rate = fuel_flow_total * (1 - fc_fuel_ratio)

# Load Drill Bit Sizes
drill_bits_df = pd.read_csv(r"Injectors\Heatsink_Injector\bidc_bit_sizes.csv")
drill_bits = np.sort(pd.to_numeric(drill_bits_df["Bits in millimeters"], errors='coerce').dropna() * 0.001) # Convert mm to m

results = []

# Main Optimization Loop
for num_holes in range(10, 100, 2): #increment by 2

    # Calculate theoretical hole diameter
    area_OX = ox_flow_rate / (discharge_coef * np.sqrt(2 * ox_density * ox_pressure_drop)) # Standard Orifice Equation
    hole_diameter = 2 * np.sqrt(area_OX / (np.pi * num_holes))                             # Diameter derived from area of a circle times the number of holes

    # Snap to nearest drill bit size
    idx = np.argmin(np.abs(hole_diameter - drill_bits))
    real_diam_OX = drill_bits[idx]
    real_area_OX = num_holes * np.pi * (real_diam_OX / 2)**2
    
    
    # Calculate velocities
    vel_OX = ox_flow_rate / (real_area_OX * ox_density) # continuity equation
    
    annular_thk = (np.pi * ox_density * real_diam_OX) / (4 * fuel_density * OF_ratio**2) # NASA SP-8089 Annular Thickness Calculation
    area_FUEL = np.pi * ((shaft_radius + annular_thk)**2 - shaft_radius**2)
    vel_FUEL = fuel_flow_rate / (area_FUEL * fuel_density)
    
    # Momentum Ratios
    TMR = (ox_flow_rate * vel_OX) / (fuel_flow_rate * vel_FUEL)
    blockage_factor = (num_holes * real_diam_OX) / (np.pi * shaft_diameter)
    LMR = TMR / blockage_factor
    
    # Only store configurations within target LMR range
    if (target_TMR_min <= TMR <= target_TMR_max) and (target_LMR_min <= LMR <= target_LMR_max):
        spray_angle = np.degrees(2 * 0.7 * np.arctan(2 * LMR))  # From Flow characteristics of a pintle injector element Eq(3)
        
        results.append({
            'num_holes': num_holes,
            'hole_diameter_in': real_diam_OX / in_to_m,  # Convert back to inches
            'hole_diameter_mm': real_diam_OX * 1000,     # Diameter in mm (convert m to mm)
            'annular_thickness': annular_thk / in_to_m,
            'LMR': LMR,
            'TMR': TMR,
            'blockage_factor': blockage_factor,
            'spray_angle_degrees': spray_angle,
            'vel_OX': vel_OX,
            'vel_FUEL': vel_FUEL,
            'area_OX_in': real_area_OX/ in_to_m**2,  # Convert area to inches^2
            'area_FUEL_in': area_FUEL/ in_to_m**2,  # Convert area to inches^2
        })
# Output Results
if results:
    target_TMR_mid = (target_TMR_min + target_TMR_max) / 2
    results.sort(key=lambda x: abs(x['LMR'] - target_TMR_mid))
    
    # Create DataFrame and save to Excel
    results_df = pd.DataFrame(results)
    results_df.to_excel('optimized_injector_configs.xlsx', index=False)
    print(f"Found {len(results)} valid configurations. Top 3:")
    top3 = results_df.head(3).round(3).reset_index(drop=True)
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
    plt.plot(results_df['num_holes'], results_df['TMR'], 'bo-')
    plt.xlabel('Number of Holes')
    plt.ylabel('TMR')
    plt.title('TMR vs Hole Count')
    plt.grid(True)
    
    # LMR vs Hole Count
    plt.subplot(1, 2, 2)
    plt.plot(results_df['num_holes'], results_df['LMR'], 'ro-')
    plt.xlabel('Number of Holes')
    plt.ylabel('LMR')
    plt.title('LMR vs Hole Count')
    plt.grid(True)
    
    # Add horizontal lines showing target LMR range
    plt.axhline(y=target_LMR_min, color='gray', linestyle='--')
    plt.axhline(y=target_LMR_max, color='gray', linestyle='--')
    
    plt.tight_layout()
    
    # Save and show plot
    plt.savefig('hole_count_vs_momentum_ratios.png', dpi=300)
    plt.show()
else:
    if plot == 1:
        print("Not enough data points to generate meaningful plots")