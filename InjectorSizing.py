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
m_dot_ox = m_dot_fuel * OF            # Oxidizer mass flow [kg/s]
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
ox_rho = CP.PropsSI ("D", "T", ox_temp, "P", inlet_P, "Oxygen")
fuel_rho = CP.PropsSI("D", "T", fuel_temp, "P", inlet_P, "Ethanol")

# Available Drill Bit Sizes
drills_list = pd.read_csv(r"drill_sizes.csv")
drills = np.sort(pd.to_numeric(drills_list["Bits in millimeters"], errors='coerce').dropna() * 0.001) # Convert mm to m
results = []





