# This code provides a function you can use to find the pressure drop across a feed system in order to find the 
# pressure at the injector by subtracting this dP from tank pressures

# Major Losses to friction:
# Uses the Darcy-Weisbach equation, using the Weisbach friction factor
# dP [Pa] = f * (L/D) * ((rho*v^2) / 2)

# Minor losses from valves/fittings:
# Derived from the Bernoulli equation where pressure drop (dP) is proportional to the square of the velocity:
# dP = K * (1/2 * rho * v^2).
#
# The Flow Coefficient (Cv) is an industry-standard value representing the flow of 60°F water in Gallons Per Minute (GPM) 
# that results in a 1 psi pressure drop across the component.
#
# Governing Equation: dP = SG * (Q / Cv)^2
# dP: Pressure drop across the component [psi]
# SG: Specific Gravity (density_fluid / density_water)
# Q:  Volumetric flow rate [gallons per minute (GPM)]
# Cv: Flow Coefficient of the valve/fitting
#
# If a manufacturer provides a dimensionless loss coefficient (K):
# Cv = (D^2 / sqrt(K)) * 29.84  [D = internal diameter in inches]

# Useful links:
# https://en.wikipedia.org/wiki/Darcy%E2%80%93Weisbach_equation
# https://www.pressure-drop.online/

# Imports:
import numpy as np


# Variable Definitions:
# m_dot = mass flow rate of propellant
# rho = density of propellant in tank
# mu = dynamic viscosity of the fluid
# L = length of the run line
# D = inner diameter of the run line
# epsilon = represents the surface roughness of the inside of the pipe. Can find online or in HalfCat Sim
# Cv_total = flow coefficient

def calculate_pressure_drop(m_dot, rho, mu, L, D, epsilon, Cv_total):
    # Calculate Flow Geometry
    area = np.pi * (D / 2)**2
    vel = m_dot / (rho * area)
    
    # Find Reynolds Number
    Re = (rho * vel * D) / mu
    
    # Find Friction Factor (f) - Haaland Equation
    if Re < 2300:
        f = 64 / Re  # Laminar
    else:
        # Turbulent approximation
        f = (-1.8 * np.log10((epsilon/D / 3.7)**1.11 + 6.9/Re))**-2
        
    # Calculate Major Loss (Pa)
    dP_line = f * (L / D) * (0.5 * rho * vel**2)
    
    # Calculate Minor Loss via Cv (Converted to SI units or calculated in Imperial)
    # Note: HalfCatSim uses Imperial for Cv. 1 gpm = 6.309e-5 m^3/s
    Q_gpm = (m_dot / rho) * 15850.3  # Convert kg/s to GPM
    SG = rho / 1000
    dP_valves_psi = SG * (Q_gpm / Cv_total)**2
    dP_valves = dP_valves_psi * 6894.76  # Convert psi to Pa
    
    return dP_line + dP_valves