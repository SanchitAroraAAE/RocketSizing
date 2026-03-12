import numpy as np


# REFERENCES
# Major Losses: https://en.wikipedia.org/wiki/Darcy%E2%80%93Weisbach_equation
# Pressure Drop Calculator: https://www.pressure-drop.online/

def calculate_pressure_drop(m_dot, rho, mu, L_in, D_in, epsilon_in, Cv_total):
    """
    Calculates total pressure drop across a feed system combining major 
    (friction) and minor (valves/fittings) losses.
    
    Parameters:
    m_dot (float): Mass flow rate [kg/s]
    rho (float):   Fluid density [kg/m^3]
    mu (float):    Dynamic viscosity [Pa*s]
    L_in (float):  Pipe length [inches]
    D_in (float):  Pipe internal diameter [inches]
    epsilon_in (float): Surface roughness [inches]
    Cv_total (float): Combined flow coefficient of valves/fittings
    
    Returns:
    float: Total pressure drop [Pa]
    """
    
    # UNIT CONVERSIONS
    in_to_m = 0.0254
    psi_to_pa = 6894.76
    mps3_to_gpm = 15850.3  # m^3/s to Gallons Per Minute
    
    L = L_in * in_to_m
    D = D_in * in_to_m
    epsilon = epsilon_in * in_to_m
    

    # FLOW GEOMETRY & REYNOLDS NUMBER
    area = np.pi * (D / 2)**2
    vel = m_dot / (rho * area)
    reynolds = (rho * vel * D) / mu
    

    # MAJOR LOSS: DARCY-WEISBACH
    # Determine Darcy Friction Factor (f)
    if reynolds < 2300:
        # Laminar Flow
        f = 64 / reynolds
    else:
        # Turbulent Flow: Haaland Equation (approximation for Colebrook-White)
        term_a = (epsilon / D) / 3.7
        term_b = 6.9 / reynolds
        f = (-1.8 * np.log10(term_a**1.11 + term_b))**-2
        
    # Major Pressure Drop [Pa]
    dp_major = f * (L / D) * (0.5 * rho * vel**2)
    

    # MINOR LOSS: FLOW COEFFICIENT (Cv)
    # Governing Equation: dP(psi) = SG * (Q(gpm) / Cv)^2
    q_gpm = (m_dot / rho) * mps3_to_gpm
    sg = rho / 1000.0  # Specific Gravity relative to water
    
    dp_minor_psi = sg * (q_gpm / Cv_total)**2
    dp_minor = dp_minor_psi * psi_to_pa
    
    # RESULTS
    return dp_major + dp_minor