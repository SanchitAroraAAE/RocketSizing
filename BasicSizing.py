import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass


# DATA STRUCTURE
@dataclass
class RocketSizing:
    m_dot_total: float
    m_dot_fuel: float
    m_dot_ox: float
    d_c: float
    d_t: float
    d_e: float
    L_c: float
    L_n: float
    CR: float
    ER: float
    theta_n: float
    theta_e: float
    OF: float
    Pc: float
    isp: float
    thrust: float


def BasicSizing(mode="Hotfire"):
    """
    Calculates fundamental rocket engine dimensions based on target thrust
    and chamber pressure using CEA-derived values.
    """
    
    # CONVERSIONS
    lbf_to_N = 4.44822162
    psi_to_pa = 6894.76
    in_to_m = 0.0254


    # INPUTS & CONSTANTS
    thrust = 400 * lbf_to_N
    d_c = 3.25 * in_to_m        # Chamber ID
    L_star = 60 * in_to_m       # Characteristic Length (L*)
    percent_bell = 0.8          # 80% Rao nozzle
    
    # Efficiency Factors
    eta_cstar = 0.85            # Combustion efficiency
    eta_cf = 0.95               # Nozzle/Thrust coefficient efficiency
    
    # CEA Values (Nitrous/Ethanol @ 300psi)
    c_ideal = 1649.9            # Effective exhaust velocity [m/s]
    c_star_ideal = 1154.1       # Characteristic velocity [m/s]
    ER = 3.9821                 # Expansion Ratio (Ae/At)


    # MODE SELECTION
    if mode == "Hotfire":
        OF = 3.0
        Pc = 300 * psi_to_pa
    else:  # Waterflow
        OF = 1.0
        Pc = 14.7 * psi_to_pa


    # PERFORMANCE CALCULATIONS 
    # Actual performance adjusted for efficiencies
    c_actual = c_ideal * eta_cstar * eta_cf
    c_star_actual = c_star_ideal * eta_cstar
    isp = c_actual / 9.80665

    m_dot_total = thrust / c_actual
    m_dot_fuel = m_dot_total / (1 + OF)
    m_dot_ox = m_dot_fuel * OF


    # NOZZLE GEOMETRY
    A_t = (c_star_actual * m_dot_total) / Pc
    d_t = 2 * np.sqrt(A_t / np.pi)
    
    A_e = A_t * ER
    d_e = 2 * np.sqrt(A_e / np.pi)
    
    # Initial estimate of nozzle length (L_n)
    L_n = percent_bell * (np.sqrt(ER) - 1) * (d_t / 2) / (np.tan(np.deg2rad(15)))


    # CHAMBER GEOMETRY
    A_c = np.pi * (d_c / 2)**2
    CR = A_c / A_t              # Contraction Ratio
    V_c = L_star * A_t          # Chamber Volume based on L*
    L_c = V_c / A_c             # Length from injector face to throat

    # RAO BELL ANGLE EXTRAPOLATION
    # 80% Rao Nozzle empirical table
    eratio_table = [4, 5, 10, 20, 30, 40, 50, 100]
    theta_n_table = [21.5, 23.0, 26.3, 28.8, 30.0, 31.0, 31.5, 33.5]
    theta_e_table = [14.0, 13.0, 11.0, 9.0, 8.5, 8.0, 7.5, 7.0]

    if ER < eratio_table[0]:
        # Linear extrapolation for small expansion ratios
        slope_n = (theta_n_table[1] - theta_n_table[0]) / (eratio_table[1] - eratio_table[0])
        slope_e = (theta_e_table[1] - theta_e_table[0]) / (eratio_table[1] - eratio_table[0])
        
        theta_n = theta_n_table[0] + slope_n * (ER - eratio_table[0])
        theta_e = theta_e_table[0] + slope_e * (ER - eratio_table[0])
    else:
        theta_n = np.interp(ER, eratio_table, theta_n_table)
        theta_e = np.interp(ER, eratio_table, theta_e_table)


    # RESULTS
    return RocketSizing(
        m_dot_total, m_dot_fuel, m_dot_ox, d_c, d_t, d_e, L_c, L_n, CR, ER, theta_n, theta_e, OF, Pc, isp, thrust)