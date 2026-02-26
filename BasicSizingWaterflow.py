import numpy as np
import matplotlib.pyplot as plt

# https://cearun.grc.nasa.gov/cgi-bin/CEARUN/donecea3.cgi
# https://www.eucass.eu/doi/EUCASS2017-474.pdf Source used for L* and Pc value

#INPUTS
thrust = 1556.87757     # [N], 350lbf 
Pc = 101352.932016      # Chamber pressure[Pa], 14.7 psia
OF = 1                  # Oxidize-fuel ratio
d_c = 0.08255           # Chamber dia. [m] 3.25"
L_star = 1.524          # Characteristic length [m], 60"
percent_bell = 0.8      # Percent rao nozzle

#GENERAL CALCULATIONS
c = 1344.6          # Effective exhaust velocity [m/s]
c_star = 2081.0     # Characteristic exhaust velocity [m/s]
ER = 1              # Expansion ratio Ae/At
isp = c/9.81        # Specific impulse [1/s]

m_dot_total = thrust/c              # Total mass flow [kg/s]
m_dot_fuel = m_dot_total/(1+OF)     # Fuel mass flow [kg/s]
m_dot_ox = m_dot_fuel*OF            # Oxidizer mass flow [kg/s]

#NOZZLE CALCULATIONS
A_t = c_star*m_dot_total/Pc      # Throat area [m^2]
d_t = 2 * np.sqrt(A_t/np.pi)     # Throat diameter [m]
A_e = A_t * ER                   # Exit area [m^2]
d_e = 2 * np.sqrt(A_e/np.pi)     # Throat diameter [m]
L_n = percent_bell * (np.sqrt(ER)-1)*(d_t/2)/(np.tan(np.deg2rad(15)))

#CHAMBER CALCULATIONS
A_c = np.pi * (d_c/2)**2     # Chamber area [m^2]
CR = A_c/A_t                 # Contraction ratio Ac/At
V_c = L_star * A_t           # Chamber volume [m^3]
L_c = V_c/A_c                # Chamber length (from injector face to throat) [m]

#THETA CALCS
# 80% Rao Nozzle empirical data
eratio     = [4,    5,    10,   20,   30,   40,   50,   100]
theta_n_80 = [21.5, 23.0, 26.3, 28.8, 30.0, 31.0, 31.5, 33.5]
theta_e_80 = [14.0, 13.0, 11.0, 9.0,  8.5,  8.0,  7.5,  7.0]

# Manual linear extrapolation for ER < 4
if ER < eratio[0]:
    # Calculate slope (m) between the first two points
    m_n = (theta_n_80[1] - theta_n_80[0]) / (eratio[1] - eratio[0])
    m_e = (theta_e_80[1] - theta_e_80[0]) / (eratio[1] - eratio[0])
    
    # y = y1 + m(x - x1)
    theta_n_deg = theta_n_80[0] + m_n * (ER - eratio[0])
    theta_e_deg = theta_e_80[0] + m_e * (ER - eratio[0])
else:
    theta_n_deg = np.interp(ER, eratio, theta_n_80)
    theta_e_deg = np.interp(ER, eratio, theta_e_80)

theta_n = (theta_n_deg)     # Theta N value [deg]
theta_e = (theta_e_deg)     # Theta E value [deg]

#PRINT OUTPUTS
print(f"Total mass flow: {m_dot_total: .3f} kg/s")
print(f"Fuel mass flow: {m_dot_fuel: .3f} kg/s")
print(f"Oxidizer mass flow: {m_dot_ox: .3f} kg/s\n")

print(f"Chamber diameter: {d_c: .6f} m")
print(f"Throat diameter: {d_t: .6f} m")
print(f"Exit diameter: {d_e: .6f} m")
print(f"Chamber length: {L_c: .6f} m")
print(f"Nozzle length; {L_n: .6f} m\n")

print(f"Contraction Ratio: {CR: .3f}")
print(f"Theta N : {theta_n: .6f} deg")
print(f"Theta E: {theta_e: .6f} deg\n")



