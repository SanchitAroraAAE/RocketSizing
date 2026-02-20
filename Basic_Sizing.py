import numpy as np
import matplotlib.pyplot as plt

# https://cearun.grc.nasa.gov/cgi-bin/CEARUN/donecea3.cgi

#INPUTS
thrust = 1556.87757 #[N], 350lbf 
Pc = 2413165.048 #chamber pressure[Pa], 350 psia
OF = 3 #oxidize-fuel ratio

#CALCULATIONS
c = 1687.5 #effective exhaust velocity [m/s]
c_star = 1156.7 #characteristic exhaust velocity [m/s]
exp_Ratio = 4.4659 #expansion ratio Ae/At
isp = c/9.81 #specific impulse [1/s]

m_dot_total = thrust/c #total mass flow [kg/s]
m_dot_fuel = m_dot_total/(1+OF) #fuel mass flow [kg/s]
m_dot_ox = m_dot_fuel*OF #oxidizer mass flow [kg/s]

A_t = c_star*m_dot_total/Pc #throat area [m^2]
d_t = (2 * np.sqrt(A_t/np.pi))*100 #throat diameter [cm]
A_e = A_t * exp_Ratio #exit area [m^2]
d_e = (2 * np.sqrt(A_e/np.pi))*100 #throat diameter [c]

#PRINT OUTPUTS
print(f"Total mass flow: {m_dot_total: .3f} kg/s")
print(f"Fuel mass flow: {m_dot_fuel: .3f} kg/s")
print(f"Oxidizer mass flow: {m_dot_ox: .3f} kg/s")

print(f"Throat diameter: {d_t: .3f} cm")
print(f"Exit diameter: {d_e: .3f} cm")

