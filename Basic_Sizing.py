import numpy as np
import matplotlib.pyplot as plt

# https://cearun.grc.nasa.gov/cgi-bin/CEARUN/donecea3.cgi

#INPUTS
thrust = 1556.87757 #[N], 350lbf 
Pc = 1723689.3233 #chamber pressure[Pa], 250 psia
OF = 3 #oxidize-fuel ratio
CR = 4 #contraction ratio Ac/At
L_star = 1.0 #characteristic length [m]

#GENERAL CALCULATIONS
c = 1603.9 #effective exhaust velocity [m/s]
c_star = 1151.0 #characteristic exhaust velocity [m/s]
ER = 3.4826 #expansion ratio Ae/At
isp = c/9.81 #specific impulse [1/s]

m_dot_total = thrust/c #total mass flow [kg/s]
m_dot_fuel = m_dot_total/(1+OF) #fuel mass flow [kg/s]
m_dot_ox = m_dot_fuel*OF #oxidizer mass flow [kg/s]

#NOZZLE CALCULATIONS
A_t = c_star*m_dot_total/Pc #throat area [m^2]
d_t = 2 * np.sqrt(A_t/np.pi) #throat diameter [m]
A_e = A_t * ER #exit area [m^2]
d_e = 2 * np.sqrt(A_e/np.pi) #throat diameter [m]

#CHAMBER CALCULATIONS
A_c = CR*A_t
d_c = 2*np.sqrt(A_c/np.pi) #chamber diameter [m]
V_c = L_star * A_t #chamber colume [m^3]
L_c = V_c/A_c #chamber length (from injector face to throat) [m]

#PRINT OUTPUTS
print(f"Total mass flow: {m_dot_total: .3f} kg/s")
print(f"Fuel mass flow: {m_dot_fuel: .3f} kg/s")
print(f"Oxidizer mass flow: {m_dot_ox: .3f} kg/s")

print(f"Chamber diameter: {d_c: .6f} m")
print(f"Throat diameter: {d_t: .6f} m")
print(f"Exit diameter: {d_e: .6f} m")
print(f"Chamber lengthL {L_c: .6f} m")



