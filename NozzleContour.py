import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

#http://www.aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf

#INPUT PARAMETERS
ER = 3.4826 #expansion ratio
CR = 4 #contraction ratio
r_c = 0.057456/2 #chamber rad.
r_t = 0.028728/2 #throat rad.
r_e = np.sqrt(ER)*r_t #exit rad.
L_star = 1.0 #characteristic length
L_c = 0.25 #chamber length

percent_bell = 0.8 #Percent Rao nozzle. 80% standard
theta_n = 21 * (np.pi/180) #copied from PSP-Liq CMS due to similar exp. ratios (I'm lazy)
theta_e = 15 * (np.pi/180) #copied from PSP-Liq CMS due to similar exp. ratios
convergence_angle = 37.5


#CALCULATIONS
#Rao converging section
theta_convergence = np.linspace(-3/4*np.pi, -np.pi/2)
x_convergence = 1.5 * r_t * np.cos(theta_convergence)
y_convergence = 1.5 * r_t * np.sin(theta_convergence) + 1.5 * r_t + r_t

A_c = np.pi * (r_c**2)
A_t = np.pi * (r_t**2)
A_e = np.pi * (r_e**2)

#Rao diverging section
theta_divergence = np.linspace(-np.pi/2, theta_n-(np.pi/2), 200)
x_divergence = 0.382 * r_t * np.cos(theta_divergence)
y_divergence = 0.382 * r_t * np.sin(theta_divergence) + 0.382 * r_t + r_t

#Finding E
x_E = percent_bell * (((np.sqrt(ER)-1)*r_t)/np.tan(np.pi/12))
y_E = r_e

#Finding N
x_N = x_divergence[-1]
y_N = y_divergence[-1]

#Finding Q
m1 = np.tan(theta_n)
m2 = np.tan(theta_e)
c1 = y_N - m1*x_N
c2 = y_E - m2*x_E
x_Q =(c2 - c1)/(m1 - m2)
y_Q = (m1*c2 - m2*c1)/(m1 - m2)

#Rao bell section
t = np.linspace(0,1,30)
x_bell = ((1-t)**2) * x_N + 2*(1-t)*t*x_Q + t**2*x_E
y_bell = ((1-t))**2 * y_N + 2*(1-t)*t*y_Q + t**2*y_E

x_p1 = min(x_convergence)
y_p1 = max(y_convergence)

Len_ConvergeCombustor = (r_c - max(y_convergence)) / np.tan(np.deg2rad(convergence_angle))
x_p2 = x_p1 - Len_ConvergeCombustor
y_p2 = r_c; #Determined from r_t input and ER

#Solve for intersections at point 3
f = sp.symbols('f')

line_p1 = sp.tan(-1 * (sp.pi / 4)) * (f - x_p1) + y_p1
line_p2 = sp.tan(-1*np.deg2rad(convergence_angle)) * (f - x_p2) + y_p2

eqn = sp.Eq(line_p1, line_p2)
s = sp.solve(eqn, f)[0]

x_p3 = float(s)
y_p3 = float(sp.tan(-1*np.deg2rad(convergence_angle)) * (x_p3 - x_p2) + y_p2)


#Converging chamber bezier
x_ConvergingBezier = ((1 - t)**2 * x_p2 + 2*(1 - t)*t * x_p3 + t**2 * x_p1)
y_ConvergingBezier = ((1 - t)**2 * y_p2 + 2*(1 - t)*t * y_p3 + t**2 * y_p1)

#Combine partial coords
x_Partial = np.concatenate([x_ConvergingBezier, x_convergence, x_divergence, x_bell])
y_Partial = np.concatenate([y_ConvergingBezier, y_convergence, y_divergence, y_bell])

# Keep only unique x values
unique_vals, unique_index = np.unique(x_Partial, return_index=True)
unique_index_sorted = np.sort(unique_index)

x_Partial = x_Partial[unique_index_sorted]
y_Partial = y_Partial[unique_index_sorted]

x_fullConverge = np.concatenate([x_ConvergingBezier, x_convergence])
y_fullConverge = np.concatenate([y_ConvergingBezier, y_convergence])

# Chamber volume calculations
#V_chamber = L_star * A_t
#V_converge = abs(np.pi * np.trapezoid(np.flip(x_fullConverge), np.flip(y_fullConverge**2)))
#V_cylinder = V_chamber - V_converge
#L_c = V_cylinder / (np.pi * r_c**2)

# Final coordinate matrix
x_Chamber = np.array([
    np.min(x_ConvergingBezier) - L_c,
    np.min(x_ConvergingBezier)
])

y_Chamber = np.array([r_c, r_c])

x_Total = np.concatenate([x_Chamber, x_ConvergingBezier, x_convergence, x_divergence, x_bell])
y_Total = np.concatenate([y_Chamber, y_ConvergingBezier, y_convergence, y_divergence, y_bell])
z_Total = np.zeros(len(x_Total))

# Keep only unique x values
unique_vals, unique_index = np.unique(x_Total, return_index=True)
unique_index_sorted = np.sort(unique_index)

Total = np.column_stack([x_Total[unique_index_sorted], y_Total[unique_index_sorted], z_Total[unique_index_sorted]])

x_plot = x_Total[unique_index_sorted]
y_plot = y_Total[unique_index_sorted]

# Write to file
np.savetxt("nozzle_contour.txt", Total, delimiter="\t", fmt="%.5f")

# PLOTTING
plt.figure()
plt.plot(x_plot, y_plot, color="black")
plt.plot(x_plot, -1*y_plot, color="black")
plt.title("Nozzle Contour")
plt.xlabel("meters")
plt.ylabel("meters")
plt.show()


