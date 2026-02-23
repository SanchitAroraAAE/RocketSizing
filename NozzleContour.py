import numpy as np
import matplotlib.pyplot as plt

# http://www.aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf
# https://rrs.org/2023/01/28/making-correct-parabolic-nozzles/
# https://wikis.mit.edu/confluence/pages/viewpage.action?pageId=153816550

#INPUT PARAMETERS
ER = 3.4826 #expansion ratio
CR = 5 #contraction ratio
r_c = 0.057456/2 #chamber rad.
r_t = 0.028728/2 #throat rad.
r_e = np.sqrt(ER)*r_t #exit rad.
L_c = 0.0762 #chamber length [m]

percent_bell = 0.8 #Percent Rao nozzle. 80% standard
theta_n = np.deg2rad(20.7239)
theta_e = np.deg2rad(14.5174)
convergence_angle = 37.5


#CALCULATIONS
#Area calcs
A_c = np.pi * (r_c**2) #chamber area
A_t = np.pi * (r_t**2) #throat area
A_e = np.pi * (r_e**2) #exit area

#Length calcs
L_n = percent_bell * (np.sqrt(ER)-1)*r_t/(np.tan(np.deg2rad(15))) #nozzle length eqt. 3

#THROAT AND BELL GEOMETRY
beta = np.deg2rad(convergence_angle)

#Rao converging section
# Starts exactly at -90 degrees minus beta (to meet the cone) and ends at -90
theta_convergence = np.linspace(-np.pi/2 - beta, -np.pi/2, 50)
x_converging = 1.5 * r_t * np.cos(theta_convergence)
y_converging = 1.5 * r_t * np.sin(theta_convergence) + 1.5 * r_t + r_t #eqt. 4

#Rao diverging section
theta_divergence = np.linspace(-np.pi/2, theta_n-(np.pi/2), 200)
x_divergence = 0.382 * r_t * np.cos(theta_divergence)
y_divergence = 0.382 * r_t * np.sin(theta_divergence) + 0.382 * r_t + r_t #eqt. 5

#Finding N
x_N = x_divergence[-1]
y_N = y_divergence[-1] #found by setting angle to (theta_n -90) in eqt. 5

#Finding E
x_E = L_n #derived from eqt. 3
y_E = r_e #derived from eqt. 2 or just BasicSizing results

#Finding Q
m1 = np.tan(theta_n) #eqt. 8
m2 = np.tan(theta_e) #eqt. 8
c1 = y_N - m1*x_N #eqt. 9
c2 = y_E - m2*x_E #eqt. 9
x_Q =(c2 - c1)/(m1 - m2) #eqt. 10
y_Q = (m1*c2 - m2*c1)/(m1 - m2) #eqt. 10

#Rao bell section
t = np.linspace(0,1,30)
x_bell = ((1-t)**2) * x_N + 2*(1-t)*t*x_Q + (t**2)*x_E
y_bell = ((1-t)**2) * y_N + 2*(1-t)*t*y_Q + (t**2)*y_E

#CHAMBER TRANSITION GEOMETRY
# Define P1 (start of throat entry)
x_p1 = -1.5 * r_t * np.sin(beta) # The x-coord where the throat arc meets the cone
y_p1 = 1.5 * r_t * (1 - np.cos(beta)) + r_t

# Define p2 (End of straight chamber / Start of shoulder)
# x_p2 is calculated by finding the horizontal distance of the shoulder + the straight cone
y_sh_end = r_c - 1.5 * r_t * (1 - np.cos(beta))
x_gap = (y_sh_end - y_p1) / np.tan(beta)
x_sh_width = 1.5 * r_t * np.sin(beta)

x_p2 = x_p1 - x_gap - x_sh_width
y_p2 = r_c

# Create the Shoulder Curve (Inwards curve)
t_sh = np.linspace(0, beta, 30)
x_sh = x_p2 + 1.5 * r_t * np.sin(t_sh)
y_sh = r_c - 1.5 * r_t * (1 - np.cos(t_sh))

# Create the Straight Transition (The Cone)
x_straight = np.linspace(x_sh[-1], x_p1, 20)
y_straight = np.linspace(y_sh[-1], y_p1, 20)

# Combine the shoulder curve and the straight cone into one transition segment
x_trans = np.concatenate([x_sh, x_straight])
y_trans = np.concatenate([y_sh, y_straight])


#FINAL ARRAY COMBINING
# Chamber cylinder
x_chamber = np.array([x_p2 - L_c, x_p2])
y_chamber = np.array([r_c, r_c])

# Concatenate all parts
x_full = np.concatenate([x_chamber, x_trans, x_converging, x_divergence, x_bell])
y_full = np.concatenate([y_chamber, y_trans, y_converging, y_divergence, y_bell])

# Clean up duplicates and sort
_, unique_indices = np.unique(x_full, return_index=True)
x_plot = x_full[np.sort(unique_indices)]
y_plot = y_full[np.sort(unique_indices)]

#PLOTTING
plt.figure(figsize=(10, 4))
plt.plot(x_plot, y_plot, color="black", linewidth=2)
plt.plot(x_plot, -y_plot, color="black", linewidth=2)
plt.fill_between(x_plot, -y_plot, y_plot, color='lightgray', alpha=0.3)
plt.title(f"Rao Nozzle Profile (ER={ER:.2f})")
plt.xlabel("Length (m)")
plt.ylabel("Radius (m)")
plt.axis("equal")
plt.grid(True, linestyle='--', alpha=0.6)
plt.show()

