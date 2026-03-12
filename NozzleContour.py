import numpy as np
import matplotlib.pyplot as plt
from BasicSizing import BasicSizing


# REFERENCES 
# http://www.aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf
# https://rrs.org/2023/01/28/making-correct-parabolic-nozzles/


# INITIALIZATION
mode = "Hotfire"
sizing = BasicSizing(mode)


# INPUT PARAMETERS
r_c = sizing.d_c / 2      # Chamber radius [m]
r_t = sizing.d_t / 2      # Throat radius [m]
r_e = sizing.d_e / 2      # Exit radius [m]
L_c = sizing.L_c          # Chamber length [m]
L_n = sizing.L_n          # Nozzle length [m]
ER = sizing.ER

theta_n = np.deg2rad(sizing.theta_n)
theta_e = np.deg2rad(sizing.theta_e)
convergence_angle = 37.5  # [deg]
beta = np.deg2rad(convergence_angle)


# THROAT GEOMETRY (RAO METHOD)
# Converging Throat Arc (Radius 1.5 * r_t)
# Starts at (-90 - beta) and ends at -90 degrees
theta_conv_range = np.linspace(-np.pi/2 - beta, -np.pi/2, 50)
x_converging = 1.5 * r_t * np.cos(theta_conv_range)
y_converging = 1.5 * r_t * np.sin(theta_conv_range) + 1.5 * r_t + r_t 

# Diverging Throat Arc (Radius 0.382 * r_t)
# Starts at -90 and ends at (theta_n - 90)
theta_div_range = np.linspace(-np.pi/2, theta_n - (np.pi/2), 200)
x_divergence = 0.382 * r_t * np.cos(theta_div_range)
y_divergence = 0.382 * r_t * np.sin(theta_div_range) + 0.382 * r_t + r_t 


# BELL GEOMETRY (PARABOLIC NOZZLE)
# Point N: Start of the parabola
x_N = x_divergence[-1]
y_N = y_divergence[-1]

# Point E: Exit of the nozzle
x_E = L_n
y_E = r_e

# Point Q: Intersection of tangents from N and E
m1 = np.tan(theta_n)
m2 = np.tan(theta_e)
c1 = y_N - m1 * x_N 
c2 = y_E - m2 * x_E

x_Q = (c2 - c1) / (m1 - m2) 
y_Q = (m1 * c2 - m2 * c1) / (m1 - m2)

# Generate Bell Curve
t = np.linspace(0, 1, 30)
x_bell = ((1 - t)**2) * x_N + 2 * (1 - t) * t * x_Q + (t**2) * x_E
y_bell = ((1 - t)**2) * y_N + 2 * (1 - t) * t * y_Q + (t**2) * y_E


# CHAMBER & CONE TRANSITION
# P1: Interface between throat entry and convergence cone
x_p1 = -1.5 * r_t * np.sin(beta) 
y_p1 = 1.5 * r_t * (1 - np.cos(beta)) + r_t

# P2: Start of the shoulder (where cylinder ends)
y_sh_end = r_c - 1.5 * r_t * (1 - np.cos(beta))
x_gap = (y_sh_end - y_p1) / np.tan(beta)
x_sh_width = 1.5 * r_t * np.sin(beta)

x_p2 = x_p1 - x_gap - x_sh_width
y_p2 = r_c

# Shoulder Curve (Inward arc from cylinder to cone)
t_sh = np.linspace(0, beta, 30)
x_sh = x_p2 + 1.5 * r_t * np.sin(t_sh)
y_sh = r_c - 1.5 * r_t * (1 - np.cos(t_sh))

# Straight Transition (The Cone)
x_straight = np.linspace(x_sh[-1], x_p1, 20)
y_straight = np.linspace(y_sh[-1], y_p1, 20)

x_trans = np.concatenate([x_sh, x_straight])
y_trans = np.concatenate([y_sh, y_straight])


# FINAL ASSEMBLY
L_convergent = abs(x_p2)
L_cylindrical = max(0, L_c - L_convergent)

if L_c < L_convergent:
    print(f"Warning: L_c ({L_c}) is shorter than convergent section ({L_convergent:.4f})!")

# Cylindrical Chamber Section
x_chamber = np.array([-L_c, x_p2])
y_chamber = np.array([r_c, r_c])

# Combine all coordinate arrays
x_full = np.concatenate([x_chamber, x_trans, x_converging, x_divergence, x_bell])
y_full = np.concatenate([y_chamber, y_trans, y_converging, y_divergence, y_bell])

# Sort and remove duplicates for plotting
_, unique_indices = np.unique(x_full, return_index=True)
x_plot = x_full[np.sort(unique_indices)]
y_plot = y_full[np.sort(unique_indices)]


# PLOTTING
plt.figure(figsize=(12, 5))
plt.plot(x_plot, y_plot, color="black", linewidth=2, label="Nozzle Wall")
plt.plot(x_plot, -y_plot, color="black", linewidth=2)

# Reference Lines
plt.axvline(x=0, color="red", linestyle="--", alpha=0.5)
plt.text(0.005, r_t * 1.5, 'Throat ($x=0$)', color='red', fontweight='bold')

plt.fill_between(x_plot, -y_plot, y_plot, color='skyblue', alpha=0.1)
plt.title(f"Rao Parabolic Nozzle Profile (Expansion Ratio: {ER:.2f})")
plt.xlabel("Length [m]")
plt.ylabel("Radius [m]")
plt.axis("equal")
plt.grid(True, linestyle=":", alpha=0.7)
plt.tight_layout()
plt.show()