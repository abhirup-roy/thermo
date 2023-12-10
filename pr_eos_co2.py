import numpy as np
import matplotlib.pyplot as plt

# Constants for CO2
T_C = 304.18  # K
P_C = 73.8e5  # Pa
w = 0.2249

R = 8.314  # kJ/kmol K

# Temperature and Volume ranges
T = np.array([200, 300, 500, 1000])  # K
v_m = np.linspace(0.001, 1, 1000)


# Calulate a and b in eqn.
a = 0.45724 * R ** 2 * T_C ** 2 / P_C
b = 0.0778 * R * T_C / P_C

for temp in T:
    # Calculate alpha and pressure for given temp
    T_R = temp / T_C
    alpha = (1 + (0.37464 + 1.54226 * w - 0.26992 * w ** 2) * (1 - T_R ** 0.5)) ** 2

    P = R * temp / (v_m - b) - a * alpha / (v_m ** 2 + 2 * b * v_m - b ** 2)
    # Plot P-V diagram
    plt.plot(v_m, P, label=f'T = {temp} K')

# Plot critical temp
plt.axhline(y=P_C, linestyle = 'dashed', color="k", label = "Critical Pressure")

plt.xlabel('Specific Volume ($m^3$/kmol)')
plt.ylabel('Pressure (Pa)')
plt.legend()
plt.title('P-V Diagram - $CO_2$ (Peng-Robinson Equation of State)')
plt.grid(True)
plt.show()