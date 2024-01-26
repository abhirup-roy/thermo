import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Constants for CO2
T_C = 304.18  # K
P_C = 7.38e6  # Pa
w = 0.2249

R = 8.314  # kJ/kmol K

# # Temperature and Volume ranges
T = [T_C-60,T_C-40,T_C-20,T_C, T_C+20, T_C+40, T_C+60]  # K
V = np.linspace(0.45e-4, 3e-4, 10000)


# Calulate a and b in eqn.
a = 0.45724 * R ** 2 * T_C ** 2 / P_C
b = 0.0778 * R * T_C / P_C


for temp in T:
    # Calculate alpha and pressure for given temp
    T_R = temp / T_C
    alpha = (1 + (0.37464 + 1.54226 * w - 0.26992 * w ** 2) * (1 - T_R ** 0.5)) ** 2
    P = (R * temp / (V - b)) - ((a * alpha) / (V ** 2 + 2 * b * V - b ** 2))
    # Plot P-V diagram
    if T == T_C:
        plt.plot(V*1e6, P/1e6, label=f'T = {temp} K ($T_C$)')
    else:
        plt.plot(V*1e6, P/1e6, label=f'T = {temp} K')


df_sat = pd.read_excel(r"C:\Users\abhir\Documents\ChemE Y3\ART\Thermo\Coursework\co2_sat_data.xlsx")
m_r = 44.0095
df_sat["vf (cm3/mol)"] = df_sat["vf (m3/kg)"] * m_r * 1e3
df_sat["vg (cm3/mol)"] = df_sat["vg (m3/kg)"] * m_r * 1e3
df_sat["T (K)"] =  df_sat["T (degC)"] + 273
df_sat.drop(["vf (m3/kg)", "vg (m3/kg)", "T (degC)"], axis=1, inplace=True)
# print(df_sat.head())
# print(V)

v_f = df_sat.drop(["vg (cm3/mol)", "T (K)"], axis=1)
v_f = v_f[v_f["vf (cm3/mol)"]>=45]

v_g = df_sat.drop(["vf (cm3/mol)", "T (K)"], axis=1)
v_g = v_g[v_g["vg (cm3/mol)"]<=300]

plt.plot(v_f["vf (cm3/mol)"], v_f["P (MPa)"], linestyle = '--', color="black", label="Saturation Plot")
plt.plot(v_g["vg (cm3/mol)"], v_g["P (MPa)"], linestyle = '--', color="black")

Vc = 91.9 #cm3/mol
plt.plot(Vc,P_C*1e-6, marker="x", color="black", label="Critical Point")


# Plot critical temp
# plt.axhline(y=P_C, linestyle = 'dashed', color="k", label = "Critical Pressure")
ax = plt.gca()
ax.set_ylim(top=35, bottom=0)

plt.xlabel('Specific Volume ($cm^3$/mol)')
plt.ylabel('Pressure (MPa)')
plt.legend(loc="upper right")
plt.title('P-V Diagram - $CO_2$ (Peng-Robinson Equation of State)')
plt.grid(True)
plt.show()  