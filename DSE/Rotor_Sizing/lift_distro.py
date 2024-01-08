import pandas as pd
from DSE.Rotor_Sizing import *
from DSE import const
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

file_dir = Path(__file__).parent

def data_file():
    data = pd.read_csv(file_dir / r"clarky.csv")

    x_axis = 'Alpha'
    y_axis = 'Cl'
    z_axis = 'Cd'

    return data[x_axis], data[y_axis], data[z_axis]

## Actual calculations
def lift_coefficient(r):
    a, cl, cd = data_file()
    a_blade = 14 - 22 * r
    return np.interp(a_blade, a, cl)

def drag_coefficient(r):
    a, cl, cd = data_file()
    a_blade = 14 - 22 * r
    return np.interp(a_blade, a, cd)

def chord_length(r):
    return 0.03 - (0.06 * r)

def lift_blade(r, rpm):
    Cl = lift_coefficient(r)
    rho = 1.225
    omega = (rpm * 2 * np.pi) / 60
    c = chord_length(r)
    return 0.5 * Cl * rho * (omega * r)**2 * c

def drag_blade(r, rpm):
    Cd = 0.01092
    rho = 1.225
    omega = (rpm * 2 * np.pi) / 60
    c = chord_length(r)
    return 0.5 * Cd * rho * (omega * r) ** 2 * c

#Rotor span
lower_limit = 0.0
upper_limit = 0.5
rpm_values = np.linspace(2000, 8000, num=20)
N = 8

r_values = np.linspace(lower_limit, upper_limit, 100)

#Lift distribution
lift_values = np.zeros((len(rpm_values), len(r_values)))

total_lift_values = np.zeros(len(rpm_values))

for i, rpm in enumerate(rpm_values):
    lift_values[i, :] = [lift_blade(r, rpm) for r in r_values]
    total_lift_values[i] = np.trapz(lift_values[i, :], r_values)
    #print(f"Total lift for RPM = {rpm}:", total_lift_values[i])

    #plt.plot(r_values, lift_values[i, :], label=f'RPM = {rpm}')
    #plt.fill_between(r_values, 0, lift_values[i, :], alpha=0.2)

# plt.title('Lift distribution over blade')
# plt.xlabel('Radius [r]')
# plt.ylabel('Lift [N]')
# plt.grid()
# plt.legend()
# plt.show()

plt.figure()
plt.plot(rpm_values, total_lift_values)
plt.axhline(y=441.45/N, color='grey', linestyle='--')
plt.xlabel('RPM')
plt.ylabel('Total Lift')
plt.title('Total Lift for Different RPMs')
plt.grid(True)
plt.show()

#Drag
drag_values = np.zeros((len(rpm_values), len(r_values)))

total_drag_values = np.zeros(len(rpm_values))

for i, rpm in enumerate(rpm_values):
    drag_values[i, :] = [drag_blade(r, rpm) for r in r_values]
    total_drag_values[i] = np.trapz(drag_values[i, :], r_values)
    #print(f"Total lift for RPM = {rpm}:", total_drag_values[i])

plt.figure()
plt.plot(rpm_values, total_drag_values)
plt.xlabel('RPM')
plt.ylabel('Total Drag')
plt.title('Total Drag for Different RPMs')
plt.grid(True)
plt.show()
