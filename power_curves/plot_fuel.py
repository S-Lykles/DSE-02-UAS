import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
import seaborn as sns
import const
# import dual_phase.inputs as inputs
from power_curves.mass_frac import fuel_weight
# from aero.cl_cd import dragpolar_dual


def plot_fuel(b, S, polar, b_design=6, S_design=3.763, which='endurance'):

    fuel_w = np.zeros((len(b), len(S)))
    A = np.zeros((len(b), len(S)))
    unit_conversion = 0.45359237 / 745.7/ 60**2
    SFC = 0.5 * unit_conversion
    for j in range(len(S)):
        for i in range(len(b)):
            CL, CD = polar(b[i], S[j])
            fw = fuel_weight(CL, CD, SFC, S[j], which=which, v_cruise=const.v_cruise*1.1)
            fuel_w[j][i] = fw
            A[j][i] = b[i]**2/S[j]

    # fuel_w = pd.DataFrame(fuel_w)
    # # Create a heatmap
    original_cmap = plt.get_cmap('viridis')
    new_cmap = mcolors.LinearSegmentedColormap.from_list(
        'truncated_viridis', original_cmap(np.linspace(0.15, 1, 256))
    )
    plt.figure()
    CS = plt.contour(b, S, fuel_w, levels=8, colors='k')
    CSA = plt.contour(b, S, A, levels=8, colors='white', linestyles='dashed', linewidths=0.5)
    plt.clabel(CS, inline=1, fontsize=10, colors='k')
    plt.clabel(CSA, inline=1, fontsize=10, colors='white')
    plt.contourf(b, S, fuel_w, levels=100, cmap=new_cmap)
    plt.scatter(b_design, S_design, c='r', marker='x', label='Design point')
    # sns.heatmap(fuel_w, annot=False, xticklabels=b, yticklabels=S, cmap='viridis')
    plt.xlabel('b [m]')
    plt.ylabel('S [$m^2$]')
    plt.legend()
    # plt.title('Fuel weight')
    plt.show()


def plot_cd0(b, S, CD0func, b_design=6, S_design=3.763):
    cd0 = np.zeros((len(b), len(S)))
    unit_conversion = 0.45359237 / 745.7/ 60**2
    SFC = 0.5 * unit_conversion
    for i in range(len(b)):
        for j in range(len(S)):
            cd0[i][j] = CD0func(b[i], S[j])


    # fuel_w = pd.DataFrame(fuel_w)
    # # Create a heatmap
    original_cmap = plt.get_cmap('viridis')
    new_cmap = mcolors.LinearSegmentedColormap.from_list(
        'truncated_viridis', original_cmap(np.linspace(0.15, 1, 256))
    )
    plt.figure()
    CS = plt.contour(b, S, cd0, levels=8, colors='k')
    plt.clabel(CS, inline=1, fontsize=10, colors='k')
    plt.contourf(b, S, cd0, levels=100, cmap=new_cmap)
    # plt.scatter(inputs.b, inputs.S, c='r', marker='x')
    # sns.heatmap(fuel_w, annot=False, xticklabels=b, yticklabels=S, cmap='viridis')
    plt.xlabel('b [m]')
    plt.ylabel('S [$m^2$]')
    # plt.title('Fuel weight')
    plt.show()

# W = const.MTOW
# N_r = 4
# DL = np.arange(100,400+n,n)
# R = np.sqrt((W/N_r)/(DL*np.pi))
# plt.plot(DL,R)
# plt.xlabel("Disk loading [N/m^2]")
# plt.ylabel("Radius [m]")
# plt.show()


# #Example: JUMP20
# MTOW_e = 97.5*9.80665
# D = np.sqrt((87-193)**2+(204-249)**2)/np.sqrt((52-648)**2+(388-53)**2)*5.7
# bar = np.sqrt((133-342)**2+(241-332)**2)/np.sqrt((52-648)**2+(388-53)**2)*5.7
# R = D/2
# N_r = 4
# DL = (MTOW_e/N_r)/(np.pi*R**2)
# print(DL,R,bar)



# #Example: FD180P long endurance heavy VTOL fixed-wing UAV
# MTOW_e = 180*9.80665
# D = np.sqrt((349-466)**2+(538-452)**2)/np.sqrt((6-795)**2+(453-376)**2)*6.5
# R = D/2
# N_r = 4
# DL = (MTOW_e/N_r)/(np.pi*R**2)
# print(DL,R)

