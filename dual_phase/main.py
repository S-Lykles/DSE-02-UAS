import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import const
import power_curves.inputs as i
from power_curves.mass_frac import fuel_weight
from aero.cl_cd import dragpolar_dual

n = 0.5
b = np.arange(2,6+n,n)
S = np.arange(1,5+n,n)
print(len(b), len(S))

fuel_w = np.zeros((len(b), len(S)))
unit_conversion = 0.45359237 / 745.7/ 60**2
SFC = 0.5 * unit_conversion
for i in range(len(b)):
    for j in range(len(S)):
        CL, CD = dragpolar_dual(b[i], S[j], CL_start=0.01)
        fw = fuel_weight(CL, CD, SFC, S[j], which='endurance', v_cruise=const.v_cruise*1.1)
        fuel_w[i][j] = fw


fuel_w = pd.DataFrame(fuel_w)
# # Create a heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(fuel_w, annot=True, xticklabels=b, yticklabels=S, cmap='viridis')
plt.xlabel('b')
plt.ylabel('S')
plt.title('Fuel weight')
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

