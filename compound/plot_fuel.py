import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import const
import power_curves.inputs as i
from power_curves.mass_frac import fuel_weight
from aero.compound_helicopter import dragpolar_comp

n = 0.5
b = np.arange(2,6+n,n)
S = np.arange(1,5+n,n)
print(len(b), len(S))

fuel_w = np.zeros((len(b), len(S)))
unit_conversion = 0.45359237 / 745.7/ 60**2
SFC = 0.5 * unit_conversion
for i in range(len(b)):
    for j in range(len(S)):
        CL, CD = dragpolar_comp(b[i], S[j], CL_start=0.01)
        fw = fuel_weight(CL, CD, SFC, S[j], which='payload', v_cruise=const.v_cruise*1.1)
        fuel_w[i][j] = fw


fuel_w = pd.DataFrame(fuel_w)
# # Create a heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(fuel_w, annot=True, xticklabels=b, yticklabels=S, cmap='viridis')
plt.xlabel('b')
plt.ylabel('S')
plt.title('Fuel weight')
plt.show()