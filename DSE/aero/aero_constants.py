from DSE.aero.xfoil_plots import *
from math import pi, atan
import numpy as np

# Obtain list for Cl/Cd and find max for cruise condition
cl_list = np.array([cl_23012_airfoil])
cd_list = np.array([cd_23012_airfoil])
c = cl_list/cd_list

# Get index from cruise condition
index_cruise_airfoil = np.argmax(c)
alpha_cruise_airfoil = alpha_23012_airfoil[index_cruise_airfoil]
index_cruise_wing = np.where(alpha_23012_wing == alpha_cruise_airfoil)
index_cruise_airfoil_h = np.where(alpha_0012_airfoil == alpha_cruise_airfoil)


alpha_0 = 0
b = 6
#c_bar = 0.598

# Cl relates to airfoil only, CL to the entire wing. Same for Cd and CD about drag

CD0_wing = 0.0077
CD_alpha_wing = (cd_23012_wing[20]-cd_23012_wing[0])/((alpha_23012_wing[20]-alpha_23012_wing[0])*pi/180)
CL_0 = cl_23012_wing[np.where(alpha_23012_wing == 0)]
CL_alpha_wing = (cl_23012_wing[20]-cl_23012_wing[0])/((alpha_23012_wing[20]-alpha_23012_wing[0])*pi/180)
Cl_alpha_v = (cl_0012_airfoil[20]-cl_0012_airfoil[0])/((alpha_0012_airfoil[20]-alpha_0012_airfoil[0])*pi/180)
Cl_alpha_h = Cl_alpha_v
CL_cruise = cl_23012_wing[index_cruise_wing]
CD_cruise = cd_23012_wing[index_cruise_wing]
print(CD_cruise)
Cl_cruise = cl_23012_airfoil[index_cruise_airfoil]
Cl_cruise_h = cl_0012_airfoil[index_cruise_airfoil_h]
CL_max = max(cl_23012_wing)
Cl_max_h = max(cl_0012_airfoil)

Cm_0_airfoil= cm_23012_airfoil[np.where(alpha_23012_airfoil == 0)]
Cm_alpha = (cm_23012_airfoil[20]-cm_23012_airfoil[0])/((alpha_23012_airfoil[20]-alpha_23012_airfoil[0])*pi/180)

c_root = 0.833
c_tip = 0.333
tc = 0.12
taper = 0.4
c_bar = c_root * (2/3) * ((1+taper+taper**2)/(1+taper))


sweep_ang_rad = atan(0.125/3)
sweep_ang_25_c_rad = 0*pi/180
sweep_ang_50_c_rad = -2.386*pi/180
S = 3.5
S_h = 1.26

AR = b**2 / S
e = 1.78*(1 - 0.045 * AR **0.68) - 0.64

# winglet data
tc_winglet = tc
taper_winglet = taper
sweep_ang_25_c_rad_winglet = 45 * const.deg2rad
S_winglet = 0.21
b_winglet = 0.45
c_root_winglet = 0.333

Cl_alpha_wing = (cl_23012_airfoil[20]-cl_23012_airfoil[0])/((alpha_23012_airfoil[20]-alpha_23012_airfoil[0])*pi/180)
print(Cl_alpha_h)