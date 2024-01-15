from DSE.aero.xfoil_plots import *
from math import pi, atan

# def aero_values():
c_bar = 0.598
Cm_0_airfoil= cm_23012_airfoil[10]
Cm_alpha = (cm_23012_airfoil[20]-cm_23012_airfoil[0])/((alpha_23012_airfoil[20]-alpha_23012_airfoil[0])*pi/180)
alpha_0 = 0
sweep_ang_25_c_rad = 0*pi/180
sweep_ang_50_c_rad = -2.386*pi/180
CL_alpha_wing = (cl_23012_wing[20]-cl_23012_wing[0])/((alpha_23012_wing[20]-alpha_23012_wing[0])*pi/180)
CD0_wing = 0.0077
# Cl_alpha_h
CL_max = max(cl_23012_airfoil)
Cl_max_h = max(cl_0012_airfoil)
CL_0 = cl_23012_wing[10]
c_root = 0.833
c_tip = 0.333
CL_cruise = cl_23012_wing[9]
Cl_cruise_h = cl_0012_airfoil[33]
b = 6
S = 3.5
S_h = 1.26
Cl_alpha_v = (cl_0012_airfoil[20]-cl_0012_airfoil[0])/((alpha_0012_airfoil[20]-alpha_0012_airfoil[0])*pi/180)
Cl_alpha_h = Cl_alpha_v
sweep_ang_rad = atan(0.125/3)

# Cm_0_airfoil =
# sweep_ang_14_c_rad =
# CL_alpha_w =
# sweep_ang_12_Ch_rad =
print(Cl_max_h)
print(Cl_cruise_h)
print("hoi")


    # return c_bar, cm_0, Cm_alpha, alpha_0, sweep_ang_25_c_rad, sweep_ang_50_c_rad, Cl_alpha_wing, CL_max, CL_0, c_root, c_tip, CL_cruise, b, S
