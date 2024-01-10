from DSE.aero.xfoil_plots import *
from math import pi

# def aero_values():
c_bar = 0.598
cm_0 = cm_23012[10]
Cm_alpha = (cm_23012[20]-cm_23012[0])/(alpha_23012[20]-alpha_23012[0])*pi/180
alpha_0 = 0
sweep_ang_25_c_rad = 0*pi/180
sweep_ang_50_c_rad = 0*pi/180
Cl_alpha_wing = (cl_23012[20]-cl_23012[0])/(alpha_23012[20]-alpha_23012[0])*pi/180
# Cl_alpha_h
CL_max = max(cl_23012)
CL_0 = cl_23012[10]
c_root = 0.833
c_tip = 0.333
CL_cruise = cl_23012[9]
b = 6
S = 3.5
    # return c_bar, cm_0, Cm_alpha, alpha_0, sweep_ang_25_c_rad, sweep_ang_50_c_rad, Cl_alpha_wing, CL_max, CL_0, c_root, c_tip, CL_cruise, b, S

# print(aero_values())