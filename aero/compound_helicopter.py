import numpy as np
from matplotlib import pyplot as plt
from math import pi, sqrt
from cl_cd import *
from titl_wing import *
import const

# Set up cl and cd for plot
def dragpolar_comp(b,S,d_eng=0.5,N_eng=2,Lambda=0.45,h=500,v=const.v_cruise,c=None,Sf=2,CL_start=0.,CL_end=1.2,CL_step=1000):
    if c == None:
        c = S/b
    d_eng = d_eng/1.5                    # Compared to engine diameter of tilt wing
    cd_parasite = 0.04/0.67*1.2
    cl,cd_tilt_wing, cd_prop = dragpolar_tilt_wing(b,S,h,v,c,Sf,d_eng,N_eng,Lambda,CL_start,CL_end,CL_step)
    cd0_eng = N_eng * 0.006
    cl = np.linspace(CL_start,CL_end,CL_step)
    A = b**2/S
    e = 1.78*(1-0.045*A**0.68)-0.64 
    cd_comp = cd_parasite + cd_prop + cd0_eng + cl**2/(pi * A * e)
    # print("cd_profile = ",cd_profile)
    # print("cd_parasite2 =",cd_parasite)
    # print("cd_0 =",cd_parasite+cd_profile)
    return cl, cd_comp


# Other method, do not use

# def profile_drag(h,v,c,Sf):
#     # Obtain Xfoil data
#     data_23012 = np.loadtxt(r"C:\Users\Florian Tjepkema\Documents\Aerospace Engineering\AeroSpace 2023-2024\DSE\DSE-02-UAS\aero\naca23012.txt")
#     cd_pressure = data_23012[4, 3]                              # Pressure drag coefficient at alpha = 2
#     cd0_fus,cf,Re, rho, T, p, M = cd0_fuselage(h,v,c,Sf)
#     cd_profile = cd_pressure + cf
#     return cd_profile

# # Calculate parasite drag coefficient with equations from literature
# def parasite_drag(MTOM,S):
#     d_q_airframe_wing = 1.6*(2.2046*MTOM/1000)**(2/3)                            # ft^2 d/q for fuselage and wing
#     d_q_rotors = 0.4*(2.2046*MTOM/1000)**(2/3)                                   # ft^2 d/q for rotors
#     cd_interference = 0.003                                                      # [-] cd 
#     cd_parasite = 0.0929* (d_q_airframe_wing + d_q_rotors)/S + cd_interference   #conversion to m2, then make it unitless by /S        [-]
#     return cd_parasite


