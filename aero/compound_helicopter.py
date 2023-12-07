import numpy as np
from matplotlib import pyplot as plt
from math import pi, sqrt
from cl_cd import *
import const




def profile_drag(h,v,c,Sf):
    # Obtain Xfoil data
    data_23012 = np.loadtxt(r"C:\Users\Florian Tjepkema\Documents\Aerospace Engineering\AeroSpace 2023-2024\DSE\DSE-02-UAS\aero\naca23012.txt")
    cd_pressure = data_23012[4, 3]                              # Pressure drag coefficient at alpha = 2
    cd0_fus,cf,Re, rho, T, p, M = cd0_fuselage(h,v,c,Sf)
    cd_profile = cd_pressure + cf
    return cd_profile

# Calculate parasite drag coefficient with equations from literature
def parasite_drag(MTOM,S):
    d_q_airframe_wing = 1.6*(2.2046*MTOM/1000)**(2/3)                            # ft^2 d/q for fuselage and wing
    d_q_rotors = 0.4*(2.2046*MTOM/1000)**(2/3)                                   # ft^2 d/q for rotors
    cd_interference = 0.003                                                      # [-] cd 
    cd_parasite = 0.0929* (d_q_airframe_wing + d_q_rotors)/S + cd_interference   #conversion to m2, then make it unitless by /S        [-]
    return cd_parasite

# Set up cl and cd for plot
def dragpolar_comp(b,S,h=500,v=const.v_cruise,c=None,Sf=2):
    if c == None:
        c = S/b
    cd_profile = profile_drag(h,v,c,Sf)
    cd_parasite = parasite_drag(160,S)
    cl = np.linspace(-0.4,1.5,150)
    A = b**2/S
    e = 1.78*(1-0.045*A**0.68)-0.64 
    cd_comp = cd_parasite + cd_profile + cl**2/(pi * A * e)
    # print("cd_profile = ",cd_profile)
    # print("cd_parasite =",cd_parasite)
    # print("cd_0 =",cd_parasite+cd_profile)
    return cl, cd_comp



