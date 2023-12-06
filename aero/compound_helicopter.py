import numpy as np
from math import pi, sqrt

# Calculate parasite drag coefficient with equations from literature
def parasite_drag(MTOM,S):
    d_airframe_wing = 1.6*(2.2046*MTOM/1000)**(2/3)                         # ft^2
    d_rotors = 0.4*(2.2046*MTOM/1000)**(2/3)                                # ft^2
    d_interference = 0.003                                                  # [-]
    cd_parasite = 0.0929* (d_airframe_wing + d_rotors)/S + d_interference   #conversion to m2, then make it unitless by /S        [-]
    return cd_parasite, d_airframe_wing, d_rotors, d_interference

# Set up cl and cd for plot
def dragpolar_heli(b,S,Cl_start=0.2, Cl_end=0.8, Cl_step=1000):
    cd_parasite, d_airframe_wing, d_rotors, d_interference = parasite_drag(160,2)
    cl = np.linspace(Cl_start,Cl_end, Cl_step)
    A = b**2/S
    e = 0.78
    cd_heli = cd_parasite + cl**2/(pi * A * e)
    return cl, cd_heli

def v_stall(w, s, rho, cl):
    v = np.sqrt(w*2/(s*rho*cl))
    return v



