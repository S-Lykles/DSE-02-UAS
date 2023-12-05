import numpy as np
from matplotlib import pyplot as plt
from math import pi, sqrt

# Calculate parasite drag coefficient with equations from literature
def parasite_drag(MTOM,S):
    d_airframe_wing = 1.6*(2.2046*MTOM/1000)**(2/3)                         # ft^2
    d_rotors = 0.4*(2.2046*MTOM/1000)**(2/3)                                # ft^2
    d_interference = 0.003                                                  # [-]
    cd_parasite = 0.0929* (d_airframe_wing + d_rotors)/S + d_interference   #conversion to m2, then make it unitless by /S        [-]
    return cd_parasite, d_airframe_wing, d_rotors

# Set up cl and cd for plot
def dragpolar_heli(b,S):
    cd_parasite, d_airframe_wing, d_rotors = parasite_drag(160,2)
    cl = np.linspace(-1.5,1.5,150)
    A = b**2/S
    e = 0.78
    cd_heli = cd_parasite + cl**2/(pi * A * e)
    return cl, cd_heli




