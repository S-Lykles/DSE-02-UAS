import numpy as np
from DSE.aero.aero_constants import *
from DSE.structures.center_of_gravity import *

Lambda_LE = np.arctan(0.25*(c_root - c_tip)/(0.5*b)) # sweep leading edge in rad
Lambda_TE = np.arctan(0.75*(c_root - c_tip)/(0.5*b)) # sweep trailing edge in rad

D_vertical = 1 # placeholder, will change when vertical rotor sizing is done
D_vtol = 0.5
taper = 0.4
MAC = c_root * (2/3) * ((1+taper+taper**2)/(1+taper))
y_mac = 0.5 * b * (1 + 2 * taper) / (3 + 3 * taper)

x_cg = 1.5 #class_two_cg_estimation(True, False, False, False)[1][0]
print('xcg', x_cg)
x_MAC_025c = x_cg - 0.4
x_wing_LE_rootchord = x_MAC_025c - 0.25*MAC - y_mac * np.tan(Lambda_LE)
print(x_wing_LE_rootchord)

x_wing_MAC = 1 #most front point of the wing wrt datum
clearance = 0.05 # rotor clearance with fuselage
b_ht = D_vertical + D_vtol

def rotor_locations(x_wing_LE_rootchord, b_ht, Lambda_LE, Lambda_TE, D_vtol, clearance):
    x_PF = x_wing_LE_rootchord + 0.5*b_ht*np.tan(Lambda_LE) - (D_vtol/2 + clearance)/np.cos(Lambda_LE)
    x_PR = x_wing_LE_rootchord + c_root + 0.5*b_ht*np.tan(Lambda_TE) + (D_vtol/2 + clearance)/np.cos(Lambda_TE)
    return x_PF, x_PR
print(rotor_locations(x_wing_LE_rootchord, b_ht, Lambda_LE, Lambda_TE, D_vtol, clearance))
