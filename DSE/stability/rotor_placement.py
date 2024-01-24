import numpy as np
from DSE import const
from DSE.aero.aero_constants import *
# from DSE.structures.center_of_gravity import *

Lambda_LE = np.arctan(0.25*(c_root - c_tip)/(0.5*b)) # sweep leading edge in rad
Lambda_TE = np.arctan(0.75*(c_root - c_tip)/(0.5*b)) # sweep trailing edge in rad

D_vertical = 0.55*2 # placeholder, will change when vertical rotor sizing is done
D_vtol = 0.5*2

taper = c_tip/c_root
MAC = c_bar
y_mac = 0.5 * b * (1 + 2 * taper) / (3 + 3 * taper)

x_LEMAC = const.xlemac
x_wing_LE_rootchord = x_LEMAC #- y_mac * np.tan(Lambda_LE)
print(x_wing_LE_rootchord)

clearance = 0.2 # rotor clearance
b_ht = D_vertical + D_vtol

def rotor_locations(x_wing_LE_rootchord=x_wing_LE_rootchord, b_ht=b_ht, Lambda_LE=Lambda_LE, Lambda_TE=Lambda_TE, D_vtol=D_vtol, clearance=clearance):
    x_PF = x_wing_LE_rootchord + 0.5*b_ht*np.tan(Lambda_LE) - (D_vtol/2 + clearance)/np.cos(Lambda_LE)
    x_PR = x_wing_LE_rootchord + c_root +0.5*b_ht*np.tan(Lambda_TE) + (D_vtol/2 + clearance)/np.cos(Lambda_TE)
    return x_PF, x_PR

Rotor_Front_X, Rotor_Rear_X = rotor_locations()
print(Rotor_Front_X, Rotor_Rear_X , x_wing_LE_rootchord, x_LEMAC)

print('R1_xcg, R2_xcg',rotor_locations())
print('R_xcg mean',(rotor_locations()[0]+rotor_locations()[1])/2)
