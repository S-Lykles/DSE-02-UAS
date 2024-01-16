import numpy as np
from DSE import const
from DSE.Locations import locations
from DSE.aero import aero_constants



rho = const.rho0
S = aero_constants.S
b = aero_constants.b
c_bar = aero_constants.c_bar
m = const.total_mass

Ixx = -9999 # placeholder, input from structures
Iyy = -9999 # placeholder, input from structures
Ixz = -9999 # placeholder, input from structures
# Kxz = -9999
Jxy = Ixz/(m*b**2)
Ky_2 = Iyy/(m*c_bar**2)
mu_c = m/(rho*S*c_bar)
mu_b = m/(rho*S*b)
Kx_2 = Ixx/(m*b**2)

T1 = -9999 # placeholder, input from propulsion
T2 = -9999 # placeholder, input from propulsion
T3 = -9999 # placeholder, input from propulsion
T4 = -9999 # placeholder, input from propulsion
Lw = -9999 # placeholder, input from propulsion
Lh = -9999 # placeholder, input from propulsion
q_rad= -9999 # placeholder, input for control

l_fr, l_aft, l_acw,l_h,h_p,h_acw,h_h,z_h,X_lemac, Xcg, Zac, Zh = locations()

vtol=True
if vtol:
    CXalpha = 0
    CZalpha = 0
    Cmalpha = 0
    CXalphadott = 0
    CZalphadott = 0
    Cmalphadott = 0
    CZq = (T1+T2+T3+T4)*np.sin(q_rad)
    Cnr = -9999
    Cmq = -9999
    CYr = 0
    Clr = -9999

else:
    CXalpha = - CD_alpha
    CZalpha = - aero_constants.CL_alpha_wing - aero_constants.Cl_alpha_h 
    Cmalpha = aero_constants.CL_alpha_wing * l_acw / aero_constants.c - aero_constants.Cl_alpha_h * l_h / aero_constants.c - CD_alpha_w * Zac / aero_constants.c   
    CXalphadott = 0
    CZalphadott = - CL_alphadott_w - CL_alphadott_h
    Cmalphadott = - CL_alphadott_w * l_acw / c- CL_alphadott_h * l_h / aero_constants.c - CDalphadott * Zac / c
    CZq = (Lw+Lh)*np.sin(q_rad)
    Cnr = -9999
    Cmq = -9999
    CYr = -9999
    Clr = -9999




# Db = -9999
# Dc = -9999
