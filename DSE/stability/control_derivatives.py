import numpy as np
from DSE import const
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

vtol=True
if vtol:
    CZq = (T1+T2+T3+T4)*np.sin(q_rad)
    Cnr = -9999
    Cmq = -9999
    CYr = 0
    Clr = -9999

else:
    CZq = (Lw+Lh)*np.sin(q_rad)
    Cnr = -9999
    Cmq = -9999
    CYr = -9999
    Clr = -9999




# Db = -9999
# Dc = -9999
