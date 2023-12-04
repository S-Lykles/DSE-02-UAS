import numpy as np
from .inputs import *
from .rotor_tool import *

D_v = 0.04 * W

T_hover = ((1 + D_v / W) * W)

HP = ((T_hover * np.sqrt(DL / (2 * rho))) / 550) * HPtoWatt

v1 = np.sqrt(DL / (2 * rho))

Delta_P_climb = (W / 550) * ((vc / 2) + np.sqrt((vc / 2) ** 2 + v1 ** 2) - v1) * HPtoWatt

Clim_P = HP + Delta_P_climb

print('hover thrust', T_hover)
print('hover power', HP)
print('delta p for climb', Delta_P_climb)
print('climb power', Clim_P)

