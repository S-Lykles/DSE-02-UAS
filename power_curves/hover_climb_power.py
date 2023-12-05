import numpy as np
from .inputs import *
from .rotor_tool import *

D_v = 0.04 * W #Vertical drag

T_hover = ((1 + D_v / W) * W) #Hover thrust

HP = ((T_hover * np.sqrt(DL / (2 * rho))) / 550) * HPtoWatt  #Hover power in Watt (not horse power tristan)

v1 = np.sqrt(DL / (2 * rho)) #induced velocity

Delta_P_climb = ((W / 550) * ((vc / 2) + np.sqrt(((vc / 2) ** 2 + (v1 ** 2) - v1) ))) * HPtoWatt #Extra power required to climb at 2m/s in Watt

Clim_P = HP + Delta_P_climb #Power required to climb in watt


#print('hover power', HP)
#print('climb power', Clim_P)
