import numpy as np
from inputs import *
from main import *
from rotor_tool import *


Dv = 0.04 * W

T_hover = (1 + D_v/W) * W

HP = (T_hover * np.sqrt(DL / 2*rho)) / 550

print(T_hover)
print(HP * HPtoWatt)

