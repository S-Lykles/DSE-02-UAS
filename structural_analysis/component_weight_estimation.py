import numpy as np
from parameters_weight_estimations import *

# Weight estimation of main wing
W_wing = 0.04674*MTOW**0.397 * S**0.36 * n_ult**0.397 * A**1.712

# Weight estimation of empennage
W_h = (3.184*MTOW**0.887 *S_h**0.101 *A_h**0.138)/(174.04*t_rh**0.223)
W_v = (1.68*MTOW**0.567 *S_v**1.249 *A_v**0.482)/(639.95*t_rv**0.747 *(np.cos(chord_sweep_angle))**0.882)

W_emp = W_h + W_v

# Weight estimation of fuselage
W_f = 14.86*MTOW**0.144 *(l/d)**0.778 *l**0.383

# Weight estimation of the nacelle
W_n = 0.24 * P_max
