import numpy as np
from parameters_weight_estimations import *

# Weight estimation of main wing
W_wing = 0.04674*MTOW**0.397 * S**0.36 * n_ult**0.397 * A**1.712

# Weight estimation of empennage
# W_h = (3.184*MTOW**0.887 *S_h**0.101 *A_h**0.138)/(174.04*t_rh**0.223)
# W_v = (1.68*MTOW**0.567 *S_v**1.249 *A_v**0.482)/(639.95*t_rv**0.747 *(np.cos(chord_sweep_angle))**0.882)
#
# W_emp = W_h + W_v

W_emp = 0.023*MTOW

# Weight estimation of fuselage
W_f = 14.86*MTOW**0.144 *(l/d)**0.778 *l**0.383

# Weight estimation of the nacelle
W_nac = 0.24 * P_cruise_max

# Weight estimation of a non-retractable landing gear

W_m_lg = 0.013*MTOW + 0.362*W_L**0.417 *l_sm**0.183
W_n_lg = 0.0013*MTOW + 0.007157*W_L**0.749 *n_ult_l *l_sn**0.788
# Summation of both nose and main landing gear (decreased
W_lg = 1.24 + W_m_lg + W_n_lg

# Summation to generate a final estimation for the structural weight

W_struc = W_wing + W_emp + W_f + W_nac + W_lg

# Summation to generate estimation for the propulsion system
W_prop = N_electric*(W_electro_motor + W_rotor) + N_gas*(W_gas_motor + W_rotor) + W_fuel_sys + W_battery + W_generator

#Summaion to generate an estimation for the avionics system
W_avionics = W_missioncomputer + W_nav_sys + W_flt_ctrl

OEW = W_prop + W_struc + W_avionics + payload_sup
MTOW = OEW + W_fuel
payload_range = 50 + 160-MTOW
print("The payload range for the supply mission is [50,", payload_range "].")