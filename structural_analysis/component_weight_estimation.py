import numpy as np
from parameters_weight_estimations import *
def class_two_dual_phase( ):
    #Weight estimation dual phase:
    # Weight estimation of main wing in kg(full cantilever wing assumed)
    W_wing = (0.04674*MTOW**0.397 * S**0.36 * n_ult**0.397 * A**1.712)/2.20462262

    # Weight estimation of empennage (traditional empennage assumed)
    # W_h = (3.184*MTOW**0.887 *S_h**0.101 *A_h**0.138)/(174.04*t_rh**0.223) /2.20462262
    # W_v = (1.68*MTOW**0.567 *S_v**1.249 *A_v**0.482)/(639.95*t_rv**0.747 *(np.cos(chord_sweep_angle))**0.882)/ 2.20462262
    #
    # W_emp = W_h + W_v

    W_emp = 0.023*MTOW /2.20462262

    # Weight estimation of fuselage (high wing configuration assumed)
    W_f = (14.86*MTOW**0.144 *(l/d)**0.778 *l**0.383)/2.20462262

    # Weight estimation of the nacelle (horizontally opposed engines assumed)
    W_nac = (0.24 * P_cruise_max)/2.20462262

    # Weight estimation of a non-retractable landing gear
    W_m_lg = (0.013*MTOW + 0.362*W_L**0.417 *l_sm**0.183)/2.20462262
    W_n_lg = (0.0013*MTOW + 0.007157*W_L**0.749 *n_ult_l *l_sn**0.788)/2.20462262

    # Summation of both nose and main landing gear (decreased constant factor by factor 5 as design will be lighter)
    W_lg = 1.24 + W_m_lg + W_n_lg

    # Summation to generate a final estimation for the structural weight
    W_struc = W_wing + W_emp + W_f + W_nac + W_lg

    # Summation to generate estimation for the propulsion system
    W_prop = N_electric*(W_electro_motor + W_rotor) + N_gas*(W_gas_motor + W_rotor) + W_fuel_sys + W_battery + W_generator

    #Summation to generate an estimation for the avionics system
    W_avionics = W_missioncomputer + W_nav_sys + W_flt_ctrl

    OEW = W_prop + W_struc + W_avionics + payload_sup
    MTOW = OEW + W_fuel
    payload_range = 50 + 160-MTOW
    print("The payload range of the dual phase configuration for the supply mission is [50,", payload_range "].")
    return



#Weight estimation compound helicopter
#Rotor weight estimation

W_R = (3.45*10**(-4)* (R*C*N_blades* V_tip**2 *(t_avg + 0.21))**0.89)/2.20462262

#Tail weight estimation
#Tail rotor weight estimation
W_tr = 7.4*10**(-4) * MTOW /2.20462262

#Horizontal tail stabilizer weight estimation
W_hs = 6.9*10**(-4)* MTOW**1.2 /2.20462262

W_tail = W_tr + W_hs

#Weight estimation of the main body (internal cargo assumed)
W_body = 0.058*R*MTOW**0.67*n_ult**0.335 /2.20462262

#Weight estimation landing gear and support structure (skid-type assumed)
W_LG = 0.03 *MTOW /2.20462262

#Weight estimation of the wing (full cantilever wings assumed)
W_wing = 0.04674*MTOW**0.397 * S_ch**0.36 * n_ult**0.397 * A_ch**1.712 /2.20462262


W_struc = W_tr + W_hs + W_tail + W_body + W_LG + W_wing

#Weight estimation of propulsion group (shaft-driven assumed)
#Weight estimation of engine
W_E = 36.4 * P_hov_max**0.31 /2.20462262

#Weight estimation of main rotor drive system
W_T = 0.2 * (P_hov_max*5250/rpm)**0.787 /2.20462262

#Weight estimation of propulsion system accessories
W_PA = 0.52*W_E /2.20462262

#Weight estimation of tail rotor drive system
W_trd = 2.46*10**(-3)*R**3.248 /2.20462262


W_prop_sys = W_E + W_T + W_PA + W_trd

W_avionics = W_missioncomputer + W_nav_sys + W_flt_ctrl

OEW = W_prop_sys + W_avionics + W_R + W_struc
MTOW = OEW + W_fuel
payload_range = 50 + 160-MTOW
print("The payload range of the compound helicopter for the supply mission is [50,", payload_range "].")