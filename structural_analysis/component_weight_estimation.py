import numpy as np
from parameters_weight_estimations import *

# Class 2 weight estimation dual phase using cessna method from Roskam book:
def class_two_dual_phase(MTOW):
    # Weight estimation of main wing in kg(full cantilever wing assumed)
    kw_to_hp = 1.341022
    lb_to_kg = 2.20462262
    m_to_ft = 0.3048
    W_wing = (0.04674 * MTOW**0.397 * S**0.36 * n_ult**0.397 * AR**1.712) / lb_to_kg

    W_emp = 0.023 * MTOW / lb_to_kg

    # Weight estimation of fuselage (high wing configuration assumed)
    W_f = (14.86*MTOW**0.144 *(l/perimeter)**0.778 *l**0.383) / lb_to_kg

    # Weight estimation of the nacelle (horizontally opposed engines assumed)
    W_nac = (0.24 * P_cruise_max* kw_to_hp) / lb_to_kg

    # Weight estimation of a non-retractable landing gear
    W_m_lg = (0.013*MTOW + 0.362*W_L**0.417 * l_sm**0.183) / lb_to_kg
    W_n_lg = (0.0013*MTOW + 0.007157*W_L**0.749 * l_sn**0.788) / lb_to_kg

    # Summation of both nose and main landing gear (decreased constant factor by factor 5 as design will be lighter)
    W_lg = 1.24 + W_m_lg + W_n_lg

    # Summation to generate a final estimation for the structural weight
    W_struc = W_wing + W_emp + W_f + W_nac + W_lg

    # Summation to generate estimation for the propulsion system
    W_prop = (20 + N_electric * W_rotor_electric + N_gas * W_rotor_gas + W_generator + W_fuel_sys)
               #+ W_fuel_sys + W_battery)

    # Summation to generate an estimation for the avionics system
    W_avionics = W_missioncomputer + W_nav_sys + W_flt_ctrl

    OEW = W_prop + W_struc + W_avionics / lb_to_kg + payload_sup / lb_to_kg
    MTOW = OEW + W_fuel / lb_to_kg
    payload_range = 50 + 160 - MTOW

    print("The weight for the structures is:", W_struc)
    print("The weight for the propulsion subsystem is:", W_prop)
    print("The weight for the avionics subsystem is:", W_avionics/lb_to_kg)
    print("The weight for the fuel is:", W_fuel / lb_to_kg)
    print("The weight for the payload is:", payload_sup/lb_to_kg)
    print(MTOW)
    return #W_prop, W_avionics/lb_to_kg, payload_sup/lb_to_kg, W_struc, W_fuel/lb_to_kg, OEW, MTOW



# Class 2 weight estimation compound helicopter using "Weight estimation of helicopter design":

def class_two_compound_helicopter(MTOW):
    lb_to_kg = 2.20462262
    kw_to_hp = 1.341022
    m_to_ft = 0.3048

    # Estimation of the empennage weight
    W_emp = 0.023 * MTOW / lb_to_kg

    # Weight estimation of the main body (internal cargo assumed)
    W_body = 0.058 * R * MTOW ** 0.67 * n_ult ** 0.335 / lb_to_kg

    # Weight estimation landing gear and support structure (skid-type assumed)
    W_lg = 0.03 * MTOW / lb_to_kg

    # Weight estimation of the wing (full cantilever wings assumed)
    W_wing = 0.04674 * MTOW ** 0.397 * S_ch ** 0.36 * n_ult ** 0.397 * AR_ch ** 1.712 / lb_to_kg


    W_struc = W_emp + W_body + W_lg + W_wing

    # Weight estimation of propulsion group
    W_R = (9.56 * 10 ** (-4) * (((sigma * MTOW * V_tip ** 2) / (DL)) * (sigma * R * t_c) / (N_blades) + 0.067) ** 0.89) / lb_to_kg

    W_prop = (20 + N_electric * W_rotor_electric + W_R + W_generator + W_fuel_sys)

    # # Weight estimation of engine
    # W_motor = 36.4 * (P_hov_max * kw_to_hp) ** 0.31 / lb_to_kg
    # # Weight estimation of main rotor drive system
    # W_drivetrain = 0.2 * (P_hov_max * kw_to_hp * 5250/rpm) ** 0.787 / lb_to_kg
    # # Weight estimation of propulsion system accessories
    # W_prop_acc = 0.52 * W_motor
    # # Weight estimation of tail rotor drive system
    # W_trd = 2.46 * 10 ** (-3) * R ** 3.248 / lb_to_kg
    # W_prop_sys = W_motor + W_drivetrain + W_prop_acc + W_trd

    W_avionics = (W_missioncomputer + W_nav_sys + W_flt_ctrl) / lb_to_kg

    OEW = W_prop + W_avionics + W_struc + payload_sup / lb_to_kg
    MTOW = OEW + W_fuel /lb_to_kg
    payload_range = 50 + 160-MTOW

    print("The weight for the structures is:", W_struc)
    print("The weight for the propulsion subsystem is:", W_prop)
    print("The weight for the avionics subsystem is:", W_avionics)
    print("The weight for the fuel is:", W_fuel / lb_to_kg)
    print("The weight for the payload is:", payload_sup/lb_to_kg)
    return MTOW

# Class 2 weight estimation for the tilt-wing configuration using cessna method from Roskam book:
# Notes: did the weight of the main wing twice without a credible source, still requires a source!!


def class_two_tilt_wing(MTOW):
    # Weight estimation of main wing in kg(full cantilever wing assumed)
    lb_to_kg = 2.20462262
    kw_to_hp = 1.341022
    W_wing = (0.04674 * MTOW ** 0.397 * S ** 0.36 * n_ult ** 0.397 * AR ** 1.712) / lb_to_kg

    # Weight estimation of empennage (traditional empennage assumed)
    # W_h = (3.184*MTOW**0.887 *S_h**0.101 *A_h**0.138)/(174.04*t_rh**0.223) /2.20462262
    # W_v = (1.68*MTOW**0.567 *S_v**1.249 *A_v**0.482)/(639.95*t_rv**0.747 *(np.cos(chord_sweep_angle))**0.882)/
    # 2.20462262
    #
    # W_emp = W_h + W_v

    W_emp = 0.023 * MTOW / lb_to_kg

    # Weight estimation of fuselage (high wing configuration assumed)
    W_f = (14.86 * MTOW ** 0.144 * (l / perimeter)**0.778 * l**0.383) / lb_to_kg

    # Weight estimation of the nacelle (horizontally opposed engines assumed)
    W_nac = (0.24 * P_cruise_max*kw_to_hp) / lb_to_kg

    # Weight estimation of a non-retractable landing gear
    W_m_lg = (0.013 * MTOW + 0.362 * W_L**0.417 * l_sm**0.183) / lb_to_kg
    W_n_lg = (0.0013 * MTOW + 0.007157 * W_L**0.749 * l_sn**0.788) / lb_to_kg

    # Summation of both nose and main landing gear (decreased constant factor by factor 5 as design will be lighter)
    W_lg = 1.24 + W_m_lg + W_n_lg

    # Weight estimation of mechanical installation for tilting the wing
    # W_mech =

    # Summation to generate a final estimation for the structural weight
    W_struc = (2*W_wing + W_emp + W_f + W_nac + W_lg)/lb_to_kg #+ W_mech

    # Summation to generate estimation for the propulsion system
    W_prop = (N_electric * (W_electro_motor + W_rotor_electric) + N_gas * (W_gas_motor + W_rotor_gas)
              + W_fuel_sys + W_battery + W_generator)

    # Summation to generate an estimation for the avionics system
    W_avionics = (W_missioncomputer + W_nav_sys + W_flt_ctrl) / lb_to_kg

    OEW = W_prop + W_avionics + W_struc + payload_sup / lb_to_kg
    MTOW = OEW + W_fuel /lb_to_kg
    payload_range = 50 + 160 - MTOW

    print("The weight for the structures is:", W_struc)
    print("The weight for the propulsion subsystem is:", W_prop)
    print("The weight for the avionics subsystem is:", W_avionics)
    print("The weight for the fuel is:", W_fuel /lb_to_kg)
    print("The weight for the payload is:", payload_sup)
    return MTOW

# Class 2 weight estimation for the tailsitter configuration based cessna method from Roskam book:
#Notes: empennage/landing gear more interesting, need input for those
def class_2_tailsitter(MTOW, ):
    # Weight estimation of main wing in kg(full cantilever wing assumed)
    lb_to_kg = 2.20462262
    kw_to_hp = 1.341022
    W_wing = (0.04674 * MTOW ** 0.397 * S ** 0.36 * n_ult ** 0.397 * A ** 1.712) / lb_to_kg

    # Weight estimation of fuselage (high wing configuration assumed)
    W_f = (14.86 * MTOW ** 0.144 * (l / d) ** 0.778 * l ** 0.383) / lb_to_kg

    # Weight estimation of the nacelle (horizontally opposed engines assumed)
    W_nac = (0.24 * P_cruise_max*kw_to_hp) / lb_to_kg

    # Weight estimation of a non-retractable landing gear using a comparable design
    W_lg = 0.02 * MTOW / lb_to_kg

    # Summation to generate a final estimation for the structural weight
    W_struc = W_wing + W_f + W_nac + W_lg

    # Summation to generate estimation for the propulsion system
    W_prop = (N_electric * (W_electro_motor + W_rotor_electric) + N_gas * (W_gas_motor + W_rotor_gas)
              + W_fuel_sys + W_battery + W_generator)

    # Summation to generate an estimation for the avionics system
    W_avionics = W_missioncomputer + W_nav_sys + W_flt_ctrl

    OEW = W_prop + W_struc + W_avionics + payload_sup
    MTOW = OEW + W_fuel
    payload_range = 50 + 160 - MTOW
    print("The payload range of the tailsitter configuration for the supply mission is [50,", payload_range, "].")
    pass