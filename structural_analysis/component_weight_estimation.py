import numpy as np


# Class 2 weight estimation dual phase using cessna method from Roskam book:
def class_two_dual_phase(MTOW, l_sm, l_sn, W_rotor):
    # Weight estimation of main wing in kg(full cantilever wing assumed)
    lb_to_kg = 2.20462262
    W_wing = (0.04674 * MTOW**0.397 * S**0.36 * n_ult**0.397 * A**1.712) / lb_to_kg

    W_emp = 0.023 * MTOW / lb_to_kg

    # Weight estimation of fuselage (high wing configuration assumed)
    W_f = (14.86*MTOW**0.144 *(l/d)**0.778 *l**0.383) / lb_to_kg

    # Weight estimation of the nacelle (horizontally opposed engines assumed)
    W_nac = (0.24 * P_cruise_max) / lb_to_kg

    # Weight estimation of a non-retractable landing gear
    W_m_lg = (0.013*MTOW + 0.362*W_L**0.417 * l_sm**0.183) / lb_to_kg
    W_n_lg = (0.0013*MTOW + 0.007157*W_L**0.749 * l_sn**0.788) / lb_to_kg

    # Summation of both nose and main landing gear (decreased constant factor by factor 5 as design will be lighter)
    W_lg = 1.24 + W_m_lg + W_n_lg

    # Summation to generate a final estimation for the structural weight
    W_struc = W_wing + W_emp + W_f + W_nac + W_lg

    # Summation to generate estimation for the propulsion system
    W_prop = (N_electric * (W_electro_motor + W_rotor) + N_gas * (W_gas_motor + W_rotor)
              + W_fuel_sys + W_battery + W_generator)

    # Summation to generate an estimation for the avionics system
    W_avionics = W_missioncomputer + W_nav_sys + W_flt_ctrl

    OEW = W_prop + W_struc + W_avionics + payload_sup
    MTOW = OEW + W_fuel
    payload_range = 50 + 160 - MTOW
    print("The payload range of the dual phase configuration for the supply mission is [50,", payload_range, "].")
    return OEW



# Class 2 weight estimation compound helicopter using "Weight estimation of helicopter design":

def class_two_compound_helicopter(R, C, N_blades, V_tip, t_avg, MTOW, S_ch, A_ch, rpm):
    # Rotor weight estimation
    lb_to_kg = 2.20462262
    kw_to_hp = 1.341022
    W_R = (3.45 * 10**(-4) * (R * C * N_blades * V_tip**2 * (t_avg + 0.21))**0.89) / lb_to_kg

    # Tail weight estimation
    # Tail rotor weight estimation
    W_tr = 7.4 * 10 ** (-4) * MTOW / lb_to_kg

    # Horizontal tail stabilizer weight estimation
    W_hs = 6.9 * 10 ** (-4) * MTOW ** 1.2 / lb_to_kg

    W_tail = W_tr + W_hs

    # Weight estimation of the main body (internal cargo assumed)
    W_body = 0.058 * R * MTOW ** 0.67 * n_ult ** 0.335 / lb_to_kg

    # Weight estimation landing gear and support structure (skid-type assumed)
    W_lg = 0.03 * MTOW / lb_to_kg

    # Weight estimation of the wing (full cantilever wings assumed)
    W_wing = 0.04674 * MTOW ** 0.397 * S_ch ** 0.36 * n_ult ** 0.397 * A_ch ** 1.712 / lb_to_kg


    W_struc = W_tr + W_hs + W_tail + W_body + W_lg + W_wing

    # Weight estimation of propulsion group (shaft-driven assumed)
    # Weight estimation of engine
    W_motor = 36.4 * P_hov_max ** 0.31 / lb_to_kg

    # Weight estimation of main rotor drive system
    W_drivetrain = 0.2 * (P_hov_max * kw_to_hp * 5250/rpm) ** 0.787 / lb_to_kg

    # Weight estimation of propulsion system accessories
    W_prop_acc = 0.52 * W_motor / lb_to_kg

    # Weight estimation of tail rotor drive system
    W_trd = 2.46 * 10 ** (-3) * R ** 3.248 / lb_to_kg


    W_prop_sys = W_motor + W_drivetrain + W_prop_acc + W_trd

    W_avionics = W_missioncomputer + W_nav_sys + W_flt_ctrl

    OEW = W_prop_sys + W_avionics + W_R + W_struc
    MTOW = OEW + W_fuel
    payload_range = 50 + 160-MTOW
    print("The payload range of the compound helicopter for the supply mission is [50,", payload_range, "].")
    pass

# Class 2 weight estimation for the tilt-wing configuration using cessna method from Roskam book:
# Notes: did the weight of the main wing twice without a credible source, still requires a source!!


def class_two_tilt_wing(MTOW, l_sm, l_sn):
    # Weight estimation of main wing in kg(full cantilever wing assumed)
    lb_to_kg = 2.20462262
    W_wing = (0.04674 * MTOW ** 0.397 * S ** 0.36 * n_ult ** 0.397 * A ** 1.712) / lb_to_kg

    # Weight estimation of empennage (traditional empennage assumed)
    # W_h = (3.184*MTOW**0.887 *S_h**0.101 *A_h**0.138)/(174.04*t_rh**0.223) /2.20462262
    # W_v = (1.68*MTOW**0.567 *S_v**1.249 *A_v**0.482)/(639.95*t_rv**0.747 *(np.cos(chord_sweep_angle))**0.882)/
    # 2.20462262
    #
    # W_emp = W_h + W_v

    W_emp = 0.023 * MTOW / lb_to_kg

    # Weight estimation of fuselage (high wing configuration assumed)
    W_f = (14.86 * MTOW ** 0.144 * (l / d)**0.778 * l**0.383) / lb_to_kg

    # Weight estimation of the nacelle (horizontally opposed engines assumed)
    W_nac = (0.24 * P_cruise_max) / lb_to_kg

    # Weight estimation of a non-retractable landing gear
    W_m_lg = (0.013 * MTOW + 0.362 * W_L**0.417 * l_sm**0.183) / lb_to_kg
    W_n_lg = (0.0013 * MTOW + 0.007157 * W_L**0.749 * l_sn**0.788) / lb_to_kg

    # Summation of both nose and main landing gear (decreased constant factor by factor 5 as design will be lighter)
    W_lg = 1.24 + W_m_lg + W_n_lg

    # Weight estimation of mechanical installation for tilting the wing
    # W_mech =

    # Summation to generate a final estimation for the structural weight
    W_struc = 2*W_wing + W_emp + W_f + W_nac + W_lg #+ W_mech

    # Summation to generate estimation for the propulsion system
    W_prop = N_electric * (W_electro_motor + W_rotor) + N_gas * (
                W_gas_motor + W_rotor) + W_fuel_sys + W_battery + W_generator

    # Summation to generate an estimation for the avionics system
    W_avionics = W_missioncomputer + W_nav_sys + W_flt_ctrl

    OEW = W_prop + W_struc + W_avionics + payload_sup
    MTOW = OEW + W_fuel
    payload_range = 50 + 160 - MTOW
    print("The payload range of the tilt-wing configuraiton for the supply mission is [50,", payload_range, "].")
    pass

# Class 2 weight estimation for the tailsitter configuration based cessna method from Roskam book:
#Notes: empennage/landing gear more interesting, need input for those
def class_2_tailsitter(MTOW, ):
    # Weight estimation of main wing in kg(full cantilever wing assumed)
    lb_to_kg = 2.20462262
    W_wing = (0.04674 * MTOW ** 0.397 * S ** 0.36 * n_ult ** 0.397 * A ** 1.712) / lb_to_kg

    # Weight estimation of fuselage (high wing configuration assumed)
    W_f = (14.86 * MTOW ** 0.144 * (l / d) ** 0.778 * l ** 0.383) / lb_to_kg

    # Weight estimation of the nacelle (horizontally opposed engines assumed)
    W_nac = (0.24 * P_cruise_max) / lb_to_kg

    # Weight estimation of a non-retractable landing gear using a comparable design
    W_lg = 0.02 * MTOW / lb_to_kg

    # Summation to generate a final estimation for the structural weight
    W_struc = W_wing + W_f + W_nac + W_lg

    # Summation to generate estimation for the propulsion system
    W_prop = N_electric * (W_electro_motor + W_rotor) + N_gas * (
                W_gas_motor + W_rotor) + W_fuel_sys + W_battery + W_generator

    # Summation to generate an estimation for the avionics system
    W_avionics = W_missioncomputer + W_nav_sys + W_flt_ctrl

    OEW = W_prop + W_struc + W_avionics + payload_sup
    MTOW = OEW + W_fuel
    payload_range = 50 + 160 - MTOW
    print("The payload range of the tailsitter configuration for the supply mission is [50,", payload_range, "].")
    pass