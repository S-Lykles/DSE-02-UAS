import numpy as np

def EoM():
    "Moment equations around the center of gravity"
    "Ct_p = the propellor thrust coefficient"
    "Ct_fr = the rotor thrust coefficient of the front two rotors"
    "Ct_aft = the rotor thrust coefficient of the aft two rotors"
    Cm = Cm_w + Cm_fus + Cm_h + (
                CL_w * h_ac - CD_w * l_ac + (CL_h * h_h + CD_h * l_h) * (S_h / S) * (V_h / V) ** 2) * np.sin(
        alpha * deg2rad) + (
                     CL_w * l_ac + CD_w * h_ac + (CD_h * h_h - CL_h * l_h) * (S_h / S) * (V_h / V) ** 2) * np.cos(
        alpha * deg2rad) - Ct_p * h_p + 2 * Ct_fr * l_fr - 2 * Ct_aft * l_aft

