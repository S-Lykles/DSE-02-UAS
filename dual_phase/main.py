import numpy as np
import matplotlib.pyplot as plt
import const as c
import .inputs as i


def power_required_rotor():
    CD_p = 0.0087 -0.0216*alpha_m+0.4*alpha_m**2
    # CD_p = 0.011 + 0.4*alpha**2
    # CD_p = 0.009 + 0.73*alpha_m**2

    P_phov = sigma*CD_p/8*rho*(omega*R)**3*np.pi*R**2
    mu = V/(omega*R)
    P_p = P_phov * (1+4.65*mu**2)

    k_dl
    T = k_dl*W
    P_i = k*T*v_i

    P_req_r = P_p + P_i + P_par +P_loss
    return P_req_r






def power_required_wing():
    return P_r_w


def fuel_weight(SFC, P):
    W_f = SFC*P
    return W_f
