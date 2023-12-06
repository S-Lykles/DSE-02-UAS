import numpy as np
import matplotlib.pyplot as plt
import const as c
import dual_phase.inputs as i

W = i.W
R = i.R
alpha_m = -9999
rho = i.rho
sigma = N_bl*c*R/(np.pi*R**2)

omega = -9999
V = i.V
D_v = i.D_v
k = i.k
v_i = i.v_i
A_eq = i.A_eq


def power_required_rotor(W, R, alpha_m, rho, sigma, omega, V, D_v, k, v_i, A_eq):


    CD_p = 0.0087 -0.0216*alpha_m+0.4*alpha_m**2
    # CD_p = 0.011 + 0.4*alpha**2
    # CD_p = 0.009 + 0.73*alpha_m**2



    P_phov = sigma*CD_p/8*rho*(omega*R)**3*np.pi*R**2
    mu = V/(omega*R)
    P_p = P_phov * (1+4.65*mu**2)


    k_dl = (1+D_v/W)
    T = k_dl*W
    P_i = k*T*v_i


    P_par = A_eq/2 * rho * V**3

    P_mr = P_p + P_i + P_par

    # P_loss = (3% to 6%) of P_mr
    P_loss = 3/100*P_mr

    P_req_r = P_mr +P_loss
    return P_req_r






def power_required_wing():
    return P_r_w


def fuel_weight(SFC, P):
    W_f = SFC*P
    return W_f

