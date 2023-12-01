import numpy as np
import matplotlib.pyplot as plt
from inputs import *



def rotor_sizing_tool(DL, N):
    #rotor sizing
    R               = np.sqrt(W/N/DL/pi)
    V_tip           = 150*(2*R)**0.171
    D_v             = 0.04*W
    k_dl            = 1 + D_v/W
    omega           = V_tip/R
    Vne             = 1.1* V_max
    mu_Vne          = 1.1*Vne/(omega*R)
    Advance_ratio   = V_max / V_tip

    #level flight
    T_level         = W*k_dl
    C_T_level       = T_level/ (rho*pi*R**2*omega**2*R**2)
    sig_level       = C_T_level/C_T_sig

    #turning flight
    n_z             = 1 / np.cos(psi_rad)
    T_turn          = W * k_dl * n_z
    C_T_turn        = T_turn / (rho * pi * R**2 * (omega*R)**2)
    sig_turn        = C_T_turn / C_T_sig

    # T_gust = n_z*k_dl*
    # C_T_gust = T_gust/ (Vne*pi*R**2*omega**2*R**2)

    sig_max = max(sig_level, sig_turn)

    return R, D_v, omega, T_level, sig_max




def generate_Preq_rotor(A_eq, R, D_v, omega, T_level, sig_max, t_start, t_end, step):
    v = np.linspace(t_start, t_end, step)

    #Profile drag
    advanced_ratio = v / (omega*R)

    C_t = W / (rho * np.pi * R**2 * (omega*R)**2)
    Cl_bar = 6.6*(C_t / sig_max)
    alpha_m = Cl_bar / Cl_alpha_rot

    CD_p1 = 0.0087 - 0.0216*alpha_m + (0.4 * alpha_m **2)
    CD_p2 = 0.011 + 0.4*(alpha_m**2)
    CD_p3 = 0.009 + 0.73*(alpha_m**2)
    C_D_p = (CD_p1 + CD_p3 + CD_p2) / 3

    P_hov = (1/8)*sig_max*C_D_p*rho*(omega*R)**3*np.pi*(R**2)

    P_profile_drag_arr = P_hov * (1 + 4.65*(advanced_ratio**2))

    #induced drag power
    v_ih = np.sqrt(W / (2*rho*np.pi*R**2))
    v_ibar = 1 / v
    v_i = v_ibar * v_ih
    P_induced_1_arr = k * T_level * v_i

    #parsite power
    P_parasite_arr = 0.5 * rho * A_eq * (v**3)

    #power loss
    P_loss_arr = 0.04*(P_profile_drag_arr + P_induced_1_arr + P_parasite_arr)

    #total power
    P_tot_req_level_arr = P_profile_drag_arr + P_induced_1_arr + P_parasite_arr + P_loss_arr

    return P_tot_req_level_arr, v

