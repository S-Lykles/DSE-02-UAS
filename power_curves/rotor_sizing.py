from Inputs_Preq_rotorcraft import *



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
    n_z             = 1 / cos(psi_rad)
    T_turn          = W * k_dl * n_z
    C_T_turn        = T_turn / (rho * pi * R**2 * (omega*R)**2)
    sig_turn        = C_T_turn / C_T_sig


    # T_gust = n_z*k_dl*
    # C_T_gust = T_gust/ (Vne*pi*R**2*omega**2*R**2)

    sig_max = max(sig_level, sig_turn)

    return R, D_v, omega, T_level, sig_max