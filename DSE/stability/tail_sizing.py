import numpy as np
from const import T, gamma, R

# default value
T = 288.15
gamma = 1.4
R = 287
def tail_sizing(tail_config, V, R, gamma, T):

    SM = 0.05 #PLACEHOLDER


    dL_A_h = CL_alpha_A_h *da/2*rho*V**2*S
    dL_h = CL_alpha_A_h * (da-de)/2*rho*V_h**2*S_h

    C_m = C_m0 + (CL_h + CL_A_h*S_h*V_h/(S*V))*(x_cg-x_n)

    x_np_bar = x_ac_bar +CL_alpha_h/C_L_alpha_A_h * (1-de/da) * (S_h * l_h)/(S*c_bar)*Vh_V_2

    if tail_config == 'uselage-mounted':
        Vh_V_2 = 0.85
    elif tail_config == 'uselage-mounted':
        Vh_V_2 = 0.95
    elif tail_config == 'uselage-mounted'
        Vh_V_2 = 1
    else:
        print('Check the values of tail configuration')

    A_h = b_h/S_h
    V_h = np.sqrt(Vh_V_2)*V
    a = np.sqrt(R*gamma*T)
    M = V_h/a
    beta = np.sqrt(1-M**2)

    CL_alpha_h = 2*np.pi*A_h/(2+np.sqrt(4+(A_h*beta/eta)**2*(1+np.tan(sweep_05_Ch_rad)**2/beta**2)))
    CL_alpha_A_h = Cl_alpha_w * (1+2.15*b_f/b)*S_net/S+np.pi/2*b_f**2/S

    x_cg_bar = x_ac_bar + Cl_alpha_h/Cl_alpha_A_h*(1-de/da)*S_h*l_h/(S*c_bar)*Vh_V_2-SM

    Sh_S = 1/(CL_alpha_h/CL_alpha_A_h*(1-de/da)*l_h/c_bar*Vh_V_2)*x_cg_bar-(x_ac_bar-SM)/(CL_alpha_h/CL_alpha_A_h*(1-de/da)*l_h/c_bar*Vh_V_2)
    S_h = S*Sh_S
