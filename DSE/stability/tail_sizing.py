import numpy as np
from const import T, gamma, R

# default values
T = 288.15
gamma = 1.4
R = 287
eta = 0.95
dz_h = z_t-z_w
# A_h = b_h/S_h
def tail_sizing(S, b, c_bar, b_f, h_f, l_fn, sweep_ang_14_c_rad, CL_alpha_w, S_net, eta, A_h, l_h,dz_h, sweep_ang_015_Ch_rad, c_root, c_tip, tail_config, V, R, gamma, T, rho):

    SM = 0.05 #PLACEHOLDER

    if tail_config == 'fuselage-mounted':
        Vh_V_2 = 0.85
    elif tail_config == 'fin-mounted':
        Vh_V_2 = 0.95
    elif tail_config == 't-tail';
        Vh_V_2 = 1
    else:
        print('Check the values of tail configuration')


    # Sh_S = 1/(CL_alpha_h/CL_alpha_A_h*(1-de_da)*l_h/c_bar*Vh_V_2)*x_cg_bar-(x_ac_bar-SM)/(CL_alpha_h/CL_alpha_A_h*(1-de_da)*l_h/c_bar*Vh_V_2)
    # S_h = S*Sh_S

    # dL_A_h = CL_alpha_A_h *da/2*rho*V**2*S
    # dL_h = CL_alpha_A_h * (da-de)/2*rho*V_h**2*S_h
    # C_m = C_m0 + (CL_h + CL_A_h*S_h*V_h/(S*V))*(x_cg-x_n)

    m_tv = 2*dz_h/b
    r = 2*l_h/b

    K_ev = (0.1124+0.1265*sweep_ang_14_c_rad+0.1766*sweep_ang_14_c_rad**2)/r**2+0.1024/r+2
    K_ev0 = 0.1124/r**2+0.1024/r+2

    de_da = K_ev/K_ev0 * (r/(r**2+m_tv**2)*0.4876/np.sqrt(r**2+0.6319+m_tv**2)+(1+(r**2/(r**2+0.7915+5.0734*m_tv**2))**0.3113)*(1-np.sqrt(m_tv**2/(1+m_tv**2))))*CL_alpha_w/(np.pi*A)



    V_h = np.sqrt(Vh_V_2)*V
    a = np.sqrt(R*gamma*T)
    M = V_h/a
    beta = np.sqrt(1-M**2)
    CL_alpha_h = 2*np.pi*A_h/(2+np.sqrt(4+(A_h*beta/eta)**2*(1+np.tan(sweep_ang_015_Ch_rad)**2/beta**2)))
    CL_alpha_A_h = CL_alpha_w * (1+2.15*b_f/b)*S_net/S+np.pi/2*b_f**2/S

    x_ac_fc1 = -1.8/CL_alpha_A_h*b_f*h_f*l_fn/(S*c_bar)
    c_g = S/b
    lambd = c_tip/c_root
    x_ac_fc2 =  0.273/(1+lambd) * b_f*c_g*(b-b_f)/(c_bar**2*(b+2.15*b_f))*np.tan(sweep_ang_14_c_rad)

    x_ac_w = 0.3  #PLACEHOLDER, the value shall be taken from graph E-10, lecture 7
    x_ac_bar = x_ac_w + x_ac_fc1 + x_ac_fc2

    Sh_S = np.arange(0, 1.1, 0.1)
    x_np_bar = x_ac_bar +CL_alpha_h/CL_alpha_A_h * (1-de_da) * Sh_S*l_h/c_bar*Vh_V_2
    x_cg_bar = x_ac_bar + CL_alpha_h/CL_alpha_A_h*(1-de_da)*Sh_S*l_h/c_bar*Vh_V_2-SM

    return Sh_S, x_np_bar, x_cg_bar