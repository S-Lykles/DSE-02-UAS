import numpy as np
from DSE import const
from DSE.stability.EoM import distance_stability

# default values
T = 288.15
#gamma = 1.4
#R = 287
eta = 0.95
# A_h = b_h/S_h

def horizontal_tail_sizing(S=3.5, b=6, c_bar=0.619, b_f, h_f, l_fn, sweep_ang_14_c_rad, CL_alpha_w, S_net, eta, A_h, sweep_ang_015_Ch_rad, c_root, c_tip, V, R, gamma, T, rho, tail_config= 't-tail'):
    l_h  = distance_stability()[3]
    dz_h = distance_stability()[7]
    SM = 0.05 #PLACEHOLDER

    if tail_config == 'fuselage-mounted':
        Vh_V_2 = 0.85
    elif tail_config == 'fin-mounted':
        Vh_V_2 = 0.95
    elif tail_config == 't-tail':
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

    return Sh_S, x_np_bar, x_cg_bar,Vh_V_2

def vertical_tail_size_1(tail_config='t-tail'):
    'Sv_bv is still the coupled ratio of vertical tail span and surface area of the both sections. Sv1_bv1 is the coupled ratio of the vertical tail and span of one of the the vertical tail sections.

    Bp =        # number of blades of the pusher propellor
    lp =  distance_stability()[4]      # the height difference between propellor and c.g.
    Dp =        # Diameter of the propellor
    lf =        # length of the fuselage
    hf_max =     # the maximum height of the fuselage
    hf1 =       # height of the fuselage at 25% length
    hf2 =       # height of the fuselage at 75% length
    bf1 =       # width of the fuselage at 25% length
    bf2 =       # width of the fuselage at 75% length
    lcg =       # distance from nose to cg
    Sfy =       # side surface area of the fuselage

    if tail_config == 'fuselage-mounted':
        C_eta_beta_i =  0.024
    elif tail_config == 'fin-mounted':
        C_eta_beta_i =  0.012
    elif tail_config == 't-tail':
        C_eta_beta_i = -0.017

    k_beta = 0.3 * (lcg/lf) + 0.75* hf_max/lf -0.105
    C_eta_beta_fuse = -1*k_beta * (Sfy*lf)/(S*b)* np.sqrt(hf1/hf2) * (bf2/bf1)**(1/3)
    C_eta_beta_prop = -0.053 *Bp *(lp*Dp**2)/(S*b)
    C_eta_beta = C_eta_beta_fuse +C_eta_beta_prop + C_eta_beta_i
    print ('the value of C_eta_beta is:', C_eta_beta)
    print ('In order to determine the vertical tail surface, the tail volume of the tail has to be found from the graph using C_eta_beta')
    return C_eta_beta

def vertical_tail_surface(V_tail_volume=0.055):

    lv = distance_stability()[8]
    S =  #aero import
    b =  #aero import

    Sv = V_tail_volume *S*b/lv
    return Sv


def control_surface_sizing(c_bar=0.619,Cm_0=-0.111,Cm_alpha=-0.029,alpha=0,alpha_0=0,CL_alpha_h= 0.12,bh_be=1):
    # speed range ( Stall <-> Max + safety margin)

    # c_bar =
    # Cm_0 =
    # Cm_alpha =
    # alpha =
    # alpha_0 =
    # CL_alpha_h =
    Sh_S = horizontal_tail_sizing()[0]
    Vh_V_2 = horizontal_tail_sizing()[3]

    delta = 25  # Elevator deflection range ( -25 <-> 25 degrees)
    l_h = distance_stability()[3]

    Cm_delta_el = -1*(Cm_0 + Cm_alpha*(alpha - alpha_0)) / (delta)
    Tau_el = -1*Cm_delta_el / CL_alpha_h * (bh_be) * 1/Sh_S * c_bar/l_h * (1/Vh_V_2)**2

    return Tau_el, Cm_delta_el

