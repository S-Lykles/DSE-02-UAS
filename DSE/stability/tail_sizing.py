import numpy as np
from DSE import const
from DSE.Locations import locations
from DSE.stability.Loaddiag import load_diagram_plot
from DSE.structures.center_of_gravity import class_two_cg_estimation

# default values
T = 288.15
#gamma = 1.4
#R = 287
eta = 0.95
# A_h = b_h/S_h

def horizontal_tail_sizing(S=3.5, b=6, c_bar=0.619, b_f, hf_max, l_fn, sweep_ang_14_c_rad, CL_alpha_w, S_net, eta, A_h, sweep_ang_015_Ch_rad, c_root, c_tip, V, R, gamma, T, rho, tail_config= 't-tail'):
    l_h  = locations()[3]
    dz_h = locations()[7]
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

    x_ac_fc1 = -1.8/CL_alpha_A_h*b_f*hf_max*l_fn/(S*c_bar)
    c_g = S/b
    lambd = c_tip/c_root
    x_ac_fc2 =  0.273/(1+lambd) * b_f*c_g*(b-b_f)/(c_bar**2*(b+2.15*b_f))*np.tan(sweep_ang_14_c_rad)

    x_ac_w = 0.3  #PLACEHOLDER, the value shall be taken from graph E-10, lecture 7
    x_ac_bar = x_ac_w + x_ac_fc1 + x_ac_fc2

    Sh_S = np.arange(0, 1.1, 0.1)
    x_np_bar = x_ac_bar +CL_alpha_h/CL_alpha_A_h * (1-de_da) * Sh_S*l_h/c_bar*Vh_V_2
    x_cg_bar = x_ac_bar + CL_alpha_h/CL_alpha_A_h*(1-de_da)*Sh_S*l_h/c_bar*Vh_V_2-SM

    return Sh_S, x_np_bar, x_cg_bar,Vh_V_2

def vertical_tail_size():
    'Assumptions made in this code: (Cl_alpha=1.0, CL_w=1.3,  Taper ratio=1/0.4, The outputs where generated on 10-1-2024 becarefull. '

    # base imports.
    deg2rad = const.deg2rad
    b =  #aero import
    S =  #aero import
    l_fus =  #struct import
    eta = 0.95 # vertical tail efficiency
    b_max = # struct import maximum width
    Cl_alpha = 1.0  # the cl_alpha of the vertical tail at cruise speed
    CL_w =  1.3 # CL of the wing tail at cruise speed
    Xcg  = class_two_cg_estimation()[1][0] #locations import
    AR_w = b**2 / S
    Vtrans = # The transition velocity
    v_v =  # the maximum perpendicual gust velocity

    # initial starting values these are assumtions
    lv = 3
    tail_volume = 0.055    # This tail volume is based on literature study for single small propellor aircraft.
    C_eta_beta = 0.058     # This moment coefficient is based on literature study for single small propellor aircraft, this is coupled to the tail volume.
    taper_v = 1 / 0.40     ### This is a assumption

    # intergration space
    AR_v = np.arange(0.5, 2, 0.1)
    beta = 30  # this should still be changed to   np.atan(v_v / Vtrans) * const.rad2deg
    sweep_v = np.arange(0, 45, 1)

    # Empty list set
    Surface = []
    span = []
    Moment_arm = []

    Sv = tail_volume * S * b / lv

    for p in range(len(AR_v)):
        Surface_k = []
        span_k = []
        moment_arm_k = []
        for j in range(len(sweep_v)):
            for k in range(100):
                bv = np.sqrt(AR_v[p] * Sv)

                Cv = 2 / (1 + taper_v) * Sv / bv
                Cv_bar = 2 / 3 * Cv * ((1 + taper_v + taper_v ** 2) / (1 + taper_v))

                # Updated values
                X_LEMAC_v = bv / 6 * ((1 + 2 * taper_v) / (1 + taper_v)) * np.tan(sweep_v[j] * deg2rad)
                lv = 6 - Xcg - Cv + X_LEMAC_v + 0.25 * Cv_bar

                # Update C_eta_beta
                sweep_05_cord_v = sweep_v[j]  # for now a constant sweep is assumed
                Cl_v_alpha = (Cl_alpha * AR_w) / (2 + np.sqrt(4 + (AR_v[p] * beta * deg2rad / eta) * (
                            1 + (np.tan(sweep_05_cord_v * deg2rad) / (beta * deg2rad)) ** 2)))

                C_eta_beta_w = CL_w ** 2 / (4 * np.pi * AR_w)  # + CL_h**2 / (4*np.pi*AR_h)* (Sh*bh) / (S*b)
                new = (np.pi * l_fus * b_max ** 2) / 3
                C_eta_beta_fuse = -2 / (S * b) * new

                # Update tail surface
                Sv = (C_eta_beta - C_eta_beta_fuse - C_eta_beta_w) / Cl_v_alpha * (S * b) / lv

                # List of all values
            Surface_k.append(Sv)
            span_k.append(bv)
            moment_arm_k.append(lv)
        Surface.append(Surface_k)
        span.append(span_k)
        Moment_arm.append(moment_arm_k)

    plot = False
    if plot == True:
        fig, (ax, ay, az) = plt.subplots(1, 3)
        cp = ax.contourf(sweep_v, AR_v, Surface)
        fig.colorbar(cp)
        ax.set_title('Vertical tail surface (Sv)')
        ax.set_xlabel('Sweep angle vertical tail')
        ax.set_ylabel('Aspect ratio vertical tail')

        cy = ay.contourf(sweep_v, AR_v, span)
        fig.colorbar(cy)
        ay.set_title('Span of a single vertical tail plaine (bv)')
        ay.set_xlabel('Sweep angle vertical tail')
        ay.set_ylabel('Aspect ratio vertical tail')

        cz = az.contourf(sweep_v, AR_v, Moment_arm)
        fig.colorbar(cz)
        az.set_title('Moment arm of the vertical tail plaine (lv)')
        az.set_xlabel('Sweep angle vertical tail')
        az.set_ylabel('Aspect ratio vertical tail')

    # Values for now just so there is a output!!
    Sv_1 = 0.625
    bv_1 = 0.6
    lv_1 = 2.3
    AR_v_1 = 1.125
    Sweep_v_1 = 30 #deg
    return Sv_1, bv_1, lv_1, AR_v_1, Sweep_v_1,taper_v

def elevator_surface_sizing(c_bar=0.619,Cm_0=-0.111,Cm_alpha=-0.029,alpha=0,alpha_0=0,CL_alpha_h= 0.12,bh_be=1):
    # speed range ( Stall <-> Max + safety margin)

    # c_bar =       # aero import
    # Cm_0 =        # aero import
    # Cm_alpha =    # aero import
    # alpha =       # general list
    # alpha_0 =     # aero import
    # CL_alpha_h =  # aero import
    Sh_S = horizontal_tail_sizing()[0]
    Vh_V_2 = horizontal_tail_sizing()[3]

    delta = 25  # Elevator deflection range ( -25 <-> 25 degrees)
    l_h = locations()[3]

    Cm_delta_el = -1*(Cm_0 + Cm_alpha*(alpha - alpha_0)) / (delta)
    Tau_el = -1*Cm_delta_el / CL_alpha_h * (bh_be) * 1/Sh_S * c_bar/l_h * (1/Vh_V_2)**2

    return Tau_el, Cm_delta_el

def rudder_surface_sizing(C_L_v_alpha = 0.1, S_v = vertical_tail_surface(), l_v, S, b, V_cross, V_trans, S_fus_side, X_AreaCent_fus, rho, C_d_y = 0.8, V_max):
    """Function to determine minimum rudder chord based on desired crosswind to correct for.

    Imports respectively vertical tail lift curve slope, vertical tail surface area, moment arm for vertical tail, wing surface area, wing span, max cross wind speed,
    transition airspeed, fuselage side area, area center of fuselage, air density, fuselage sideways drag coefficient (assumed to be 0.8), max airspeed (necessary for the max tail load)"""

    # Method from: O., A.-S., R., A., and H. S., H., “An Educational Rudder Sizing Algorithm for Utilization in Aircraft Design Software,” Tech. Rep. 10, 2018
    # Still need Cno

    # Typical Cn_Beta values 0.04-0.11/rad for subsonic single engine aircraft (SEAD lecture 9)
    Rat_br_bv = 0.9  # Ratio of vertical tail fitted with rudder !!!!!Check if this is not in conflict with max deflected elevator!!!!!!

    # Determining crosswind force first
    C_g_fwd, C_g_aft = min(load_diagram_plot()[1]), max(load_diagram_plot()[1])
    V_total = np.sqrt(V_trans**2 + V_cross**2)
    Fus_dist_fwd, Fus_dist_aft = X_AreaCent_fus - C_g_fwd, X_AreaCent_fus - C_g_aft
    F_crosswind = 0.5 * rho * V_cross**2 * S_fus_side * C_d_y
    Sideslip_ang = np.arctan(V_cross/V_trans)
    q = 0.5 * rho * V_total**2

    # Determining control derivatives
    C_n_Beta = K_f_2 * C_L_v_alpha * (1 - dsigma_dbeta) * eta_v * ((l_v * S_v)/(b * S))
    C_y_Beta = -1 * K_f_1 * C_L_v_alpha * (1 - dsigma_dbeta) * eta_v * (S_v/S)
    Rat_cr_cv = np.linspace(0.15, 0.40, 500)  # Ratio of rudder (mean aerodynamic) chord to total elevator (mean aerodynamic) chord
    Tau_rudder = 1.129 * Rat_cr_cv ** 0.4044 - 0.1772  # O., A.-S., R., A., and H. S., H., “An Educational Rudder Sizing Algorithm for Utilization in Aircraft Design Software,” Tech. Rep. 10, 2018
    # If Tau_rudder is larger than 1, redesign required
    C_n_delta_r = -1 * C_L_v_alpha * ((S_v * l_v) / (S * b)) * Tau_rudder * Rat_br_bv  # Same source
    C_y_delta_r = C_L_v_alpha * eta_v * Tau_rudder * Rat_cr_cv * (S_v/S)

    # Solving equations 22 and 23 from the method source
    sigma = np.linspace(0, np.pi, 500)
    delta_r = np.linspace(0, np.pi, 500)
    eqn_22 = 0.5 * rho * V_total**2 * S * b * (C_n_0 + C_n_Beta * (Sideslip_ang - sigma) + C_n_delta_r * delta_r) + F_crosswind * d_c * np.cos(sigma)
    eqn_23 = F_crosswind - q * S * (C_y_0 + C_y_Beta * (Sideslip_ang - sigma) + C_y_delta_r * delta_r)

    int_func = F_crosswind * (1 + Fus_dist_aft * np.cos(sigma)) + q * S * (b * C_n_0 - C_y_0 + (b * C_n_Beta - C_y_Beta)*(Sideslip_ang - sigma) + (b * C_n_delta_r - C_y_delta_r) * delta_r)


    C_l_v = C_l_v_0 + C_l_v_Beta * Sideslip_ang + C_l_v_delta_r * delta_r
    L_v = q * S_v * C_l_v
    Rudder_load = l_v * L_v

    return delta_r, Rat_cr_cv, Rudder_load

def aileron_surface_sizing():



    return