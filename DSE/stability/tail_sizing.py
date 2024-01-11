import numpy as np
from DSE import const
from DSE.Locations import locations
from DSE.stability.Loaddiag import load_diagram_plot
from DSE.structures.center_of_gravity import class_two_cg_estimation
from DSE.aero import aero_constants
import matplotlib.pyplot as plt

# default values
T = const.T0
gamma = const.gamma
R = const.R
rho =const.rho0
V = const.V_min # change this value to cruise speed
eta = 0.95 #airfoil efficiency
z_t = 0.5 #PLACEHOLDER, please change it ASAP
z_w = 0 #PLACEHOLDER, please change it ASAP
dz_h = z_t-z_w
S = aero_constants.S
b = aero_constants.b
c_bar = aero_constants.c_bar
# Cl_alpha_h = aero_constants.Cl_alpha_h
CL_max = aero_constants.CL_max
# l_f = aero_constants.l_f
CL_0 = aero_constants.CL_0
# sweep_ang_rad = aero_constants.sweep_ang_rad
# Cm_0_airfoil = aero_constants.Cm_0_airfoil
# b_f = aero_constants.b_f
# hf_max = aero_constants.hf_max
# l_fn = aero_constants.l_fn
# sweep_ang_14_c_rad = aero_constants.sweep_ang_14_c_rad
# CL_alpha_w = aero_constants.CL_alpha_w
# sweep_ang_12_Ch_rad = aero_constants.sweep_ang_12_Ch_rad
c_root = aero_constants.c_root
c_tip = aero_constants.c_tip
c_tip = aero_constants.c_tip
c_tip = aero_constants.c_tip
c_tip = aero_constants.c_tip


# ATTENTION!!!!!  CL_alpha are represented in 1/radians
# def tail_sizing(S, b, c_bar, b_f, h_f, l_fn, sweep_ang_25_c_rad, CL_max, CL_alpha_w, S_net, eta, A_h, l_h,dz_h, sweep_ang_50_c_rad, c_root, c_tip, tail_config, V, R, gamma, T, rho):
# def horizontal_tail_sizing(eta, V, R, gamma, T, rho, S, b, c_bar, Cl_alpha_h = 2*np.pi,  CL_max=1.5, l_f=2, CL_0=0.3, sweep_ang_rad = 0, Cm_0_airfoil = 0.012, b_f=0.8, hf_max=0.6, l_fn=1.1, sweep_ang_25_c_rad=0, CL_alpha_w=5.787, sweep_ang_50_c_rad=0, c_root=0.6, c_tip=0.5):

def horizontal_tail_sizing(eta = 0.95, V = const.v_cruise, R = const.R, gamma = const.gamma, T = const.T0, rho = const.rho0, S  = aero_constants.S, b = aero_constants.b, c_bar = aero_constants.c_bar, Cl_alpha_h = aero_constants.Cl_alpha_h,  CL_max = aero_constants.S, l_f = 2, CL_0 = aero_constants.CL_0, sweep_ang_rad = aero_constants.sweep_ang_rad, Cm_0_airfoil = aero_constants.Cm_0_airfoil, b_f = 0.8, hf_max = 0.8, l_fn = 0.7, CL_alpha_w = aero_constants.CL_alpha_wing, sweep_ang_25_c_rad = aero_constants.sweep_ang_25_c_rad, sweep_ang_50_c_rad = aero_constants.sweep_ang_50_c_rad, c_root = aero_constants.c_root, c_tip = aero_constants.c_tip):
    l_h  = locations()[3]
    dz_h = locations()[7]
    SM = 0.05 #PLACEHOLDER
    A = b**2/S
    S_net = S - b_f*c_root
    A_h = 2/3*A # Initial guess value for the aspect ratio of the wing is

    # if tail_config == 'fuselage-mounted':
    #     Vh_V_2 = 0.85
    # elif tail_config == 'fin-mounted':
    #     Vh_V_2 = 0.95
    # elif tail_config == 't-tail':
    #     Vh_V_2 = 1
    # else:
    #     print('Check the values of tail configuration')

    Vh_V_2 = 1

    m_tv = 2*dz_h/b
    r = 2*l_h/b

    K_ev = (0.1124+0.1265*sweep_ang_25_c_rad+0.1766*sweep_ang_25_c_rad**2)/r**2+0.1024/r+2
    K_ev0 = 0.1124/r**2+0.1024/r+2

    de_da = K_ev/K_ev0 * (r/(r**2+m_tv**2)*0.4876/np.sqrt(r**2+0.6319+m_tv**2)+(1+(r**2/(r**2+0.7915+5.0734*m_tv**2))**0.3113)*(1-np.sqrt(m_tv**2/(1+m_tv**2))))*CL_alpha_w/(np.pi*A)
    V_h = np.sqrt(Vh_V_2)*V
    a = np.sqrt(R*gamma*T)
    M = V_h/a
    beta = np.sqrt(1-M**2)
    i = 'True'
    if i=='True':
        # Cl_alpha_h is taken as 2pi, but it should be an input from aero
        CL_alpha_h = Cl_alpha_h*A_h/(2+np.sqrt(4+(A_h*beta/eta)**2*(1+np.tan(sweep_ang_50_c_rad)**2/beta**2)))
        CL_alpha_A_h = CL_alpha_w * (1+2.15*b_f/b)*S_net/S+np.pi/2*b_f**2/S
        i = 'False'
        x_ac_fc1 = -1.8/CL_alpha_A_h*b_f*hf_max*l_fn/(S*c_bar)
        c_g = S/b
        lambd = c_tip/c_root
        x_ac_fc2 =  0.273/(1+lambd) * b_f*c_g*(b-b_f)/(c_bar**2*(b+2.15*b_f))*np.tan(sweep_ang_25_c_rad)

        x_ac_w = 0.3  #PLACEHOLDER, the value shall be taken from graph E-10, lecture 7 (can be set as input from graph according to wing design)
        x_ac_bar = x_ac_w + x_ac_fc1 + x_ac_fc2

        Sh_S = np.arange(0, 1.1, 0.1)
        x_np_bar = x_ac_bar +CL_alpha_h/CL_alpha_A_h * (1-de_da) * Sh_S*l_h/c_bar*Vh_V_2
        x_cg_bar = x_ac_bar + CL_alpha_h/CL_alpha_A_h*(1-de_da)*Sh_S*l_h/c_bar*Vh_V_2-SM


        # controllability

        # if tail_ability == 'full moving tail':
        #     CL_h = -1
        # elif tail_ability == 'adjustable tail':
        #     CL_h = -0.8
        # elif tail_ability == 'fixed tail':
        #     CL_h = -0.35*A_h**(1/3)
        # else:
        #     print('Check the values of tail ability')

        CL_h = -0.35 * A_h ** (1 / 3)
        x_ac_bar_x = x_ac_bar #PLACEHOLDER

        CL_A_h = 1.1 #INCORRECT VALUE, it is an input, currently just an assumtion

        delta_f_Cm_ac = 0  # INCORRECT VALUE, it is assumed 0 as there may not be any flaps
        delta_nac_Cm_ac = 0  # INCORRECT VALUE, it is assumed 0 as there may

        delta_fus_Cm_ac = -1.8 * (1 - 2.5 * b_f / l_f) * np.pi * b_f * hf_max * l_f / (
                    4 * S * c_bar) * CL_0 / CL_alpha_A_h
        Cm_ac_w = Cm_0_airfoil * (A * np.cos(sweep_ang_rad) ** 2 / (A + 2 * np.cos(sweep_ang_rad)))

        Cm_ac = Cm_ac_w + delta_f_Cm_ac + delta_fus_Cm_ac + delta_nac_Cm_ac

        x_cg_bar_c = x_ac_bar_x - Cm_ac / CL_A_h + CL_h / CL_A_h * Sh_S * l_h / c_bar * Vh_V_2
        sr =[]
        for element in range(len(x_cg_bar)):

            x_cg_bar_sel =  x_cg_bar[element] #double check this codeline
            x_np_bar_sel = x_np_bar[x_cg_bar==x_cg_bar_sel] #double check this codeline

            delta_x_cg_bar =np.abs(x_cg_bar[element] - x_cg_bar_c[element])
            surface_ratio = ((delta_x_cg_bar + x_np_bar_sel-x_cg_bar_sel- Cm_ac/CL_max) / ((CL_alpha_h/CL_alpha_A_h*(1-de_da)-CL_h/CL_max)*Vh_V_2*l_h/c_bar))
            # print(x_np_bar_sel-x_cg_bar_sel, delta_x_cg_bar)

            sr.append(surface_ratio)
    return  sr


t = horizontal_tail_sizing()#eta, V, R, gamma, T, rho)
print('test,', t)


def vertical_tail_size_1(l_fus=2,eta=0.95,b_max=0.7,b=aero_constants.b,S=aero_constants.S,CL_w=aero_constants.CL_max,Cl_alpha=aero_constants.CL_alpha_wing,Xcg=class_two_cg_estimation(True, False, False,  False)[1][0],deg2rad=const.deg2rad):
    """Sv_bv is still the coupled ratio of vertical tail span and surface area of the both sections.
     Sv1_bv1 is the coupled ratio of the vertical tail and span of one of the the vertical tail sections.
     Assumptions made during these calculations:
     Currently there is no interaction between the wing,body,horizontal tail and vertical tail (d_sigma / d_beta = 0)
     The induced velocity interaction between the tail-less aircraft and vertical tail is assumed to be 1=( V_hv / V)**2
     The location of the CG and fuselage length where assumed on 10/01/24 and can therefore differ from the current design.
     The tail volume and Yawing moment coefficient used in this calculation where based on a literature study on small single propellor aircraft.
     A small taper ratio was used during the sizing this was based on a literature study regression.
     The sweep angle was keep constant allong the cord of the tail for now, there is a option to change this within the code. (sweep_05_cord_v)"""

    # base imports.
    AR_w = b ** 2 / S

    # initial starting values
    lv = 3
    tail_volume = 0.055   # Literature study regression coupled to C_eta_beta
    C_eta_beta = 0.058    # Literature study regression coupled to tail volume
    taper_v = 1 / 0.40    # Literature study regression on small propellor aircraft.


    # intergration space
    AR_v = np.arange(0.5, 2, 0.1)
    beta = 30
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
        fig.colorbar(cp)  # Add a colorbar to a plot
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

    #return

def vertical_tail_size():
    """Preliminary vertical tail size, it can change ask Bas before using. Last update 10/01/24."""
    Sv = 1.250/2
    bv = 1.2
    l_v = 2.3
    AR_v = 1.125
    Sweep_angle_v = 30 #deg
    taper_v = 1/0.4
    print('These values can be subject to change and where last updated on 10/01/24. Before using ask Bas, if there are any changes due to changing design parameters.')
    return Sv,bv,l_v,AR_v,Sweep_angle_v,taper_v

def elevator_surface_sizing(l_h=locations()[3],c_bar=aero_constants.c_bar,Cm_0=aero_constants.Cm_0_airfoil,Cm_alpha=aero_constants.Cm_alpha,alpha=0,alpha_0=aero_constants.alpha_0,CL_alpha_h= 0.12,bh_be=1):
    # speed range ( Stall <-> Max + safety margin)
    Sh_S = horizontal_tail_sizing()[0]
    Vh_V_2 = horizontal_tail_sizing()[3]

    delta = 25  # Elevator deflection range ( -25 <-> 25 degrees)

    Cm_delta_el = -1*(Cm_0 + Cm_alpha*(alpha - alpha_0)) / delta
    Tau_el = -1*Cm_delta_el / CL_alpha_h * (bh_be) * 1/Sh_S * c_bar/l_h * (1/Vh_V_2)**2

    return Tau_el, Cm_delta_el

#def rudder_surface_sizing(S_v, l_v, S, b, V_cross, V_trans, S_fus_side, X_AreaCent_fus, rho, C_L_v_alpha = 0.1, C_d_y = 0.8):
def rudder_surface_sizing( l_v,  V_cross, V_trans, S_fus_side, X_AreaCent, rho, V_max, C_L_v_alpha = 0.1, S_v = vertical_tail_size()[0], S = aero_constants.S, b = aero_constants.b, C_d_y = 0.8, dsigma_dbeta = 0.0, eta_v = 0.95, C_n_0 = 0.0, C_y_0 = 0.0):
    """Function to determine minimum rudder chord based on desired crosswind to correct for.

    !!!Currently the vertical tail span that is fitted with a rudder is assumed to be 90% of the total span, when an elevator chord is determined, it must be made sure that elevator and rudder do not collide at maximum deflection!!!"""

    # Method from: O., A.-S., R., A., and H. S., H., “An Educational Rudder Sizing Algorithm for Utilization in Aircraft Design Software,” Tech. Rep. 10, 2018

    # Typical Cn_Beta values 0.04-0.11/rad for subsonic single engine aircraft (SEAD lecture 9)
    Rat_br_bv = 0.9  # Ratio of vertical tail fitted with rudder !!!!!Check if this is not in conflict with max deflected elevator!!!!!!
    S_v_total = 2 * S_v # Correction for having 2 vertical stabilisers

    # Determining crosswind force first
    C_g_fwd, C_g_aft = min(load_diagram_plot()[1]), max(load_diagram_plot()[1])
    V_total = np.sqrt(V_trans**2 + V_cross**2)
    Fus_dist_fwd, Fus_dist_aft = X_AreaCent - C_g_fwd, X_AreaCent - C_g_aft
    F_crosswind = 0.5 * rho * V_cross**2 * S_fus_side * C_d_y
    beta = np.arctan(V_cross/V_trans)
    q = 0.5 * rho * V_total**2

    # Determining control derivatives
    C_n_Beta = K_f_2 * C_L_v_alpha * (1 - dsigma_dbeta) * eta_v * ((l_v * S_v_total)/(b * S))
    C_y_Beta = -1 * K_f_1 * C_L_v_alpha * (1 - dsigma_dbeta) * eta_v * (S_v_total/S)
    Rat_cr_cv = np.linspace(0.15, 0.40, 500)  # Ratio of rudder (mean aerodynamic) chord to total elevator (mean aerodynamic) chord
    Tau_rudder = 1.129 * Rat_cr_cv ** 0.4044 - 0.1772  # O., A.-S., R., A., and H. S., H., “An Educational Rudder Sizing Algorithm for Utilization in Aircraft Design Software,” Tech. Rep. 10, 2018
    # If Tau_rudder is larger than 1, redesign required
    C_n_delta_r = -1 * C_L_v_alpha * ((S_v * l_v) / (S * b)) * Tau_rudder * Rat_br_bv  # Same source
    C_y_delta_r = C_L_v_alpha * eta_v * Tau_rudder * Rat_cr_cv * (S_v/S)

    # Solving equations 22 and 23 from the method source
    sigma = np.linspace(0, np.pi/2, 500)
    delta_r = np.linspace(0, np.pi/6, 500)
    eqn_22 = q * S * b * (C_n_0 + C_n_Beta * (beta - sigma) + C_n_delta_r * delta_r) + F_crosswind * d_c * np.cos(sigma)
    eqn_23 = F_crosswind - q * S * (C_y_0 + C_y_Beta * (beta - sigma) + C_y_delta_r * delta_r)

    eqn_22_abs = abs(eqn_22)
    eqn_23_abs = abs(eqn_23)


    # int_func = F_crosswind * (1 + Fus_dist_aft * np.cos(sigma)) + q * S * (b * C_n_0 - C_y_0 + (b * C_n_Beta - C_y_Beta)*(beta - sigma) + (b * C_n_delta_r - C_y_delta_r) * delta_r)


    # C_l_v = C_l_v_0 + C_l_v_Beta * beta + C_l_v_delta_r * delta_r
    # L_v = q * S_v * C_l_v
    # Rudder_load = l_v * L_v

    return delta_r_final, Rat_cr_cv_final # , Rudder_load

def aileron_surface_sizing():



    return