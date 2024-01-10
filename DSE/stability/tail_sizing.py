import numpy as np
from DSE import const
# from DSE.stability.EoM import distance_stability
from DSE.Locations import locations
from DSE.stability.Loaddiag import load_diagram_plot

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


# ATTENTION!!!!!  CL_alpha are represented in 1/radians
# def tail_sizing(S, b, c_bar, b_f, h_f, l_fn, sweep_ang_14_c_rad, CL_max, CL_alpha_w, S_net, eta, A_h, l_h,dz_h, sweep_ang_12_Ch_rad, c_root, c_tip, tail_config, V, R, gamma, T, rho):

def horizontal_tail_sizing( eta, V, R, gamma, T, rho, S=3.5, b=6, c_bar=0.619, Cl_alpha_h = 2*np.pi,  CL_max=1.5, l_f=2, CL_0=0.3, sweep_ang_rad = 0, Cm_0_airfoil = 0.012, b_f=0.8, hf_max=0.6, l_fn=1.1, sweep_ang_14_c_rad=0, CL_alpha_w=5.787, sweep_ang_12_Ch_rad=0, c_root=0.6, c_tip=0.5):
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

    K_ev = (0.1124+0.1265*sweep_ang_14_c_rad+0.1766*sweep_ang_14_c_rad**2)/r**2+0.1024/r+2
    K_ev0 = 0.1124/r**2+0.1024/r+2

    de_da = K_ev/K_ev0 * (r/(r**2+m_tv**2)*0.4876/np.sqrt(r**2+0.6319+m_tv**2)+(1+(r**2/(r**2+0.7915+5.0734*m_tv**2))**0.3113)*(1-np.sqrt(m_tv**2/(1+m_tv**2))))*CL_alpha_w/(np.pi*A)
    V_h = np.sqrt(Vh_V_2)*V
    a = np.sqrt(R*gamma*T)
    M = V_h/a
    beta = np.sqrt(1-M**2)
    i = 'True'
    if i=='True':
        # Cl_alpha_h is taken as 2pi, but it should be an input from aero
        CL_alpha_h = Cl_alpha_h*A_h/(2+np.sqrt(4+(A_h*beta/eta)**2*(1+np.tan(sweep_ang_12_Ch_rad)**2/beta**2)))
        CL_alpha_A_h = CL_alpha_w * (1+2.15*b_f/b)*S_net/S+np.pi/2*b_f**2/S
        i = 'False'
        x_ac_fc1 = -1.8/CL_alpha_A_h*b_f*hf_max*l_fn/(S*c_bar)
        c_g = S/b
        lambd = c_tip/c_root
        x_ac_fc2 =  0.273/(1+lambd) * b_f*c_g*(b-b_f)/(c_bar**2*(b+2.15*b_f))*np.tan(sweep_ang_14_c_rad)

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

        sr =[]
        for element in range(len(x_cg_bar)):

            x_cg_bar_sel =  x_cg_bar[element] #double check this codeline
            x_np_bar_sel = x_np_bar[x_cg_bar==x_cg_bar_sel] #double check this codeline

            delta_f_Cm_ac = 0 #INCORRECT VALUE, it is assumed 0 as there may not be any flaps
            delta_nac_Cm_ac = 0 #INCORRECT VALUE, it is assumed 0 as there may

            delta_fus_Cm_ac = -1.8*(1-2.5*b_f/l_f)*np.pi*b_f*hf_max*l_f/(4*S*c_bar)*CL_0/CL_alpha_A_h
            Cm_ac_w = Cm_0_airfoil * (A*np.cos(sweep_ang_rad)**2/(A+2*np.cos(sweep_ang_rad)))


            Cm_ac = Cm_ac_w + delta_f_Cm_ac + delta_fus_Cm_ac + delta_nac_Cm_ac


            x_cg_bar_c = x_ac_bar_x - Cm_ac/CL_A_h + CL_h/CL_A_h *Sh_S*l_h/c_bar*Vh_V_2


            delta_x_cg_bar =np.abs(x_cg_bar[element] - x_cg_bar_c[element])
            surface_ratio = ((delta_x_cg_bar + x_np_bar_sel-x_cg_bar_sel- Cm_ac/CL_max) / ((CL_alpha_h/CL_alpha_A_h*(1-de_da)-CL_h/CL_max)*Vh_V_2*l_h/c_bar))
            # print(x_np_bar_sel-x_cg_bar_sel, delta_x_cg_bar)

            sr.append(surface_ratio)
    return  sr


t = horizontal_tail_sizing(eta, V, R, gamma, T, rho)
print('test,', t)


# def vertical_tail_size_1(tail_config='t-tail'):
#     # 'Sv_bv is still the coupled ratio of vertical tail span and surface area of the both sections. Sv1_bv1 is the coupled ratio of the vertical tail and span of one of the the vertical tail sections.
def vertical_tail_size():
    # base imports.
    deg2rad = np.pi / 180
    b = 6
    S = 3.5
    l_fus = 6
    eta = 0.9
    b_max = 0.7
    Cl_alpha = 1.0  # the cl_alpha of the vertical tail at transitional velocity
    CL = 1.3  # at transitional velocity
    Xcg = 3.5147906877402817
    AR_w = b ** 2 / S

    # initial starting values
    lv = 3
    tail_volume = 0.055
    C_eta_beta = 0.058
    taper_v = 1 / 0.40
    step_size = 50

    Bp =    3     # number of blades of the pusher propellor
    lp =  locations()[4]      # the height difference between propellor and c.g.
    Dp =    2*0.5    # Diameter of the propellor
    lf =   2     # length of the fuselage
    hf_max = 0.9   # the maximum height of the fuselage
    # hf1 =       # height of the fuselage at 25% length
    # hf2 =       # height of the fuselage at 75% length
    # bf1 =       # width of the fuselage at 25% length
    # bf2 =       # width of the fuselage at 75% length
    # lcg =       # distance from nose to cg
    # Sfy =       # side surface area of the fuselage
    # intergration space
    AR_v = np.arange(0.5, 2, 0.1)
    beta = 30
    sweep_v = np.arange(0, 45, 1)
    aa = np.linspace(0, l_fus, num=step_size)

    # Empty list set
    Surface = []
    span = []
    Moment_arm = []

    Sv = tail_volume * S * b / lv

def vertical_tail_surface(V_tail_volume=0.055):
    print('first run "vertical_tail_surface_1".')
    lv = locations()[8]
    # S =  #aero import
    # b =  #aero import
    # c_h_tip = # import horizontal tail surface.
    AR_v = 1.25 #current base from lit review (Sadraey, M., Aircraft Design: A Systems Engineering Approach, 2012) #import
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

                C_eta_beta_w = CL ** 2 / (4 * np.pi * AR_w)  # + CL_h**2 / (4*np.pi*AR_h)* (Sh*bh) / (S*b)
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

def rudder_surface_sizing(S_v, l_v, S, b, V_cross, V_trans, S_fus_side, X_AreaCent_fus, rho, C_L_v_alpha = 0.1, C_d_y = 0.8):
    # Method from: O., A.-S., R., A., and H. S., H., “An Educational Rudder Sizing Algorithm for Utilization in Aircraft Design Software,” Tech. Rep. 10, 2018

    # Typical Cn_Beta values 0.04-0.11/rad for subsonic single engine aircraft (SEAD lecture 9)
    Rat_br_bv = np.linspace(0.7, 1.0, 500)  # Ratio of vertical tail fitted with rudder !!!!!Check if this is not in conflict with max deflected elevator!!!!!!

    # Determining crosswind force first
    C_g_fwd, C_g_aft = min(load_diagram_plot()[1]), max(load_diagram_plot()[1])
    V_total = np.sqrt(V_trans**2 + V_cross**2)
    Fus_dist_fwd, Fus_dist_aft = X_AreaCent_fus - C_g_fwd, X_AreaCent_fus - C_g_aft
    F_crosswind = 0.5 * rho * V_cross**2 * S_fus_side * C_d_y
    Sideslip_ang = np.arctan(V_cross/V_trans)
    # C_n_Beta =
    # C_y_Beta =
    Rat_cr_cv = np.linspace(0.15, 0.40, 500)  # Ratio of rudder (mean aerodynamic) chord to total elevator (mean aerodynamic) chord
    Tau_rudder = 1.129 * Rat_cr_cv ** 0.4044 - 0.1772  # O., A.-S., R., A., and H. S., H., “An Educational Rudder Sizing Algorithm for Utilization in Aircraft Design Software,” Tech. Rep. 10, 2018
    # If Tau_rudder is larger than 1, redesign required
    C_n_delta_r = -1 * C_L_v_alpha * ((S_v * l_v) / (S * b)) * Tau_rudder * Rat_br_bv  # Same source
    C_y_delta_r = C_L_v_alpha * eta_v * Tau_rudder * Rat_cr_cv * (S_v/S)


    return

def aileron_surface_sizing():



    return