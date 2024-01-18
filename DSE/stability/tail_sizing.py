import numpy as np
from DSE import const
from DSE.Locations import locations
from DSE.stability.Loaddiag import load_diagram_plot
from DSE.structures.center_of_gravity import class_two_cg_estimation
from DSE.aero import aero_constants
from loading_diagram import loading_diagram_extremes
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

# ATTENTION!!!!!  CL_alpha are represented in 1/radians
# def tail_sizing(S, b, c_bar, b_f, h_f, l_fn, sweep_ang_25_c_rad, CL_max, CL_alpha_w, S_net, eta, A_h, l_h,dz_h, sweep_ang_50_c_rad, c_root, c_tip, tail_config, V, R, gamma, T, rho):
# def horizontal_tail_sizing(eta, V, R, gamma, T, rho, S, b, c_bar, Cl_alpha_h = 2*np.pi,  CL_max=1.5, l_f=2, CL_0=0.3, sweep_ang_rad = 0, Cm_0_airfoil = 0.012, b_f=0.8, hf_max=0.6, l_fn=1.1, sweep_ang_25_c_rad=0, CL_alpha_w=5.787, sweep_ang_50_c_rad=0, c_root=0.6, c_tip=0.5):

def find_intersection(x1, y1, x2, y2):
    # Calculate slopes
    m1 = (y1[1] - y1[0]) / (x1[1] - x1[0])
    m2 = (y2[1] - y2[0]) / (x2[1] - x2[0])

    # Calculate y-intercepts
    b1 = y1[0] - m1 * x1[0]
    b2 = y2[0] - m2 * x2[0]

    # Check if lines are parallel
    if m1 == m2:
        print("Lines are parallel and don't intersect.")
        return None

    # Calculate intersection point
    x_intersect = (b2 - b1) / (m1 - m2)
    y_intersect = m1 * x_intersect + b1

    return x_intersect, y_intersect

def horizontal_tail_sizing(eta = 0.95, V = const.v_cruise, R = const.R, gamma = const.gamma, T = const.T0, S  = aero_constants.S, b = aero_constants.b, c_bar = aero_constants.c_bar, Cl_alpha_h = aero_constants.Cl_alpha_h,  CL_max = aero_constants.S, l_f = 2, CL_0 = aero_constants.CL_0, sweep_ang_rad = aero_constants.sweep_ang_rad, Cm_0_airfoil = aero_constants.Cm_0_airfoil, b_f = 0.8, hf_max = 0.8, l_fn = 0.7, CL_alpha_w = aero_constants.CL_alpha_wing, sweep_ang_25_c_rad = aero_constants.sweep_ang_25_c_rad, sweep_ang_50_c_rad = aero_constants.sweep_ang_50_c_rad, c_root = aero_constants.c_root, c_tip = aero_constants.c_tip):
    l_h  = locations()[3]
    dz_h = locations()[7]
    #print('lift rate',Cl_alpha_h)
    SM = 0.05 #PLACEHOLDER
    A = b**2/S
    S_net = S - b_f*c_root
    # A_h = 2/3*A # Initial guess value for the aspect ratio of the wing is
    # A_hh = np.arange(0.1, 10.6, 0.1)
    A_hh = np.arange(6.7, 6.9, 0.1)
    Vh_V_2 = 1 # tail configuration is assumed to be t-tail

    m_tv = 2*dz_h/b
    r = 2*l_h/b
    K_ev = (0.1124+0.1265*sweep_ang_25_c_rad+0.1766*sweep_ang_25_c_rad**2)/r**2+0.1024/r+2
    K_ev0 = 0.1124/r**2+0.1024/r+2

    de_da = K_ev/K_ev0 * (r/(r**2+m_tv**2)*0.4876/np.sqrt(r**2+0.6319+m_tv**2)+(1+(r**2/(r**2+0.7915+5.0734*m_tv**2))**0.3113)*(1-np.sqrt(m_tv**2/(1+m_tv**2))))*CL_alpha_w/(np.pi*A)


    V_h = np.sqrt(Vh_V_2)*V
    a = np.sqrt(R*gamma*T)
    M = V_h/a
    beta = np.sqrt(1-M**2)
    sr = []
    d_xcg = []
    dcg=0
    check_y =[]
    for i in range(len(A_hh)):
       # print(i,len(A_hh))
        # print(i,len(A_hh))
        if A_hh[i] == 6.8:
            PloT=True
        else:
            PloT=False
        A_h = A_hh[i]
        # # A_h = 2/3*36/3.5
        # # Cl_alpha_h is taken as 2pi, but it should be an input from aero
        # # CL_alpha_h = 2*np.pi*A_h/(2+np.sqrt(4+(A_h*beta/eta)**2*(1+np.tan(sweep_ang_50_c_rad)**2/beta**2)))
        CL_alpha_h = Cl_alpha_h*A_h/(2+np.sqrt(4+(A_h*beta/eta)**2*(1+np.tan(sweep_ang_50_c_rad)**2/beta**2)))
        # CL_alpha_h =CL_alpha_h/10

        CL_alpha_A_h = CL_alpha_w * (1+2.15*b_f/b)*S_net/S+np.pi/2*b_f**2/S

        x_ac_bar_fc1 = -1.8/CL_alpha_A_h*b_f*hf_max*l_fn/(S*c_bar)
        c_g = S/b # mean geometric chord
        lambd = c_tip/c_root
        x_ac_bar_fc2 =  0.273/(1+lambd) * b_f*c_g*(b-b_f)/(c_bar**2*(b+2.15*b_f))*np.tan(sweep_ang_25_c_rad)
        x_ac_bar_w = 0.25  #PLACEHOLDER,. the value shall be taken from graph E-10, lecture 7 (can be set as input from graph according to wing design)
        # x_ac_bar_w = 0.3+l_fn/c_bar  #PLACEHOLDER, the value shall be taken from graph E-10, lecture 7 (can be set as input from graph according to wing design)
        x_ac_bar = x_ac_bar_w + x_ac_bar_fc1 + x_ac_bar_fc2

        Sh_S = np.arange(0, 0.8, 0.1)
        x_np_bar = x_ac_bar + CL_alpha_h/CL_alpha_A_h * (1-de_da) * Sh_S*l_h/c_bar*Vh_V_2
        x_cg_bar = x_ac_bar + CL_alpha_h/CL_alpha_A_h * (1-de_da) * Sh_S*l_h/c_bar*Vh_V_2 - SM

        CL_h = -0.35 * A_h ** (1 / 3) # It is assumed that tail_ability is 'fixed tail'

        CL_A_h = 0.3 #INCORRECT VALUE, it is an input, currently just an assumption (aircraft less tail lift coefficient)

        delta_f_Cm_ac = 0  # INCORRECT VALUE, it is assumed 0 as there may not be any flaps
        delta_nac_Cm_ac = 0  # INCORRECT VALUE, it is assumed 0 as there is no engine mounted in the wing

        delta_fus_Cm_ac = -1.8 * (1 - 2.5 * b_f / l_f) * np.pi * b_f * hf_max * l_f / (
                    4 * S * c_bar) * CL_0 / CL_alpha_A_h
        Cm_ac_w = Cm_0_airfoil * (A * np.cos(sweep_ang_rad) ** 2 / (A + 2 * np.cos(sweep_ang_rad)))

        Cm_ac = Cm_ac_w + delta_f_Cm_ac + delta_fus_Cm_ac + delta_nac_Cm_ac

        x_cg_bar_c = x_ac_bar - Cm_ac / CL_A_h + CL_h / CL_A_h * Sh_S * l_h / c_bar * Vh_V_2

        x_cg_bar_min, x_cg_bar_max, xlemac_lf, slope_lemac = loading_diagram_extremes()

        y_offset = np.arange(0,10,0.001)
        y_datum_shift=0
        y_diff=5

        x_i, y_i = find_intersection([x_cg_bar[0], x_cg_bar[5]], [Sh_S[0], Sh_S[5]], [x_cg_bar_c[0], x_cg_bar_c[5]], [Sh_S[0], Sh_S[5]])
        for i in range(len(y_offset)):
            y_shift = y_offset[i]

            x_1, y_1 = find_intersection([x_cg_bar[0], x_cg_bar[5]], [Sh_S[0], Sh_S[5]], [x_cg_bar_max[0], x_cg_bar_max[5]], [xlemac_lf[0]-y_shift, xlemac_lf[5]-y_shift])
            x_2, y_2 = find_intersection([x_cg_bar_c[0], x_cg_bar_c[5]], [Sh_S[0], Sh_S[5]], [x_cg_bar_min[0], x_cg_bar_min[5]], [xlemac_lf[0]-y_shift, xlemac_lf[5]-y_shift])

            if (y_i> y_1 )or (y_i> y_2):
                continue
            else:
                y_dist=np.abs(y_1-y_2)
                if y_dist<y_diff:
                    y_datum_shift=y_shift
                    y_diff = y_dist
                    x1, x2, y1, y2 = x_1, x_2, y_1, y_2

        if x1>0 and x2>0:
            check_y.append(np.abs(y1-y2))
        else:
            check_y.append(1000)

        delta_x_cg_bar = x1-x2
        # print('dcg',x1,x2)
        d_xcg.append(delta_x_cg_bar)


        # PloT = False
        # if x1>0 and x2>0:
        #     PloT = True
        # else:
        #     PloT = False


        # surface_ratio = ((delta_x_cg_bar + SM- Cm_ac/CL_max) / ((CL_alpha_h/CL_alpha_A_h*(1-de_da)-CL_h/CL_max)*Vh_V_2*l_h/c_bar))
        # print('plot sr_h:',surface_ratio)
        surface_ratio =np.maximum(y1, y2)
        sr.append(surface_ratio)
        # print('sur.r.:', surface_ratio+y_datum_shift)


        if delta_x_cg_bar> dcg:
            dcg = delta_x_cg_bar
            sf= S*surface_ratio
            arh=A_h

        # print()
        xlemac_lf = [x-y_datum_shift for x in xlemac_lf]
        # xlemac_lf=xlemac_lf-1.57
        PloT = False
        if PloT == True:
            print(surface_ratio*S/2.3,surface_ratio)
            plt.plot(x_cg_bar, Sh_S, label ='stab')
            # plt.plot(x_np_bar, Sh_S, label ='neutral')
            plt.plot(x_cg_bar_c, Sh_S, label ='cont')
            plt.plot(x_cg_bar_max, xlemac_lf, label ='xcg max')
            plt.plot(x_cg_bar_min, xlemac_lf, label ='xcg min')
            plt.plot(x_cg_bar_min, surface_ratio*np.ones(np.shape(x_cg_bar_min)), label ='xcg  optimum')
            plt.plot(x_cg_bar_min, np.minimum(y1,y2)*np.ones(np.shape(x_cg_bar_min)), label ='xcg  minoptimum')
            # plt.xlim(0,1)
            # plt.ylim(0,0.1)
            plt.ylabel('Sh/S')
            plt.xlabel('% Xcg of MAC')
            plt.title('Horizontal tail [AR=6.8]')
            plt.legend()
            plt.show()

    Sh_area = [x*S for x in sr]

    Sh_area0, A_hh0 = np.meshgrid(Sh_area, A_hh)
    span_h = np.sqrt(Sh_area0*A_hh0)
    # chord=Sh_area0/span_h
    N = 400

    # fig, ax = plt.subplots()
    # cp = ax.contourf(Sh_area0, A_hh0, span_h, N)
    # cbar = fig.colorbar(cp)
    # contour_line = ax.contour(Sh_area0, A_hh0, span_h, levels=[2.3], colors='r')
    # ax.legend([contour_line.collections[0]], ['Span = 2.3 m'], loc='upper right')
    # # Find the intersection points
    # contour_data = contour_line.collections[0].get_paths()[0].vertices
    # # intersection_point = contour_data[np.argmin(np.abs(contour_data[:, 1] - 6.8))]
    # #
    # # # Plot a point at the intersection,
    # # ax.plot(intersection_point[0], intersection_point[1], 'o', color='pink', label='Intersection')
    #
    # ax.set_title('Span of the horizontal tail [m]')
    # ax.set_xlabel('Surface area of horizontal tail [m^2]')
    # ax.set_ylabel('Aspect ratio horizontal tail')
    # plt.show()
    # print(np.shape(Sh_area),np.shape(d_xcg))
    #
    # plt.scatter(Sh_area, d_xcg)
    # plt.xlabel('Minimum optimal surface area of the horizontal tail')
    # plt.ylabel('Optimum CG range')
    # plt.show()

    plot2=False
    if plot2 == True:
        fig, (ax, ay, az) = plt.subplots(1, 3)
        cp = ax.contourf(Sh_area0,A_hh0,  span_h, N)
        fig.colorbar(cp)  # Add a colorbar to a plot
        contour_line = ax.contour(Sh_area0,A_hh0,  span_h, levels=[2.3], colors='r', label='Span = 2.3 m')
        ax.set_title('Span of the horizontal tail')
        ax.set_xlabel('Aspect ratio horizontal tail')
        ax.set_ylabel('Surface area of horizontal tail')
        plt.show()
        #
        cy = ay.contourf(A_hh, Sh_area,  d_xcg, N)
        ay.clabel(cy, inline=True, fontsize=10)
        fig.colorbar(cy)
        ay.set_title('Stability margin [m]')
        ax.set_xlabel('Aspect ratio horizontal tail')
        ax.set_ylabel('Surface area of horizontal tail')

    # Sh = S*surface_ratio
    Sh = S*np.min(sr)

    #print('AR', 2.3**2/Sh, Sh/2.3)
    return Sh, x_cg_bar, x_cg_bar_c, surface_ratio,de_da
print("hoi",horizontal_tail_sizing()[4])
test_print = False
if test_print ==True:
    a, b, c, d, e = horizontal_tail_sizing()
    print('test 1', a)
    print('test 2', b)
    print('test 3', c)
    print('test 4', d)

def elevator_surface_sizing(l_h=locations()[3],Sh=horizontal_tail_sizing()[0],CL_0=aero_constants.CL_0,CL_alpha=aero_constants.CL_alpha_wing,c_bar=aero_constants.c_bar,CL_l = 0.91,Cm_0=aero_constants.Cm_0_airfoil,Cm_alpha=aero_constants.Cm_alpha,alpha=0,alpha_0=aero_constants.alpha_0,CL_alpha_h= 0.12,bh_be=1):
    # speed range ( Stall <-> Max + safety margin)
    eta_h = 0.9
    delta = 25*np.pi/180  # Elevator deflection range ( -25 <-> 25 degrees)

    Vh_vol = (Sh*l_h) / (S*c_bar)

    Tau_el = (Cm_0 *CL_alpha + (CL_l - CL_0)*Cm_alpha ) / ( delta*eta_h* ( -CL_alpha*CL_alpha_h*Vh_vol*(1/bh_be) - Cm_alpha* CL_alpha_h * Sh/S*(1/bh_be)))

    S_a_S = np.linspace(0.0, 0.7, 100)
    tau_a = -6.624 * (S_a_S) ** 4 + 12.07 * (S_a_S) ** 3 - 8.292 * (S_a_S) ** 2 + 3.295 * S_a_S + 0.004942 - Tau_el

    print('The Chord-to-Chord and Surface-to-Surface ratio are determined for Tau_el = 0.46, if Tau_el is different check wiht "Jakob or Bas"')
    Ce_Ch = 0.25
    Se_Sh = Ce_Ch * (1/bh_be)
    return Tau_el[0], Ce_Ch, Se_Sh


#def rudder_surface_sizing(S_v, l_v, S, b, V_cross, V_trans, S_fus_side, X_AreaCent_fus, rho, C_L_v_alpha = 0.1, C_d_y = 0.8):
def rudder_surface_sizing( V_cross, V_trans, S_fus_side, X_AreaCent, rho, V_max, C_L_v_alpha = 4.5, S_v =0.1811301138015332, l_v =  2.452651830421447, S = aero_constants.S, b = aero_constants.b, C_d_y = 0.8, dsigma_dbeta = 0.0, eta_v = 0.95, C_n_0 = 0.0, C_y_0 = 0.0, K_f_1 = 0.75, K_f_2 = 1.4):
    """Function to determine minimum rudder chord based on desired crosswind to correct for.

    !!!Currently the vertical tail span that is fitted with a rudder is assumed to be 90% of the total span, when an elevator chord is determined, it must be made sure that elevator and rudder do not collide at maximum deflection!!!"""

    # Method from: O., A.-S., R., A., and H. S., H., “An Educational Rudder Sizing Algorithm for Utilization in Aircraft Design Software,” Tech. Rep. 10, 2018

    # Typical Cn_Beta values 0.04-0.11/rad for subsonic single engine aircraft (SEAD lecture 9)
    Rat_br_bv = 0.9  # Ratio of vertical tail fitted with rudder !!!!!Check if this is not in conflict with max deflected elevator!!!!!!
    S_v_total = 2 * S_v # Correction for having 2 vertical stabilisers
    No_values = 500
    # Tolerance = 5e-3
    # K_f_1 between 0.65 and 0.75 for typical aircraft
    # K_f_2 between 1.3 and 1.4 for typical aircraft

    # Determining crosswind force first
    C_g_fwd, C_g_aft = 1.5, 1.5
    V_total = np.sqrt(V_trans**2 + V_cross**2)
    Fus_dist_fwd, Fus_dist_aft = abs(X_AreaCent - C_g_fwd), abs(X_AreaCent - C_g_aft)
    d_c = max(Fus_dist_fwd, Fus_dist_aft)
    F_crosswind = 0.5 * rho * V_cross**2 * S_fus_side * C_d_y
    beta = np.arctan(V_cross/V_trans)
    q = 0.5 * rho * V_total**2

    # Determining control derivatives
    C_n_Beta = K_f_2 * C_L_v_alpha * (1 - dsigma_dbeta) * eta_v * ((l_v * S_v_total)/(b * S))
    C_y_Beta = -1 * K_f_1 * C_L_v_alpha * (1 - dsigma_dbeta) * eta_v * (S_v_total/S)
    Rat_cr_cv = 0.14  # Ratio of rudder (mean aerodynamic) chord to total elevator (mean aerodynamic) chord
    Tau_rudder = 1.129 * Rat_cr_cv ** 0.4044 - 0.1772  # O., A.-S., R., A., and H. S., H., “An Educational Rudder Sizing Algorithm for Utilization in Aircraft Design Software,” Tech. Rep. 10, 2018
    # If Tau_rudder is larger than 1, redesign required, should not happen with given range for chord ratios
    C_n_delta_r = -1 * C_L_v_alpha * ((S_v * l_v) / (S * b)) * Tau_rudder * Rat_br_bv  # Same source
    C_y_delta_r = C_L_v_alpha * eta_v * Tau_rudder * Rat_cr_cv * (S_v/S)

    # Solving equations 22 and 23 from the method source
    delta_r = np.linspace(0, np.pi/3, No_values)
    delta_r_final = 40 * (np.pi / 180)

    while delta_r_final > (30 * (np.pi / 180)):
        Rat_cr_cv += 0.01
        Tau_rudder = 1.129 * Rat_cr_cv ** 0.4044 - 0.1772  # O., A.-S., R., A., and H. S., H., “An Educational Rudder Sizing Algorithm for Utilization in Aircraft Design Software,” Tech. Rep. 10, 2018
        # If Tau_rudder is larger than 1, redesign required, should not happen with given range for chord ratios
        C_n_delta_r = -1 * C_L_v_alpha * ((S_v * l_v) / (S * b)) * Tau_rudder * Rat_br_bv  # Same source
        C_y_delta_r = C_L_v_alpha * eta_v * Tau_rudder * Rat_cr_cv * (S_v/S)

        sigma = (C_y_0/C_y_Beta) + beta + ((C_y_delta_r/C_y_Beta) * delta_r) - (F_crosswind / (q * S * C_y_Beta))
        eqn_22 = q * S * b * (C_n_0 + C_n_Beta * (beta - sigma) + C_n_delta_r * delta_r) + F_crosswind * d_c * np.cos(sigma)

        k = np.min(np.where(eqn_22 < 0.0))
        delta_r_final = delta_r[k]

    if Rat_cr_cv > 0.4:
        print('Rudder is too big')

    # Rat_cr_cv_final = min(Rat_cr_cv[np.where(deflections < 30)])
    # int_func = F_crosswind * (1 + Fus_dist_aft * np.cos(sigma)) + q * S * (b * C_n_0 - C_y_0 + (b * C_n_Beta - C_y_Beta)*(beta - sigma) + (b * C_n_delta_r - C_y_delta_r) * delta_r)


    # C_l_v = C_l_v_0 + C_l_v_Beta * beta + C_l_v_delta_r * delta_r
    # L_v = q * S_v * C_l_v
    # Rudder_load = l_v * L_v
    #print(Rat_cr_cv)

    return Rat_cr_cv # , Rudder_load

rudder_surface_sizing(10.0, 43.0, 10.0, 1.0, 1.225, 120.0)

def aileron_surface_sizing(V_trans, roll_rate = 0.2618, span_wise_inner_frac = 0.5, aileron_chord_frac = 0.27, deflection_up = 20.0, deflection_down = 20.0, b = aero_constants.b, CL_alpha_w = aero_constants.CL_alpha_wing, CD_0_wing = aero_constants.CD0_wing, taper_w = aero_constants.c_tip / aero_constants.c_root, C_root = aero_constants.c_root, S = aero_constants.S):
    """Aileron sizing for roll rate requirement, assumes a roll rate now, assumes that the aileron starts immediately beyond the propeller projection on the wing, deflects no more than 20 degrees, does not cause adverse yaw and takes up 27% of wing chord. Can all be changed if required"""
    # Let's do this! давай!

    # Typical values:
    # S_a_S (aileron area to wing area) 0.05 - 0.1
    # b_a_b (aileron span to wing span) 0.2 - 0.3
    # C_a_C (aileron chord to wing chord) 0.15 - 0.25
    # delta_a_max (aileron max deflection) 30 deg or 25 deg

    V_man = 1.3 * V_trans # Requirement, I think
    span_wise_outer_frac  = span_wise_inner_frac + 0.01
    span_wise_inner, span_wise_outer = span_wise_inner_frac * (b / 2), span_wise_outer_frac * (b / 2)
    b_a = span_wise_outer - span_wise_inner
    aileron_deflection_up = const.deg2rad(deflection_up) # Assuming Frise aileron (no adverse yaw hooray)
    aileron_deflection_down = const.deg2rad(deflection_down)
    deflection_aileron = 0.5 * (aileron_deflection_up + aileron_deflection_down)

    C_lp = -1 * ((CL_alpha_w + CD_0_wing) * C_root * b / (24 * S)) * (1 + 3 * taper_w)
    Cl_delta_A_req = -1 * ((roll_rate * b * C_lp) / (2 * V_man * deflection_aileron))
    Ca_inner = aileron_chord_frac * C_root * (1 + 2 * ((taper_w - 1) / b) * span_wise_inner)
    Ca_outer = aileron_chord_frac * C_root * (1 + 2 * ((taper_w - 1) / b) * span_wise_outer)
    S_a_S = b_a / S * (Ca_inner + Ca_outer)
    tau_a = -6.624 * (S_a_S) ** 4 + 12.07 * (S_a_S) ** 3 - 8.292 * (S_a_S) ** 2 + 3.295 * S_a_S + 0.004942
    Cl_delta_A = ((CL_alpha_w * tau_a * C_root) / (S * b)) * ((span_wise_outer ** 2 - span_wise_inner ** 2) + (4 / 3) * ((taper_w - 1) / b) * (span_wise_outer ** 3 - span_wise_inner ** 3))

    while Cl_delta_A < Cl_delta_A_req:
        span_wise_outer_frac += 0.001
        span_wise_outer = span_wise_outer_frac * (b / 2)
        b_a = span_wise_outer - span_wise_inner
        Ca_outer = aileron_chord_frac * C_root * (1 + 2 * ((taper_w - 1) / b) * span_wise_outer)
        S_a_S = b_a / S * (Ca_inner + Ca_outer)
        tau_a = -6.624 * (S_a_S) ** 4 + 12.07 * (S_a_S) ** 3 - 8.292 * (S_a_S) ** 2 + 3.295 * S_a_S + 0.004942
        Cl_delta_A = ((CL_alpha_w * tau_a * C_root) / (S * b)) * ((span_wise_outer ** 2 - span_wise_inner ** 2) + (4 / 3) * ((taper_w - 1) / b) * (span_wise_outer ** 3 - span_wise_inner ** 3))

    if span_wise_outer_frac > 1.0:
        print('Aileron is too big, consider decreasing inboard span-wise location!!!')

    span_wise_outer_final = span_wise_outer
    S_a_final = S_a_S * S

    return(span_wise_inner, span_wise_outer_final, S_a_final, Cl_delta_A)


#def stab_con_int_structure():
    """Keep in mind that these values are all subject to change due to interations"""
    print('!!! Carefull when using this Data, it will be subject to change due to shifts in Xcg and rotor placement. !!!')

    # vertical tail dimensions
    S_vert = vertical_tail_size()[0]
    b_vert = vertical_tail_size()[1]
    cord_root_vert = vertical_tail_size()[4]
    cord_tip_vert = vertical_tail_size()[6]*cord_root_vert
    Dist_X_leading_vert = vertical_tail_size()[2] - 0.25* vertical_tail_size()[4] - vertical_tail_size()[7]
    Sweep_angle_vert = vertical_tail_size()[5]

    # rudder dimensions
    cord_root_rudder = rudder_surface_sizing()[0]*cord_root_vert
    bv_rudder = rudder_surface_sizing()[1]*b_vert
    Location_rudder_z = 0       # location with respect to vertical tail, in the horizontal direction.

    # horizontal tail dimensions
    S_hori = horizontal_tail_sizing()[2]
    cord_root_hori = cord_tip_vert
    cord_tip_hori = cord_root_hori
    b_hori = S_hori / cord_root_hori
    Dist_X_leading_hori = Dist_X_leading_vert + b_vert * np.tan(Sweep_angle_vert*const.deg2rad)
    Sweep_angle_hori = 0 # no sweep is assumed at the moment.

    # elevator dimensions
    tau_elevator = elevator_surface_sizing()[0]
    print('using the graph in overleaf the cord-to-cord ratio can be found, and also the surface-to-surface ratio')
    print('this is based on the assumption that the span-to-span ratio is 1')

    # aileron dimensions
    cord_root_aileron = 0
    bv_aileron = 0
    Dist_X_leading_aileron = 0 # b_a - bv_aileron/2

    # forces
    print('The tail forces')
    print('The vertical force is for a single tail surface')
    Fz_tail = S_vert*rho*15**2 # extreem assumption
    Fx_tail = S_hori*rho*V**2

   # return S_vert, b_vert, cord_root_vert, cord_tip_vert, Dist_X_leading_vert, Sweep_angle_vert, cord_root_rudder, bv_rudder, S_hori,cord_root_vert, cord_tip_vert, cord_tip_hori, b_hori, Dist_X_leading_hori, Sweep_angle_hori, tau_elevator, cord_root_aileron, bv_aileron, Dist_X_leading_aileron, Fz_tail, Fx_tail, Fx_tail



def Tail_opt_DO_NOT_RUN(l_fus=2,eta=0.95,b_max=0.7,b=aero_constants.b,S=aero_constants.S,CL=[aero_constants.CL_max,aero_constants.CL_cruise],Cl_alpha_v=aero_constants.Cl_alpha_v,Xcg=1.5,deg2rad=const.deg2rad):
    """Sv_bv is still the coupled ratio of vertical tail span and surface area of the both sections.
     Sv1_bv1 is the coupled ratio of the vertical tail and span of one of the the vertical tail sections.
     Assumptions made during these calculations:
     Currently there is no interaction between the wing,body,horizontal tail and vertical tail (d_sigma / d_beta = 0)
     The induced velocity interaction between the tail-less aircraft and vertical tail is assumed to be 1=( V_hv / V)**2
     The location of the CG and fuselage length where assumed on 10/01/24 and can therefore differ from the current design.
     The tail volume and Yawing moment coefficient used in this calculation where based on a literature study on small single propellor aircraft.
     A small taper ratio was used during the sizing this was based on a literature study regression.
     The sweep angle was keep constant allong the cord of the tail for now, there is a option to change this within the code. (sweep_05_cord_v)"""
    Mission = 1  # Mission phase 0=transition, 1=Cruise
    number_vertical_tail = 2
    PRINT = True  # Printing the optimal values
    # Constraints
    min_span = 0.6
    min_AR_v = 1.0
    max_AR_v = 1.4
    max_sweep = 30

    # Resolution
    accuracy = 2
    steps = 100
    iteration = 100
    # intergration space
    AR_v = np.arange(0.5, 2.5, 0.1)
    sweep_v = np.arange(0, 45, 1)
    taper_v = np.arange(0.4, 1.1, 0.1)

    # base imports.
    # new imports
    X_rot_aft = 2.120  # FIXT IMPORTS
    Dv = 0.5  # FIXT IMPORTS
    CL_h = [1.4741, 0.337]     #[CL_h(max),CL_h(cruise)]

    # Horizontail tail dimensions
    AR_h = 6.8  # FIXT IMPORTS
    bh = 2.3  # FIXT IMPORTS
    M = 0.12  # FIXT IMPORTS

    # initial starting values
    tail_volume = 0.035 / number_vertical_tail  # FIXT IMPORTS
    C_eta_beta = 0.058  # FIXT IMPORTS

    Sh = np.sqrt(AR_h * bh)
    AR_w = b ** 2 / S
    Beta = np.sqrt(1 - M ** 2)


    # Empty list set
    Surface = []
    Span = []
    Moment_arm = []
    root_cord = []
    boom_length = []
    num = 0
    count = 0

    # C_eta_beta
    C_eta_beta_w = CL[Mission] ** 2 / (4 * np.pi * AR_w) + CL_h[Mission] ** 2 / (4 * np.pi * AR_h) * (Sh * bh) / (S * b)
    fuse_volume = 4 / 3 * np.pi * l_fus / 2 * (b_max / 2) ** 2
    C_eta_beta_fuse = -2 / (S * b) * fuse_volume

    for q in range(len(taper_v)):
        Surface_p = []
        span_p = []
        moment_arm_p = []
        Cv_r_p = []
        boom_length_p = []
        for p in range(len(AR_v)):
            Surface_k = []
            span_k = []
            moment_arm_k = []
            Cv_r_k = []
            boom_length_k = []

            for j in range(len(sweep_v)):
                Cv_r_r = []
                Surface_r = []
                span_r = []
                moment_arm_r = []
                boom_length_r = []

                bv = 0.6
                Sv = bv ** 2 / AR_v[p]
                lv = tail_volume * S * b / Sv

                Cv_r = 2 / (1 + taper_v[q]) * (Sv / bv)

                sweep_05_cord_v = sweep_v[j]  # for now a constant sweep is assumed
                X_LEMAC_v = bv / 6 * ((1 + 2 * taper_v[q]) / (1 + taper_v[q])) * np.tan(sweep_v[j] * deg2rad)
                CL_v_alpha = (Cl_alpha_v * AR_v[p]) / (2 + np.sqrt(
                    4 + (((AR_v[p] * Beta) / eta) ** 2) * (((np.tan(sweep_05_cord_v * deg2rad)) ** 2 / Beta ** 2) + 1)))

                # The range for minimum and maximum length for the moment arm
                laa = X_rot_aft - Xcg + Dv + Cv_r * 1.2  # MINIMUM MOMENT ARM
                lAA = 0 + Xcg
                la = round((laa * 1.01), accuracy)
                lA = round((5.75 - lAA), accuracy)  # MAXIMUM MOMENT ARM
                differ = (lA - la) / steps
                lb = np.arange(la, lA, differ)

                # print(count)
                for r in range(steps):
                    count = count + 1
                    for k in range(iteration):
                        # if Sv > 0 :

                        bv = np.sqrt(AR_v[p] * Sv)

                        Cv_r = 2 / (1 + taper_v[q]) * (Sv / bv)
                        Cv_bar = 2 / 3 * Cv_r * ((1 + taper_v[q] + taper_v[q] ** 2) / (1 + taper_v[q]))

                        # Updated values
                        lv = lb[r] - Cv_r + X_LEMAC_v + 0.25 * Cv_bar

                        # Update tail surface
                        # dsigma_dbeta =
                        Sv = 1 / number_vertical_tail * (C_eta_beta - C_eta_beta_fuse - C_eta_beta_w) / (CL_v_alpha) * (
                                    S * b) / lv  # * 1/(Vv_V**2)
                        num = num + 1

                    # List of all values
                    boom_length_r.append(lb[r])
                    Cv_r_r.append(Cv_r)
                    Surface_r.append(Sv)
                    span_r.append(bv)
                    moment_arm_r.append(lv)
                boom_length_k.append(boom_length_r)
                Cv_r_k.append(Cv_r_r)
                Surface_k.append(Surface_r)
                span_k.append(span_r)
                moment_arm_k.append(moment_arm_r)
            boom_length_p.append(boom_length_k)
            Cv_r_p.append(Cv_r_k)
            Surface_p.append(Surface_k)
            span_p.append(span_k)
            moment_arm_p.append(moment_arm_k)
        boom_length.append(boom_length_p)
        root_cord.append(Cv_r_p)
        Surface.append(Surface_p)
        Span.append(span_p)
        Moment_arm.append(moment_arm_p)

    PRINT = True

    # weight for optimization
    weight_surf = 1.2
    weight_AR = 0.7
    weight_span = -0.2
    weight_boom = 0.6
    weight_root_cord = 1.9
    weight_taper = -0.1

    tel = 0
    H = 10 ** 9
    # NORMILIZATION
    Norm_Surface = 1 / max(max(max(max(Surface))))
    Norm_Span = 1 / max(max(max(max(Span))))
    Norm_AR = 1 / max(AR_v)
    Norm_boom = 1 / max(max(max(max(boom_length))))
    Norm_root_cord = 1 / max(max(max(max(root_cord))))
    Norm_taper = 1 / max(taper_v)

    # OPTIMIZATION
    for q in range(len(taper_v)):
        for p in range(len(AR_v)):
            for j in range(len(sweep_v)):
                for r in range(steps):
                    if Span[q][p][j][r] > min_span:
                        if AR_v[p] > min_AR_v:

                            H_opt = weight_taper * (taper_v[q] * Norm_taper) + weight_root_cord * (
                                        root_cord[q][p][j][r] * Norm_root_cord) + weight_surf * Surface[q][p][j][
                                        r] * Norm_Surface + weight_AR * AR_v[p] * Norm_AR + weight_boom * \
                                    boom_length[q][p][j][r] * Norm_boom + weight_span * Span[q][p][j][r] * Norm_Span
                            if H_opt <= H:
                                tel = tel + 1
                                H = H_opt
                                q_opt = q
                                p_opt = p
                                j_opt = j
                                r_opt = r

    optimal_surface_v = Surface[q_opt][p_opt][j_opt][r_opt]
    optimal_span_v = Span[q_opt][p_opt][j_opt][r_opt]
    optimal_moment_arm = Moment_arm[q_opt][p_opt][j_opt][r_opt]
    optimal_root_cord = root_cord[q_opt][p_opt][j_opt][r_opt]
    optimal_sweep_v = sweep_v[j_opt]
    optimal_AR_v = AR_v[p_opt]
    optimal_boom = boom_length[q_opt][p_opt][j_opt][r_opt]
    optimal_taper = taper_v[q_opt]

    S_v = optimal_surface_v
    b_v = optimal_span_v
    l_v = optimal_moment_arm
    root_cord_v = optimal_root_cord
    sweep_v = optimal_sweep_v
    Aspect_v = optimal_AR_v
    Taper_v = optimal_taper
    cg_trailing_edge_boom = optimal_boom

    if PRINT == True:
        print('----------------------------------------------------------------')
        print('')
        print('The optimal values for a single vertical tail are:')
        print('Surface area :', optimal_surface_v, 'm^2')
        print('Span :', optimal_span_v, 'm')
        print('Moment arm :', optimal_moment_arm, 'm')
        print('Root cord :', optimal_root_cord, 'm')
        print('Sweep angle:', optimal_sweep_v, 'deg')
        print('Aspect ratio:', optimal_AR_v)
        print('boom length =', optimal_boom, 'm')
        print('taper ratio =', optimal_taper, 'm')
        print('-----------------------------------------------------------------')

    print("")
    print('-------------------------------------------------')
    print('rotor position =', X_rot_aft)
    print('update weights')

    return S_v,b_v,l_v,root_cord_v,sweep_v,Aspect_v,Taper_v,cg_trailing_edge_boom