import numpy as np
from plot_setting import *
import matplotlib.pyplot as plt
from power_curves.rotor_tool import rotor_sizing_tool, P_profile_drag, P_induced
from aero.cl_cd import dragpolar_dual
import dual_phase.inputs as inputs
import const


def plot_power_curves(DLs, bs, N, polar, CD0, S_design, k_dl=1.01, Ploss_frac=0.05, rotor_calculation=True, wing_calculation=True, name=None):
    plt.rcParams.update(report_tex)

    plt.figure(figsize=set_size(width=slidewidth*0.55, height=slideheight*0.55, subplots=(1, 1)))
    plt.axvline(const.v_cruise, label='Minimum Cruise Speed', color='k', linestyle='solid',linewidth=0.3)
    if rotor_calculation:
        for i, DL in enumerate(DLs):
        # Sizing the rotor based on W, DL, N and max velocity
            R, D_v, omega, T_level, sig_max = rotor_sizing_tool(const.MTOW, DL, N, const.v_cruise*1.5)


            # Setting up a linear space for the speed of the rotorcraft (Limits still need to be refined)
            v_rot = np.linspace(0, 35)

            # Calculating the different drag components, where for power loss this is 6% of other components (as in slides)
            P_p = P_profile_drag(v_rot, const.MTOW, N, R, omega, sig_max)
            P_i = P_induced(v_rot, DL, const.MTOW, k_dl=k_dl)
            P_par = 0.5 * const.rho0 * S_design * v_rot ** 3 * CD0
            P_loss = (P_p + P_i + P_par) * Ploss_frac

            # Calculating the total required power based on all power components
            P_req_rotor = P_p + P_i + P_par + P_loss

            # plt.plot(v_rot, P_p, label='Profile Drag')
            # plt.plot(v_rot, P_i, label='Induced Drag')
            # plt.plot(v_rot, P_par, label='Parasitic Drag')
            # plt.plot(v_rot, P_loss, label='Power losses')
            plt.plot(v_rot, P_req_rotor/1e3, label=f'DL={DL} [$\\mathrm{{N}}/\\mathrm{{m}}^2$]', c=plt.get_cmap('summer')((i+1)/(len(DLs)+2)))

    v_stall = 1000 # large number
    if wing_calculation:
        for i, (b, S) in enumerate(bs):
            CL, CD = polar(b, S)
            W1 = const.MTOW
            v1 = np.sqrt(W1 * 2 / (const.rho0 * S * CL))
            D = W1 * CD / CL
            P = D * v1 + const.P_aux
            
            plt.plot(v1[v1<55], P[v1<55]/1e3, label=f'b={b} [$\\mathrm{{m}}^2$], S={S} [$\\mathrm{{m}}^2$]', c=plt.get_cmap('autumn')((i)/(len(bs))))

            v_stall = min(v1[-1], v_stall)


    # if transitional_calculation:
    #     safety_factor = 1.1
    #     print(v_stall)
    #     v_transition = v_stall * safety_factor

    #     v_transition_band = 10  # Hey! we need to check this!!!

    plt.grid()
    # plt.subplots_adjust(right=0.7)

    plt.xlabel('Velocity [$\\mathrm{{m}}/\\mathrm{{s}}$]')
    plt.ylabel('Power [$\\mathrm{{kW}}$]')
    # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),fancybox=True)
    l = plt.legend(ncols=2)
    plt.ylim(bottom=-0.1)
    plt.xlim(left=0)
    plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
    plt.gca().grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
    plt.minorticks_on()
    plt.tight_layout()
    if name is not None:
        plt.savefig('power_plots/'+name, dpi=600)
    plt.show()
