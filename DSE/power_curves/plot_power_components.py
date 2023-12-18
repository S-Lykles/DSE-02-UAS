import numpy as np
from DSE.plot_setting import *
import matplotlib.pyplot as plt
from DSE.power_curves.rotor_tool import rotor_sizing_tool, P_profile_drag, P_induced
from DSE.aero.cl_cd import dragpolar_dual
import DSE.dual_phase.inputs as inputs
from DSE import const



def plot_power_components(DL, N, polar, CD0, S_design, k_dl=1.01, Ploss_frac=0.05, name=None):
    plt.rcParams.update(report_tex)


    plt.figure(figsize=set_size(textwidth, subplots=(1, 1)))
    # Sizing the rotor based on W, DL, N and max velocity
    R, D_v, omega, T_level, sig_max = rotor_sizing_tool(const.MTOW, DL, N, const.v_cruise*1.5)


    # Setting up a linear space for the speed of the rotorcraft (Limits still need to be refined)
    v_rot = np.linspace(0, 40)

    # Calculating the different drag components, where for power loss this is 6% of other components (as in slides)
    P_p = P_profile_drag(v_rot, const.MTOW, N, R, omega, sig_max)
    P_i = P_induced(v_rot, DL, const.MTOW, k_dl=k_dl)
    P_par = 0.5 * const.rho0 * S_design * v_rot ** 3 * CD0
    P_loss = (P_p + P_i + P_par) * Ploss_frac

    # Calculating the total required power based on all power components
    P_req_rotor = P_p + P_i + P_par + P_loss

    plt.plot(v_rot, P_p/1e3, label='Profile Power')
    plt.plot(v_rot, P_i/1e3, label='Induced Power')
    plt.plot(v_rot, P_par/1e3, label='Parasitic Power')
    plt.plot(v_rot, P_loss/1e3, label='Power losses', )
    plt.plot(v_rot, P_req_rotor/1e3, label=f'Total Power')

    plt.grid()
    plt.subplots_adjust(right=0.7)

    plt.xlabel('Velocity $[\\mathrm{{m}}/\\mathrm{{s}}]$')
    plt.ylabel('Power $[\\mathrm{{kW}}]$')
    # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),fancybox=True)
    l = plt.legend(loc='upper left', ncols=2)
    l.set_zorder(20)
    plt.ylim(0, plt.ylim()[1]*1.3)
    plt.xlim(left=0)
    plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
    plt.gca().grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
    plt.minorticks_on()
    plt.tight_layout()
    if name is not None:
        plt.savefig('power_plots/'+name, dpi=600)
    plt.show()
