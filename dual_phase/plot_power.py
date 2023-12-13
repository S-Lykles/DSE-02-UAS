import numpy as np
import matplotlib.pyplot as plt
from power_curves.rotor_tool import rotor_sizing_tool, P_profile_drag, P_induced
from power_curves.wong_tool import generate_Preq_ac, find_optimum_range_and_endurance_speed
from aero.cl_cd import dragpolar_dual
import dual_phase.inputs as inputs
import const

rotor_calculation = True
wing_calculation = True
transitional_calculation = True
Plot = True

# These values vv are assumed as of now, reminder to have a discussion about this!

# DL = 500
N = 4
b = 6
S = 3.763
k_dl = 1.01 # lower than normal because the rotors are in free air
vc = 1
Ploss_frac = 0.05

# global P_req_rotor, v_rot
# global P_req_ac, v_ac

plt.figure()
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.Set1.colors)
DLs = [300, 400, 500, 600, 700]
if rotor_calculation:
    for i, DL in enumerate(DLs):
    # Sizing the rotor based on W, DL, N and max velocity
        R, D_v, omega, T_level, sig_max = rotor_sizing_tool(const.MTOW, DL, N, const.v_cruise*1.5)

        # Getting the cl/cd from Aero codes
        cl, cd = dragpolar_dual(b, S, 0,1,1)
        cd_zero = cd[0]

        # Setting up a linear space for the speed of the rotorcraft (Limits still need to be refined)
        v_rot = np.linspace(0, 30)

        # Calculating the different drag components, where for power loss this is 6% of other components (as in slides)
        P_p = P_profile_drag(v_rot, const.MTOW, N, R, omega, sig_max)
        P_i = P_induced(v_rot, DL, const.MTOW, k_dl=k_dl)
        P_par = 0.5 * const.rho0 * S * v_rot ** 3 * cd_zero
        P_loss = (P_p + P_i + P_par) * (1+Ploss_frac)

        # Calculating the total required power based on all power components
        P_req_rotor = P_p + P_i + P_par + P_loss

        if Plot:
            # plt.plot(v_rot, P_p, label='Profile Drag')
            # plt.plot(v_rot, P_i, label='Induced Drag')
            # plt.plot(v_rot, P_par, label='Parasitic Drag')
            # plt.plot(v_rot, P_loss, label='Power losses')
            plt.plot(v_rot, P_req_rotor, label=f'DL={DL}', c=plt.get_cmap('summer')((i+1)/(len(DLs)+4)))

bs = [(6, 3.763), (5, 3.5), (4, 3.2), (3, 2.7)]
v_stall = 1000 # large number
if wing_calculation:
    for i, (b, S) in enumerate(bs):
        CL, CD = dragpolar_dual(b, S, CL_start=0.1, CL_end=1.5, CL_step=1000)
        W1 = const.MTOW
        v1 = np.sqrt(W1 * 2 / (const.rho0 * S * CL))
        D = W1 * CD / CL
        P = D * v1 + const.P_aux
        
        if Plot:
            plt.plot(v1[v1<60], P[v1<60], label=f'b={b}, S={S}', c=plt.get_cmap('autumn')((i)/(len(bs))))

        v_stall = min(v1[-1], v_stall)


# if transitional_calculation:
#     safety_factor = 1.1
#     print(v_stall)
#     v_transition = v_stall * safety_factor

#     v_transition_band = 10  # Hey! we need to check this!!!

plt.axvline(const.v_cruise, label='Cruise Speed', color='k', linestyle='--',linewidth=0.5)
plt.grid()
plt.subplots_adjust(right=0.7)

plt.xlabel('Velocity [m/s]')
plt.ylabel('Power [W]')
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),fancybox=True)
plt.legend()
plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
plt.gca().grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.7)
plt.minorticks_on()
plt.show()
