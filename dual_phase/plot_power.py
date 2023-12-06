import matplotlib.pyplot as plt
from power_curves.rotor_tool import rotor_sizing_tool, P_profile_drag, P_induced
from power_curves.wong_tool import generate_Preq_ac, find_optimum_range_and_endurance_speed
from aero.cl_cd import dragpolar
from dual_phase.inputs import *

rotor_calculation = True
ac_calculation = True
transitional_calculation = True
Plot = True

# These values vv are assumed as of now, reminder to have a discussion about this!
b = 6
S = 3.76

global P_req_rotor, v_rot
global P_req_ac, v_ac

if rotor_calculation:
    # Sizing the rotor based on W, DL, N and max velocity
    R, D_v, omega, T_level, sig_max = rotor_sizing_tool(W, DL, N, V_max)

    # Getting the cl/cd from Aero codes
    cl, cd = dragpolar(b, S, 0, 1,1)
    cd_zero = cd[0]

    # Setting up a linear space for the speed of the rotorcraft (Limits still need to be refined)
    v_rot = np.linspace(0, 60)

    # Calculating the different drag components, where for power loss this is 6% of other components (as in slides)
    P_profile_drag = P_profile_drag(v_rot, W, N, R, omega, sig_max)
    P_induced = P_induced(v_rot, DL, W)
    P_par = 0.5 * const.rho0 * S * v_rot ** 3 * cd_zero
    P_loss = 0.06 * (P_profile_drag + P_induced + P_par)

    # Calculating the total required power based on all power components
    P_req_rotor = P_profile_drag + P_induced + P_par + P_loss

    if Plot:
        plt.figure(dpi=200)
        plt.plot(v_rot, P_profile_drag, label='Profile Drag')
        plt.plot(v_rot, P_induced, label='Induced Drag')
        plt.plot(v_rot, P_par, label='Parasitic Drag')
        plt.plot(v_rot, P_loss, label='Power losses')
        plt.plot(v_rot, P_req_rotor, label='Total Power Required')

if ac_calculation:
    cl, cd = dragpolar(b, S, 0.1, 1, 1000)
    P_req_ac, v_ac = generate_Preq_ac(W, S, cd, cl, eta)

    if Plot:
        plt.plot(v_ac, P_req_ac, label='Fixed Wing')
        pass

if transitional_calculation:
    safety_factor = 1.1
    v_stall = v_ac[-1]
    print(v_stall)
    v_transition = v_stall * safety_factor

    v_transition_band = 10  # Hey! we need to check this!!!

plt.grid()
plt.legend()
plt.show()
