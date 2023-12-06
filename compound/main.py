from power_curves.wong_tool import *
from compound.inputs import *
from aero.compound_helicopter import *
from power_curves.rotor_tool import *
from const import *
import numpy as np
from scipy.interpolate import interp1d


v_sl = v_stall(W, S, rho0, cl_max)
v_trans = 1.1 * v_sl

v_range_rot = np.linspace(0, v_trans, 1000)
#v_range_wing = np.linspace(40, 80, 1000)

rw_calc = True
fw_calc = True

if rw_calc:
    cd_parasite, d_airframe_wing, d_rotors, d_interference = parasite_drag(M_gross, S)
    R, D_v, omega, T_level, sig_max, sig_min = rotor_sizing_tool(W, DL, N, V_max, psi_rad=20*const.deg2rad, C_T_sig=0.11)
    P_ind = P_induced(v_range_rot, DL, W, k=1.15, k_dl=1.04)
    P_prof = P_profile_drag(v_range_rot, W, N, R, omega, sig_max, Cl_alpha_rot=5.73)
    P_par = 0.5 * rho0 * v_range_rot**3 * (d_airframe_wing + d_rotors + d_interference)     # = aeq = Cd * S

    P_tot_rot = P_ind + P_prof + P_par

    plt.figure(dpi=600)
    #plt.plot(v_range_rot, P_ind, label='Induced drag power')
    #plt.plot(v_range_rot, P_prof, label='Profile drag power')
    #plt.plot(v_range_rot, P_par, label='Parasitic drag power')
    plt.plot(v_range_rot, P_tot_rot, label='Tot drag power')

if fw_calc:
    cd_parasite, d_airframe_wing, d_rotors, d_interference = parasite_drag(M_gross,S)
    CL, CD = dragpolar_heli(b, S, Cl_start=0.2, Cl_end=cl_max, Cl_step=1000)
    P_tot_wing, v_range_wing = generate_Preq_ac(W, S, rho, CD, CL, eff_prop)

    plt.plot(v_range_wing, P_tot_wing, label='Wing power')

plt.title('Preq vs V Comparison')
#plt.plot(transition_param, smooth_transition, label='Smooth Transition')
plt.axvline(x=v_sl, linestyle='--', color='black', label='Stall speed of wing')
plt.axvline(x=1.1*v_sl, linestyle='--', color='green', label='Transition speed')
plt.xlabel('Velocity (m/s)')
plt.ylabel('Power Requirement (W)')
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.legend()
plt.grid()
plt.show()

print('Stall speed of wing', v_stall(W, S, rho0, cl_max))
print('Power required at stall', P_tot_wing[-1])

