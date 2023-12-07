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

#for cl_max in [1.4,1.5,1.6,1.7,1.8]:
 #   v_sl = v_stall(W, S, rho0, cl_max)

if rw_calc:
    plt.figure(dpi=600)
    for DL in [120]:
        cd_parasite, d_airframe_wing, d_rotors, d_interference = parasite_drag(M_gross, S)
        R, D_v, omega, T_level, sig_max = rotor_sizing_tool(W, DL, N, V_max, psi_rad=20*const.deg2rad, C_T_sig=0.11)
        P_ind = P_induced(v_range_rot, DL, W, k=1.15, k_dl=1.04)
        P_prof = P_profile_drag(v_range_rot, W, N, R, omega, sig_max, Cl_alpha_rot=5.73)
        P_par =( 0.5 * rho0 * v_range_rot**3 * cd_parasite*S )/eff_prop#(d_airframe_wing + d_rotors + d_interference)     # = aeq = Cd * S
        P_loss = 0.05*(P_ind + P_prof + P_par)
        P_tot_rot = P_ind + P_prof + P_par + P_loss

        DP = delta_p_climb(vc, W)
        climb_power = DP + P_tot_rot[0]
        print(f'Climb power at DL={DL}', climb_power)

        #plt.plot(v_range_rot, P_ind, label='Induced drag power')
        #plt.plot(v_range_rot, P_prof, label='Profile drag power')
        #plt.plot(v_range_rot, P_par, label='Parasitic drag power')
        plt.plot(v_range_rot, P_tot_rot, label=f'Rotor power @ DL={DL}')
        plt.scatter(0.2, climb_power, marker='x', label=f'Climb power (Vc=1m/s) @ DL={DL}')

if fw_calc:
    for S,b in [(3.70,6), (2.57, 5), (1.64, 4)]:
        cd_parasite, d_airframe_wing, d_rotors, d_interference = parasite_drag(M_gross,S)
        CL, CD = dragpolar_heli(b, S, Cl_start=0.2, Cl_end=cl_max, Cl_step=1000)
        P_tot_wing, v_range_wing = generate_Preq_ac(W, S, CD, CL, eff_prop)

        plt.plot(v_range_wing, P_tot_wing, label=f'Wing power b={b}, s={S}')
        print(f'Power required at stall (S,b)={S,b}', P_tot_wing[-1])

plt.title('Power required curve for compound helicopter')
#plt.plot(transition_param, smooth_transition, label='Smooth Transition')
plt.axvline(x=41.111, linestyle='--', color='black', label='Cruise Speed')
#plt.axvline(x=1.1*v_sl, linestyle='--', color='green', label='Transition speed')
#plt.axvline(x=1.15*v_sl, linestyle='--', color='black', label='Transition Complete')
plt.xlabel('Velocity (m/s)')
plt.ylabel('Power Requirement (W)')
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.legend()
plt.grid()
plt.show()

print()
print('Stall speed of wing', v_stall(W, S, rho0, cl_max))
#print('Power required at stall', P_tot_wing[-1])

print(N, W, DL)
