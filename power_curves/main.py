from inputs import *
from rotor_tool import generate_Preq_rotor, rotor_sizing_tool
from wong_tool import generate_Preq_ac
import matplotlib.pyplot as plt


#ROTOR DESIGN PARAMETERS
DL = 236
N = 4
A_eq = 1 * 0.0929

#TIME SPAN OF SIMULATION FOR ROTOR AND FIXED WING CALCULATION
t_start_rot = 5
t_end_rot   = 40

t_start_ac  = 10
t_end_ac    = 100

step        = 1000


rotor_calc = True
ac_calc = True
Plot = True


if rotor_calc:
    R, D_v, omega, T_level, sig_max = rotor_sizing_tool(DL, N)
    Preq_rotor, v_rot = generate_Preq_rotor(A_eq, R, D_v, omega, T_level, sig_max, t_start_rot, t_end_rot, step)

    if Plot:
        plt.figure(dpi=600)
        plt.plot(v_rot, Preq_rotor, label='Rotorcraft')

if ac_calc:
    Preq_ac, v_ac = generate_Preq_ac(W, rho, S, AR, e, Cd0, eff_prop, t_start_ac, t_end_ac, step)

    if Plot:
        plt.plot(v_ac, Preq_ac, label='Fixed Wing')

if Plot:
    plt.title('Preq vs V Comparison')
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Power Requirement (W)')
    plt.legend()
    plt.grid()
    plt.show()

    print("Minimum Rotor Power Requirement:", Preq_rotor.min())
    print("Minimum Fixed Wing Power Requirement:", Preq_ac.min())
