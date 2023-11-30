from inputs import *
from rotor_tool import generate_Preq_rotor, rotor_sizing_tool
from wong_tool import generate_Preq_ac
from wong_tool import find_optimum_range_and_endurance_speed
import matplotlib.pyplot as plt

#ROTOR DESIGN PARAMETERS
DL = 200
N = 1
A_eq = 0.09

#TIME SPAN OF SIMULATION FOR ROTOR AND FIXED WING CALCULATION
t_start_rot = 5
t_end_rot   = 40

t_start_ac  = 10
t_end_ac    = 100

step        = 1


rotor_calc = False
ac_calc = True
Plot = True


if rotor_calc == True:

    R, D_v, omega, T_level, sig_max = rotor_sizing_tool(DL, N)
    Preq_rotor, v_rot = generate_Preq_rotor(A_eq, R, D_v, omega, T_level, sig_max, t_start_rot, t_end_rot, step)

    if Plot == True:
        plt.figure(dpi=600)
        plt.plot(v_rot, Preq_rotor)
        plt.title('Preq vs V rotorcraft forward flight')
        plt.grid()
        #plt.legend()
        plt.show()

if ac_calc == True:
    Preq_ac, v_ac = generate_Preq_ac(W, rho, S, AR, e, Cd0, eff_prop, t_start_ac, t_end_ac, step)
    speeds = find_optimum_range_and_endurance_speed(Preq_ac, v_ac)
    print(speeds)
    if Plot == True:
        plt.figure(dpi=600)
        plt.plot(v_ac, Preq_ac)
        plt.title('Preq vs V fixed wing forward flight')
        plt.grid()
        #plt.legend()
        plt.show()