from inputs import *
from rotor_tool import generate_Preq_rotor, rotor_sizing_tool
from wong_tool import generate_Preq_ac
import matplotlib.pyplot as plt

#ROTOR DESIGN PARAMETERS
DL = 200
N = 1
A_eq = 0.09

t_start_rot = 5
t_end_rot   = 60

t_start_ac  = 30
t_end_ac    = 100

step        = 0.01


#R, D_v, omega, T_level, sig_max = rotor_sizing_tool(DL, N)
#Preq_rotor, v_rot = generate_Preq_rotor(A_eq, R, D_v, omega, T_level, sig_max, t_start_rot, t_end_rot, step)

Preq_ac, v_ac = generate_Preq_ac(W, rho, S, AR, e, Cd0, t_start_ac, t_end_ac, step)

Plot_rotor = False
Plot_ac = True

if Plot_rotor == True:
    plt.figure(dpi=600)
    plt.plot(v_rot, Preq_rotor)
    plt.title('Preq vs V rotorcraft forward flight')
    plt.grid()
    #plt.legend()
    plt.show()

if Plot_ac == True:
    plt.figure(dpi=600)
    plt.plot(v_ac, Preq_ac)
    plt.title('Preq vs V fixed wing forward flight')
    plt.grid()
    #plt.legend()
    plt.show()