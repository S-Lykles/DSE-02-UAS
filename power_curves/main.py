from .inputs import *
from .rotor_tool import generate_Preq_rotor, rotor_sizing_tool
from .wong_tool import generate_Preq_ac, find_optimum_range_and_endurance_speed
import matplotlib.pyplot as plt
from .hover_climb_power import *
from aero.cl_cd import *

#TIME SPAN OF SIMULATION FOR ROTOR AND FIXED WING CALCULATION
t_start_rot = 0
t_end_rot   = 40

t_start_ac  = 10
t_end_ac    = 80

step        = 1000

rotor_calc = True
ac_calc = True
Plot = True

single_point_v = 0
single_point_hover = HP  # Replace with the actual value for rotorcraft at single_point_v
single_point_hoverclimb = Clim_P    # Replace with the actual value for fixed-wing at single_point_v

# Initialize variables outside of if blocks
Preq_rotor = v_rot = Preq_ac = v_ac = None

if rotor_calc:
    R, D_v, omega, T_level, sig_max, sig_min = rotor_sizing_tool(DL, N)

    v_rot, P_p, P_i, P_par, _ = generate_Preq_rotor(A_eq, R, D_v, omega, T_level, sig_max, t_start_rot, t_end_rot, step)
    cl1, cd1 = dragpolar(b, S, 0, 1,1)
    CD0 = cd1[0]
    # print('CD0', CD0)
    # P_par = 0.5*rho*S*v_rot**3*CD0
    Preq_rotor = 1.04 * (P_p + P_i + P_par)

    if Plot:
        plt.figure(dpi=200)
        #plt.plot(v_rot, P_p, label='Profile Drag')
        #plt.plot(v_rot, P_i, label='Induced Drag')
        #plt.plot(v_rot, P_par, label='Parasitic Drag')
        plt.plot(v_rot, Preq_rotor, label='Total Power Required')

if ac_calc:
    # Ensure W, rho, S, AR, e, Cd0, eff_prop are defined before calling generate_Preq_ac
    Preq_ac, v_ac = generate_Preq_ac(W, S, rho, CD, CL, eff_prop)
    #Preq_ac, v_ac = generate_Preq_ac(W, rho, S, AR, e, Cd0, eff_prop, t_start_ac, t_end_ac, step)

    if Plot:
        plt.plot(v_ac, Preq_ac, label='Fixed Wing')
        pass

# Plot a single data point
#plt.scatter(single_point_v, single_point_hover, color='red', marker='o', label='Hover power')
plt.scatter(single_point_v, single_point_hoverclimb, color='blue', marker='x', label='Hover climb power')

if Plot:
    plt.title('Preq vs V Comparison')
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Power Requirement (W)')
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.legend()
    plt.grid()
    plt.savefig('s=3.7b=6.png')
    plt.show()

# Now, you can access these variables outside the if blocks
print("Minimum Rotor Power Requirement:", Preq_rotor.min() if Preq_rotor is not None else "N/A")
print("Minimum Fixed Wing Power Requirement:", Preq_ac.min() if Preq_ac is not None else "N/A")
print()
print('optimum rotor only', find_optimum_range_and_endurance_speed(Preq_rotor, v_rot))
print('optimum fixed wing', find_optimum_range_and_endurance_speed(Preq_ac, v_ac))
print()

print('Hover power clim', Clim_P)
generate_number_of_blades(R, sig_max)
