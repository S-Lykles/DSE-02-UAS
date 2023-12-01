from inputs import *
from rotor_tool import generate_Preq_rotor, rotor_sizing_tool
from wong_tool import *
import matplotlib.pyplot as plt
from Hover_Climb_Power import *

#TIME SPAN OF SIMULATION FOR ROTOR AND FIXED WING CALCULATION
t_start_rot = 3
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
    R, D_v, omega, T_level, sig_max = rotor_sizing_tool(DL, N)
    Preq_rotor, v_rot = generate_Preq_rotor(A_eq, R, D_v, omega, T_level, sig_max, t_start_rot, t_end_rot, step)

    if Plot:
        plt.figure(dpi=600)
        plt.plot(v_rot, Preq_rotor, label='Rotorcraft')

if ac_calc:
    # Ensure W, rho, S, AR, e, Cd0, eff_prop are defined before calling generate_Preq_ac
    Preq_ac, v_ac = generate_Preq_ac(W, S, rho, CD, CL, eff_prop)

    if Plot:
        plt.plot(v_ac, Preq_ac, label='Fixed Wing')

# Plot a single data point
plt.scatter(single_point_v, single_point_hover, color='red', marker='o', label='Hover power')
plt.scatter(single_point_v, single_point_hoverclimb, color='blue', marker='x', label='Hover climb power')

if Plot:
    plt.title('Preq vs V Comparison')
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Power Requirement (W)')
    plt.xlim(0)
    plt.ylim(0)
    plt.legend()
    plt.grid()
    plt.show()

# Now, you can access these variables outside the if blocks
print("Minimum Rotor Power Requirement:", Preq_rotor.min() if Preq_rotor is not None else "N/A")
print("Minimum Fixed Wing Power Requirement:", Preq_ac.min() if Preq_ac is not None else "N/A")

print('optimum rotor only', find_optimum_range_and_endurance_speed(Preq_rotor, v_rot))
print('optimum fixed wing', find_optimum_range_and_endurance_speed(Preq_ac, v_ac))
