import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumtrapz

def shear_bending_diagram(length, distributed_load_intensity, point_loads, load_factor):
    positions = np.linspace(0, length, 100)

    shear_force_distributed = np.zeros_like(positions)

    for start, end, intensity in distributed_loads:
        load_segment = np.where((positions >= start) & (positions <= end), -intensity * (positions - start), 0)
        shear_force_distributed += load_segment

    shear_force_distributed *= 9.81 * load_factor

    shear_force_point = np.zeros_like(positions)

    for point_load in point_loads:
        point_load_position, point_load_magnitude = point_load
        shear_force_point += np.where(positions >= point_load_position, -point_load_magnitude, 0)

    shear_force_point *= 9.81 * load_factor

    total_shear_force = shear_force_distributed + shear_force_point

    middle_index = len(positions) // 2

    bending_moment_positions_pos = positions[middle_index:]
    bending_moment_pos = cumtrapz(total_shear_force[middle_index:], x=bending_moment_positions_pos, initial=0)

    bending_moment_positions_neg = positions[:middle_index][::-1]
    bending_moment_neg = -cumtrapz(total_shear_force[:middle_index][::-1], x=bending_moment_positions_neg, initial=0)[
                          ::-1]

    bending_moment = np.concatenate((bending_moment_neg, bending_moment_pos))

    plt.plot(positions, total_shear_force)
    plt.xlabel('Position along Fuselage (m)')
    plt.ylabel('Shear Force (N)')
    plt.title('Shear Force Diagram for PterUAS')
    plt.grid()
    plt.show()

    plt.plot(positions, bending_moment)
    plt.xlabel('Position along Fuselage (m)')
    plt.ylabel('Bending moment (N/m)')
    plt.title('Bending Moment Diagram for PterUAS')
    plt.grid()
    plt.show()

    max_abs_shear_force = np.max(np.abs(total_shear_force))
    max_abs_bending_moment = np.max(np.abs(bending_moment))

    print(f'Valid for load factor: {load_factor}')
    print("Maximum Absolute Shear Force:", max_abs_shear_force, "N")
    print("Maximum Absolute Bending Moment:", max_abs_bending_moment, "N/m")
    print()


#variables
n_pos = 3.8
n_neg = -1.52

l_fuselage = 2
l_wing = 0.5 * l_fuselage
#l_cg = l_wing + 0.5
l_eng = l_fuselage - 0.4
l_gen = l_fuselage - 0.8
l_plmod1 = 0.25 * l_fuselage
l_plmod2 = 0.75 * l_fuselage
l_prop = l_fuselage
l_bat = 0.4
l_avionics = 0.3

w_wing = -160
w_eng = 35
w_gen = 4
w_plmod = 50
w_prop = 2
w_bat = 2
w_avionics = 1
q_fuse = 66 / l_fuselage

#input
distributed_loads = np.array([[0, 2, q_fuse]])
point_loads = np.array([[l_wing, w_wing], [l_eng, w_eng], [l_gen, w_gen], [l_plmod1, 0.5*w_plmod],
                        [l_plmod2, 0.5*w_plmod], [l_prop, w_prop], [l_bat, w_bat], [l_avionics, w_avionics]])

#diagrams
shear_bending_diagram(l_fuselage, distributed_loads, point_loads, n_pos)
shear_bending_diagram(l_fuselage, q_fuse, point_loads, n_neg)

