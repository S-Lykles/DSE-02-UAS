import numpy as np
import matplotlib.pyplot as plt
#from DSE.transition.Locations import class_two_cg_estimation

def shear_bending_diagram(length, distributed_loads, point_loads, point_moment, load_factor, title):
    positions = np.linspace(0, length, 1000)

    shear_force_distributed = np.zeros_like(positions)
    total_weighted_positions_shear = 0
    total_weights_shear = 0

    for start, end, intensity in distributed_loads:
        load_segment = np.where((positions >= start) & (positions <= end), -intensity * (positions - start), 0)
        shear_force_distributed += load_segment

        total_weighted_positions_shear += np.trapz(positions * load_segment, positions)
        total_weights_shear += np.trapz(load_segment, positions)

    shear_force_distributed *= 9.81 * load_factor

    shear_force_point = np.zeros_like(positions)

    for point_load in point_loads:
        point_load_position, point_load_magnitude = point_load
        shear_force_point += np.where(positions >= point_load_position, -point_load_magnitude, 0)

        total_weighted_positions_shear += point_load_position * point_load_magnitude
        total_weights_shear += point_load_magnitude

    shear_force_point *= 9.81 * load_factor

    total_shear_force = shear_force_distributed + shear_force_point

    point_moment_position, point_moment_magnitude = point_moment
    moment_arm = positions - point_moment_position
    point_moment_effect = -point_moment_magnitude * np.heaviside(moment_arm, 0.5)

    total_bending_moment = np.cumsum(total_shear_force) * (length / len(positions))
    total_bending_moment += point_moment_effect
    total_bending_moment -= total_bending_moment[-1]  # Subtract constant offset to ensure B.M. is zero at the end

    cg_position_shear = total_weighted_positions_shear / total_weights_shear
    cg_distance_from_start_shear = cg_position_shear - positions[0]

    plt.figure(figsize=(10, 6))

    # Plot Shear Force Diagram
    plt.subplot(2, 1, 1)
    plt.plot(positions, total_shear_force)
    plt.xlabel('Position along Fuselage (m)')
    plt.ylabel('Shear Force (N)')
    plt.title(f'Shear Force Diagram for PterUAS: {title} for n = {load_factor}')
    plt.grid()

    # Plot Bending Moment Diagram
    plt.subplot(2, 1, 2)
    plt.plot(positions, total_bending_moment)
    plt.xlabel('Position along Fuselage (m)')
    plt.ylabel('Bending Moment (Nm)')
    plt.title(f'Bending Moment Diagram for PterUAS: {title} for n = {load_factor}')
    plt.grid()

    plt.tight_layout()
    plt.show()

    max_abs_shear_force = np.max(np.abs(total_shear_force))
    max_abs_bending_moment = np.max(np.abs(total_bending_moment))

    print(f'Title: {title}, Valid for load factor: {load_factor}')
    print("Maximum Absolute Shear Force:", max_abs_shear_force, "N")
    print("Maximum Absolute Bending Moment:", max_abs_bending_moment, "Nm")
    print("CG distance from the Start of l:", cg_distance_from_start_shear, "m")
    print()
    return cg_distance_from_start_shear

#variables
n_pos = 3.8
n_neg = -1.52

l_fuselage = 2
l_wing = 0.5 * l_fuselage
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
w_prop = 0
w_bat = 2
w_avionics = 1
q_fuse = 66 / l_fuselage

## fuselage plots
#input
distributed_loads = np.array([[0, 2, q_fuse]])
point_loads = np.array([[l_wing, w_wing], [l_eng, w_eng], [l_gen, w_gen], [l_plmod1, 0.5*w_plmod],
                        [l_plmod2, 0.5*w_plmod], [l_prop, w_prop], [l_bat, w_bat], [l_avionics, w_avionics]])
moment_points_pos = np.array((l_wing, 200*n_pos))
moment_points_neg = np.array((l_wing, 200*n_neg))

#diagrams
cg_fuse = shear_bending_diagram(l_fuselage, distributed_loads, point_loads, moment_points_pos, n_pos, 'fuselage')
shear_bending_diagram(l_fuselage, distributed_loads, point_loads, moment_points_neg, n_neg, 'fuselage')

## Boom plots
n_hov = 1.2

l_front_rot = 0.01
l_prop_spacing = 0.5
l_wing = cg_fuse + l_prop_spacing
l_rear_rot = l_wing + (l_fuselage-cg_fuse) + l_prop_spacing
l_tail = 3.99
l_boom =4

MTOW = 160/2
L_prop = -(1/2) * MTOW
w_prop = 5
w_tail = 5/2
w_boom = 5
w_wing_fw = -(2*w_prop + 0.5*w_tail + w_boom)
w_wing_hov = (MTOW + (w_wing_fw))
f_emp = -0.5*(567/9.81)
q_boom = w_boom / l_boom

#hover
distributed_loads_boom_hover = np.array([[0,l_boom,q_boom]])
point_loads_boom_hover = np.array([[l_front_rot, L_prop], [l_front_rot, w_prop], [l_wing, w_wing_hov], [l_rear_rot, L_prop], [l_rear_rot, w_prop],[l_tail,w_tail]])
moment_points_boom_hover_pos = np.array((l_wing, 70*n_hov))

shear_bending_diagram(l_boom, distributed_loads_boom_hover, point_loads_boom_hover, moment_points_boom_hover_pos, n_hov, 'boom in hover')

#fw flight
distributed_loads_boom_fw = np.array([[0,l_boom,q_boom]])
point_loads_boom_fw = np.array([[l_front_rot, w_prop], [l_wing, w_wing_fw], [l_rear_rot, w_prop],[l_tail,w_tail], [l_tail,-1]])
moment_points_boom_fw_pos = np.array((l_wing, 50*n_pos))
moment_points_boom_fw_neg = np.array((l_wing, 50*n_neg))

shear_bending_diagram(l_boom, distributed_loads_boom_fw, point_loads_boom_fw, moment_points_boom_fw_pos, n_pos, 'boom in fw flight')
shear_bending_diagram(l_boom, distributed_loads_boom_fw, point_loads_boom_fw, moment_points_boom_fw_neg, n_neg, 'boom in fw flight')
