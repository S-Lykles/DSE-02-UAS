import DSE
import numpy as np
from DSE import const

"""The following code is based on the aircraft design book by Roskam: 
airplane-design-part-v-component-weight-estimation."""

""" Assumptions:
For the initial version of this code, many bold assumptions are made in class one:
See the comments behind the values in the functions!
- Weights are assumed to be in [kg]
"""
# parameters
l_ground_map = 6
l_fus = 6
# x_fus_offset = 6 - 6 No fuselage offset

fus_component_dict = {'component name': '[W, x_cg, y_cg, z_cg]',
                      'empty_fus':      [34.98, 0.45 * l_fus, 0, 0],
                      'payload_supply': [50, 0.45 * l_fus, 0, 0],
                      'payload_relay':  [20, 0.45 * l_fus, 0, 0],
                      'avionics':       [10.24, 0.2 * l_fus, 0, 0]}


def fuselage_group_cg(class_one, component_dict):
    """
    The output of the calculation is in the following form:
    calc_fuselage = [W_fus, x_cg_fus, y_cg_fus, z_cg_fus]
    """
    if class_one:
        w_fus = 34.98  # Assumed based on midterm report
        x_cg_fus = 0.45 * l_fus  # Assumed to be 0.45 of the length of the fuselage as by Roskam
        y_cg_fus = 0  # Assumed to be in the middle
        z_cg_fus = 0  # Assumed to be in the middle

    if not class_one:
        # Selecting the components that are part of the fuselage
        component_list = ['payload_supply', 'avionics']
        component_matrix = np.array([component_dict['empty_fus']])
        for component in component_list:
            component_matrix = np.concatenate((component_matrix, np.array([component_dict[component]])))

        # cg calculations
        w_components = component_matrix[:, 0]
        w_fus = np.sum(w_components)

        x_cg_comp = component_matrix[:, 1]
        y_cg_comp = component_matrix[:, 2]
        z_cg_comp = component_matrix[:, 3]

        x_cg_fus = np.sum(w_components * x_cg_comp) / w_fus
        y_cg_fus = np.sum(w_components * y_cg_comp) / w_fus
        z_cg_fus = np.sum(w_components * z_cg_comp) / w_fus

        pass

    calc_fuselage = [w_fus, x_cg_fus, y_cg_fus, z_cg_fus]
    return calc_fuselage


def wing_group_cg(class_one, x_lemac, c_bar):
    """
    The output of the calculation is in the following form:
    calc_wing = [W_wing, x_cg_wing, y_cg_wing, z_cg_wing]
    """
    if class_one:
        w_wing = 51.88
        x_cg_wing = x_lemac + 0.4 * c_bar
        y_cg_wing = 0  # Assumed to be in the middle of the fuselage
        z_cg_wing = 0  # Assumed to be in the middle of the fuselage

    calc_wing = np.array([w_wing, x_cg_wing, y_cg_wing, z_cg_wing])
    return calc_wing


def tail_group_cg(class_one):
    """
    The output of the calculation is in the following form:
    calc_wing = [W_wing, x_cg_wing, y_cg_wing, z_cg_wing]
    """
    if class_one:
        w_tail = 3.68
        x_cg_tail = 2
        y_cg_tail = 0  # Assumed to be in the middle of the fuselage
        z_cg_tail = 0  # Assumed to be in the middle of the fuselage

    calc_emp = np.array([w_tail, x_cg_tail, y_cg_tail, z_cg_tail])
    return calc_emp


# Parameters like fuel percentages can be added
def class_two_cg_estimation(empty, fuel, payload_supply, payload_relay):
    """Inputs are the different groups, these values vary with certain parameters.
    We decompose the center of gravity calculations with respect to main groups:
    -----------------------------------------------------------------------------
    Fuselage group
    Wing group
    Empennage group
    Propulsion group
    Fixed equipment group (like landing gear)
    -----------------------------------------------------------------------------
    The datum for this calculation is the start of the 6x6 area boundary
    """
    fuselage_group_array = fuselage_group_cg(False, fus_component_dict)
    wing_group_array = wing_group_cg(True, x_lemac=3, c_bar=2)
    tail_group_array = tail_group_cg(True)

    group_matrix = np.array([fuselage_group_array, wing_group_array, tail_group_array])

    weight_of_groups = group_matrix[:, 0]
    x_cg_groups = group_matrix[:, 1]
    y_cg_groups = group_matrix[:, 2]
    z_cg_groups = group_matrix[:, 3]

    w_total = np.sum(weight_of_groups)
    x_cg = np.sum(weight_of_groups * x_cg_groups) / w_total
    y_cg = np.sum(weight_of_groups * y_cg_groups) / w_total
    z_cg = np.sum(weight_of_groups * z_cg_groups) / w_total

    # Possible implementation to assign names to the values
    # Multiple array to find the offset of the total cg
    sensitivity = weight_of_groups / w_total

    # These are initial estimates of the Inertia parameters, based om the composite moment of inertia rule
    # Mass moment of inertia of the different groups have to estimated go give an accurate value
    # Ixx = np.sum(0) + np.sum(wing_group_array * ((y_cg_groups - y_cg)**2 + (z_cg_groups - z_cg)**2))
    # Iyy = np.sum(0) + np.sum(wing_group_array * ((z_cg_groups - z_cg)**2 + (x_cg_groups - x_cg) **2))
    # Izz = np.sum(0) + np.sum(wing_group_array * ((x_cg_groups - x_cg)**2 + (y_cg_groups - y_cg) **2))

    return w_total, [x_cg, y_cg, z_cg], sensitivity  # [Ixx, Iyy, Izz]


print(class_two_cg_estimation(True, False, False, False))

"""After this a weight-cg diagram can be drawn by varying weight conditions"""
