import DSE
import numpy as np
from DSE import const
#
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
                      'payload_supply': [50, 0.75 * l_fus, 0, 0],
                      'payload_relay':  [20, 0.75 * l_fus, 0, 0],
                      'avionics':       [10.24, 0.2 * l_fus, 0, 0],
                      'fuel_supply':    [13.7, 0.3 * l_fus, 0, 0],
                      'fuel_relay':     [38.7, 0.3 * l_fus, 0, 0],
                      'power_plant':    [28.5, 0.5 * l_fus, 0, 0]}

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
        component_list = ['payload_supply', 'avionics', 'power_plant']
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
        x_cg_tail = 0.9 * l_fus
        y_cg_tail = 0  # Assumed to be in the middle of the fuselage
        z_cg_tail = 0  # Assumed to be in the middle of the fuselage

    calc_emp = np.array([w_tail, x_cg_tail, y_cg_tail, z_cg_tail])
    return calc_emp


def propulsion_group_cg(class_one):
    """
    The output of the calculation is in the following form:
    calc_wing = [W_wing, x_cg_wing, y_cg_wing, z_cg_wing]
    """
    if class_one:
        n_rotors = 4
        w_prop_i = 2

        x_cg_prop_1 = x_cg_prop_2 = 0.2 * l_fus
        x_cg_prop_3 = x_cg_prop_4 = 0.8 * l_fus
        y_cg_prop = 0  # Assumed to be in the middle of the fuselage
        z_cg_prop = 0  # Assumed to be in the middle of the fuselage

        w_prop = 2 * n_rotors
        x_cg_prop = (2 * w_prop_i * x_cg_prop_1 + 2 * w_prop_i * x_cg_prop_3)/w_prop

    calc_emp = np.array([w_prop, x_cg_prop, y_cg_prop, z_cg_prop])
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
    prop_group_array = propulsion_group_cg(True)

    group_matrix = np.array([fuselage_group_array, wing_group_array, tail_group_array, prop_group_array])

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
    # Don't forget about Ixy, Ixz, Iyz

    return w_total, [x_cg, y_cg, z_cg], sensitivity  # [Ixx, Iyy, Izz]

#
# print(class_two_cg_estimation(True, False, False, False))

"""After this a weight-cg diagram can be drawn by varying weight conditions"""


#Updated version of the above code
# components_dict_old = {'component name': '[W, x_cg, y_cg, z_cg]',
#                    'WingL': [10, 1.5, 1.2, 0.85],
#                    'WingR': [10, 1.5, -1.2, 0.85],
#                    'BoomL': [5, 2.75, 1, 0.75],
#                    'BoomR': [5, 2.75, -1, 0.75],
#                    'Engine': [30, 2.5, 0, 0.725],
#                    'LiftMotorLF': [5, 0.5, 1, 0.85],
#                    'LiftMotorRF': [5, 0.5, -1, 0.85],
#                    'LiftMotorLR': [5, 3.5, 1, 0.85],
#                    'LiftMotorRR': [5, 3.5, -1, 0.85],
#                    'Generator': [4, 2.5, -0.2, 0.725],
#                    'Tailsurface': [5, 5.85, 0, 1.15],
#                    'Fuselage': [5, 1.5, 0, 0.725],
#                    'FCU': [0.1, 0.5, 0, 0.725],
#                    'Comms': [1, 0.5, 0, 0.725],
#                    'Battery': [2, 1, 0, 0.725],
#                     'Payload_pd': [50,1.5,0,0.3],
#                     'Payload_le': [20,1.5,0,0.3],
#                    'Fuel_pd': [14,1.5, 0, 0.3],
#                    'Fuel_le': [34,1.5, 0, 0.3]}


components_dict = {'component name': '[W, x_cg, y_cg, z_cg]',
                       'Air data boom':[0.142,0.23,0,0],
                        'Iridium Antenna':[0.1,0.5225,-0.027,0],
                        'Iridium Satellite communication module':[0.73,0.6875,0,0],
                        'GNSS Antenna 1':[0.011,0.56,0.00925,0],
                        'GNSS Antenna 2':[0.011,0.56,0.0185,0],
                        'GNSS Antenna 3':[0.011,0.56,0.02775,0],
                        'ADSB Transponder':[0.05,0.612,0.05,0],
                        'Landing Visual Sensor':[0.0215,0.765,0.035,-0.0275],
                        'Flight Computer + Short Range Transceiver':[0.1,0.74,-0.015,0.0425],
                        'Polar IMU + Magnetometer + GNSS':[0.17,0.645,0,0.0575],
                        'Flight Data Recorder':[0.07,0.735,0.029,0.056],
                        'WingL':[9,2.2731,-1.226,0.17],
                        'WingR':[9,2.2731,1.226,0.17],
                        # 'WingL':[9,2.0828,-1.226,0.17],
                        # 'WingR':[9,2.0828,1.226,0.17],
                        'Power Management module':[6,1.6401,0,-0.1],
                        'Emergency Battery':[2,0.8,0,0],
                        'Combustion Engine':[35,1.8401,0,-0.165],
                        'Alternator':[9,2.0701,0,-0.165],
                        'Clutch':[1,2.1201,0,-0.165],
                        'Push-prop':[0.6,2.2201,0,-0.165],
                        'ECU':[1,1.6401,0,-0.05],
                        'Fuselage structure/Shell':[3,1.845,0,-0.17],
                        'ESC1':[0.62,0.1,-1.15,-0.06],
                        'ESC2':[0.62,2.6,-1.15,-0.06],
                        'ESC3':[0.62,0.1,1.15,-0.06],
                        'ESC4':[0.62,2.6,1.15,-0.06],
                        'Emotor1':[1.68,1.8603,-1.15,-0.1],
                        'Emotor2':[1.68,3.2583,-1.15,-0.1],
                        'Emotor3':[1.68,1.8603,1.15,-0.1],
                        'Emotor4':[1.68,3.2583,1.15,-0.1],
                        'Liftprops1':[0.4,1.8603,-1.15,-0.1],
                        'Liftprops2':[0.4,3.2583,-1.15,-0.1],
                        'Liftprops3':[0.4,1.8603,1.15,-0.1],
                        'Liftprops4':[0.4,3.2583,1.15,-0.1],
                        'Tailboom L':[7,2.65,-1.15,-0.1],
                        'Tailboom R':[7,2.65,1.15,-0.1],
                        'Hor Tail':[4,5.1,0,0.79],
                        'Vert Tail 1':[1,5,-1.15,0.163333333],
                        'Vert Tail 2':[1,5,1.15,0.163333333],
                        'Fuel Pump 1':[0.07,1.215,-0.1,-0.12],
                        'Fuel Pump 2':[0.07,1.0875,0.1,-0.12],
                        'Emergency Battery':[2,0.8,0,0],
                        'Payload_pd': [50,1.5,0,0.3],
                        'Payload_le': [20,1.5,0,0.3],
                       'Fuel_pd': [14,1.5, 0, 0.3],
                       'Fuel_le': [34,1.5, 0, 0.3]
                                      }



def OEW(oew = True):

    if oew:
        # oew_data = ['WingL', 'WingR', 'BoomL', 'BoomR', 'Generator', 'Tailsurface', 'Fuselage', 'FCU', 'Comms']

        oew_data = ['Air data boom', 'Iridium Antenna', 'Iridium Satellite communication module', 'GNSS Antenna 1', 'GNSS Antenna 2', 'GNSS Antenna 3', 'ADSB Transponder', 'Landing Visual Sensor', 'Flight Computer + Short Range Transceiver', 'Polar IMU + Magnetometer + GNSS', 'Flight Data Recorder', 'WingL', 'WingR', 'Power Management module', 'Emergency Battery', 'Combustion Engine', 'Alternator', 'Clutch', 'Push-prop', 'ECU', 'Fuselage structure/Shell', 'ESC1', 'ESC2', 'ESC3', 'ESC4', 'Emotor1', 'Emotor2', 'Emotor3', 'Emotor4', 'Liftprops1', 'Liftprops2', 'Liftprops3', 'Liftprops4', 'Tailboom L', 'Tailboom R', 'Hor Tail', 'Vert Tail 1', 'Vert Tail 2']

    return oew_data


def payload(payload_supply, payload_relay, components):

    if payload_supply:
        payload_data = ['Payload_pd']
    elif payload_relay:
        payload_data = ['Payload_le']
    else:
        prop_data = []

    for element in prop_data:
        components.append(element)

    return components


def propulsion(fuel_supply, fuel_relay, components):

    if fuel_supply:
        # prop_data = ['Engine', 'Generator', 'Fuel_pd']
        prop_data = ['Fuel Pump 1', 'Fuel_pd']
        # return prop_data
    elif fuel_relay:
        # prop_data = ['Engine', 'Generator', 'Fuel_le', 'Battery']
        prop_data = ['Fuel Pump 2', 'Fuel_le','Emergency Battery']
    else:
        prop_data = []

    for element in prop_data:
        components.append(element)

    return components

def cg_estimation(component_list, component_dict):
    # component_list = ['empty_fus', 'payload_supply', 'avionics', 'power_plant']
    first_element = component_list[0]
    component_matrix = np.array([component_dict[first_element]])
    component_list = component_list[1:]

    for component in component_list:
        component_matrix = np.concatenate((component_matrix, np.array([component_dict[component]])))

    weight_of_groups = component_matrix[:, 0]
    x_cg_groups = component_matrix[:, 1]
    y_cg_groups = component_matrix[:, 2]
    z_cg_groups = component_matrix[:, 3]

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
    # Don't forget about Ixy, Ixz, Iyz

    return w_total, [x_cg, y_cg, z_cg], sensitivity  # [Ixx, Iyy, Izz]

def cg_per_mission(oew, fuel_supply, fuel_relay, payload_supply, payload_relay, components_dict=components_dict ):

    # centre_of_gravity(OEM, fuel_pd, 0, 0, 0)

    component_list = OEW(oew)
    component_list = propulsion(fuel_supply, fuel_relay, component_list)
    component_list = propulsion(payload_supply, payload_relay, component_list)

    cal_cg = cg_estimation(component_list, components_dict)

    return cal_cg
def extreme_cg_calc(cal_cg=cg_per_mission.cal_cg):
    test1 = cg_per_mission(True, False, False, False, False)
    test2 = cg_per_mission(True, False, True, False, False)
    test3 = cg_per_mission(True, False, False, False, True)
    test4 = cg_per_mission(True, False, True, False, True)
    #print('cg results',test1[1][0],test2[1][0],test3[1][0],test4[1][0])
    Xcg_max = max(test1, test2, test3, test4)
    return Xcg_max


test1 = cg_per_mission(True, False, False, False, False)
test2 = cg_per_mission(True, False, True, False, False)
test3 = cg_per_mission(True, False, False, False, True)
test4 = cg_per_mission(True, False, True, False, True)
print('cg results',test1[1][0],test2[1][0],test3[1][0],test4[1][0])


def avionics_data(components_dict=components_dict):
    avionics_list = ['Air data boom', 'Iridium Antenna', 'Iridium Satellite communication module', 'GNSS Antenna 1',
                     'GNSS Antenna 2', 'GNSS Antenna 3', 'ADSB Transponder', 'Landing Visual Sensor', 'Flight Computer + Short Range Transceiver', 'Polar IMU + Magnetometer + GNSS', 'Flight Data Recorder']

    avionics_cg = cg_estimation(avionics_list, components_dict)
    return avionics_cg


engine_w, engine_xcg = np.array(components_dict['Combustion Engine'])[:2]
generator_w, generator_xcg = np.array(components_dict['Alternator'])[:2]
avionics_w, avionics_xcg = avionics_data()[0], avionics_data()[1][0]
wingL_w, wingL_xcg = np.array(components_dict['WingL'])[:2]
wingL_w, wingL_xcg = np.array(components_dict['WingR'])[:2]
main_structure_fuselage_w, main_structure_fuselage_xcg = np.array(components_dict['Fuselage structure/Shell'])[:2]
weight_fuselage = [engine_w, generator_w, avionics_w, wingL_w, wingR_w, main_structure_fuselage_w]
xcg_fuselage = [engine_xcg, generator_xcg, avionics_xcg, wingL_xcg, wingR_xcg, main_structure_fuselage_xcg]
fuselage_w = np.sum(weight_fuselage)
fuselage_xcg = np.sum(weight_fuselage*xcg_fuselage)/fuselage_w
horizontalTail_w, horizontalTail_xcg = np.array(components_dict['Hor Tail'])[:2]
verticalTail1_w, verticalTail1_xcg = np.array(components_dict['Vert Tail 1'])[:2]
verticalTail2_w, verticalTail2_xcg = np.array(components_dict['Vert Tail 2'])[:2]
rotor_aft_1_w, rotor_aft_1_xcg = np.array(components_dict['Liftprops1'])[:2]
rotor_aft_2_w, rotor_aft_2_xcg = np.array(components_dict['Liftprops2'])[:2]
rotor_aft_3_w, rotor_aft_3_xcg = np.array(components_dict['Liftprops3'])[:2]
rotor_aft_4_w, rotor_aft_4_xcg = np.array(components_dict['Liftprops4'])[:2]

from DSE.structures import center_of gravity as cg

eng_w = cg.engine.w
# fairings
# structure
#
