import DSE
import numpy as np
from DSE import const
from DSE.aero import aero_constants
from DSE.stability import rotor_placement as rp



x_tail = 4.14
x_lemac = 2.181
mac = aero_constants.c_bar
x_ac_wing = x_lemac+0.25*mac
x_rot_front, x_rot_rear = rp.Rotor_Front_X, rp.Rotor_Rear_X
c_wing_root = aero_constants.c_root
l_engine = 0.400
l_alternator = 0.060
l_clutch = 0.032
l_prop = 0.18
l_aeroboom = 0.5
l_fus = x_lemac - l_aeroboom + l_alternator + c_wing_root + l_engine + l_clutch
print('L_fus =', l_fus, 'm')


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
                        'WingL':[9,x_ac_wing,-1.226,0.17],
                        'WingR':[9,x_ac_wing,1.226,0.17],
                        # 'WingL':[9,2.0828,-1.226,0.17],
                        # 'WingR':[9,2.0828,1.226,0.17],
                        'Power Management module':[6,x_lemac + c_wing_root,0,-0.1],
                        'Emergency Battery':[2,0.8,0,0],
                        'Combustion Engine':[35,x_lemac + c_wing_root + l_engine/2,0,-0.165],
                        'Alternator':[9,x_lemac + c_wing_root + l_engine + l_alternator/2,0,-0.165],
                        'Clutch':[1,x_lemac + c_wing_root + l_engine + l_alternator + l_clutch/2,0,-0.165],
                        'Push-prop':[0.6,x_lemac+ c_wing_root + l_engine + l_alternator + l_clutch + l_prop/2,0,-0.165],
                        'ECU':[1,x_lemac + c_wing_root,0,-0.05],
                        'Fuselage structure/Shell':[3,l_fus/2 + l_aeroboom,0,-0.17],
                        'ESC1':[0.62,x_rot_front + 0.2 ,-1.15,-0.06],
                        'ESC2':[0.62,x_rot_rear - 0.2,-1.15,-0.06],
                        'ESC3':[0.62,x_rot_front + 0.2,1.15,-0.06],
                        'ESC4':[0.62,x_rot_rear - 0.2,1.15,-0.06],
                        'Emotor1':[1.68,x_rot_front,-1.15,-0.1],
                        'Emotor2':[1.68,x_rot_rear,-1.15,-0.1],
                        'Emotor3':[1.68,x_rot_front,1.15,-0.1],
                        'Emotor4':[1.68,x_rot_rear,1.15,-0.1],
                        'Liftprops1':[0.4,x_rot_front,-1.15,-0.1],
                        'Liftprops2':[0.4,x_rot_rear,-1.15,-0.1],
                        'Liftprops3':[0.4,x_rot_front,1.15,-0.1],
                        'Liftprops4':[0.4,x_rot_rear,1.15,-0.1],
                        'Tailboom L':[7,(x_tail - x_rot_front)/2 + x_rot_front,-1.15,-0.1],
                        'Tailboom R':[7,(x_tail - x_rot_front)/2 + x_rot_front,1.15,-0.1],
                        'Hor Tail':[4,x_tail,0,0.79],
                        'Vert Tail 1':[1,x_tail,-1.15,0.163333333],
                        'Vert Tail 2':[1,x_tail,1.15,0.163333333],
                        'Fuel Pump 1':[0.07,x_ac_wing,-0.1,-0.12],
                        'Fuel Pump 2':[0.07,x_ac_wing,0.1,-0.12],
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
def extreme_cg_calc():
    test1 = cg_per_mission(True, False, False, False, False)
    test2 = cg_per_mission(True, False, True, False, False)
    test3 = cg_per_mission(True, False, False, False, True)
    test4 = cg_per_mission(True, False, True, False, True)
    #print('cg results',test1[1][0],test2[1][0],test3[1][0],test4[1][0])
    Xcg_max = max(test1[1][0], test2[1][0], test3[1][0], test4[1][0])
    return Xcg_max


print('extreme',extreme_cg_calc(), x_lemac)
def avionics_data(components_dict=components_dict):
    avionics_list = ['Air data boom', 'Iridium Antenna', 'Iridium Satellite communication module', 'GNSS Antenna 1',
                     'GNSS Antenna 2', 'GNSS Antenna 3', 'ADSB Transponder', 'Landing Visual Sensor', 'Flight Computer + Short Range Transceiver', 'Polar IMU + Magnetometer + GNSS', 'Flight Data Recorder']

    avionics_cg = cg_estimation(avionics_list, components_dict)
    return avionics_cg


engine_w, engine_xcg = np.array(components_dict['Combustion Engine'])[:2]
generator_w, generator_xcg = np.array(components_dict['Alternator'])[:2]
avionics_w, avionics_xcg = avionics_data()[0], avionics_data()[1][0]
wingL_w, wingL_xcg = np.array(components_dict['WingL'])[:2]
wingR_w, wingR_xcg = np.array(components_dict['WingR'])[:2]
main_structure_fuselage_w, main_structure_fuselage_xcg = np.array(components_dict['Fuselage structure/Shell'])[:2]
weight_fuselage = [engine_w, generator_w, avionics_w, wingL_w, wingR_w, main_structure_fuselage_w]
xcg_fuselage = [engine_xcg, generator_xcg, avionics_xcg, wingL_xcg, wingR_xcg, main_structure_fuselage_xcg]
fuselage_w = np.sum(weight_fuselage)
fuselage_xcg = np.sum(np.array(weight_fuselage)*np.array(xcg_fuselage))/fuselage_w
horizontalTail_w, horizontalTail_xcg = np.array(components_dict['Hor Tail'])[:2]
verticalTail1_w, verticalTail1_xcg = np.array(components_dict['Vert Tail 1'])[:2]
verticalTail2_w, verticalTail2_xcg = np.array(components_dict['Vert Tail 2'])[:2]
rotor_aft_1_w, rotor_aft_1_xcg = np.array(components_dict['Liftprops1'])[:2]
rotor_aft_2_w, rotor_aft_2_xcg = np.array(components_dict['Liftprops2'])[:2]
rotor_aft_3_w, rotor_aft_3_xcg = np.array(components_dict['Liftprops3'])[:2]
rotor_aft_4_w, rotor_aft_4_xcg = np.array(components_dict['Liftprops4'])[:2]

