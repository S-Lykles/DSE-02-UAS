print('Structures <3 !!!')
from DSE import const
import numpy as np
import matplotlib.pyplot as plt
import csv
from pathlib import Path

def aero_data_to_numpy(file_name):
    file_dir = Path(__file__).parent
    data = np.loadtxt(file_dir/f"{file_name}")
    data = np.transpose(data)
    y = data[0] #+ 1.15
    chord = data[1]
    cl = data[3] #* 0.5
    cd = (data[4] + data[5]) #* 0.5
    cm_geom = data[6] #*0.5
    return y, chord, cl, cd, cm_geom


def find_nearest_point(point, array):
    difference_array = array - point
    index = np.argmin(abs(difference_array))
    value = array[index]
    return value, index


def distribution_from_aero_data(start, stop, step, data_name, point_range):
    y, chord, cl, cd, cm_geom = aero_data_to_numpy(data_name)
    cl_distribution = np.interp(point_range, y, cl)
    cd_distribution = np.interp(point_range, y, cd)
    cm_distribution = np.interp(point_range, y, cm_geom)
    chord_distribution = np.interp(point_range, y, chord)
    max_th = 0.12 * chord_distribution
    print(const.v_cruise)
    lift_distribution = cl_distribution * 0.5 * const.rho0 * const.v_cruise**2 * (chord_distribution)
    drag_distribution = cd_distribution * 0.5 * const.rho0 * const.v_cruise**2 * (chord_distribution)
    torque_distribution = cm_distribution * 0.5 * const.rho0 * const.v_cruise**2 * (chord_distribution ** 2)
    # Lift only for the start and stop caps
    i_start = find_nearest_point(start, point_range)[1]
    i_stop = find_nearest_point(stop, point_range)[1]

    lift_distribution = lift_distribution[i_start:i_stop + 1]
    torque_distribution = torque_distribution[i_start:i_stop + 1]

    empty_array = np.zeros(len(point_range))
    empty_array[i_start:i_stop + 1] = lift_distribution
    lift_distribution = empty_array

    empty_array2 = np.zeros(len(point_range))
    empty_array2[i_start:i_stop + 1] = torque_distribution
    torque_distribution = empty_array2

    return lift_distribution*step, drag_distribution*step, torque_distribution, max_th


def load_distribution(start, stop, step, m1, m2, type, point_range):
    # matching selected points
    start, start_i = find_nearest_point(start, point_range)
    stop, stop_i = find_nearest_point(stop, point_range)

    # creating distribution
    if type == 'linear':
        load_array = np.zeros([len(point_range)])
        a = (m2 - m1)/(stop - start)
        points = np.arange(0, stop - start, step)
        distribution = a * points + m1
        load_array[start_i:stop_i] = distribution
        return load_array*step


def point_load(point, magnitude, point_range):
    """
    This code returns a point load, as distributed on a certain data point
    """
    load_lst = np.zeros(len(point_range))
    nearest_point = find_nearest_point(point, point_range)[1]
    load_lst[nearest_point] = magnitude
    # print(load_lst)
    return load_lst


def combined_loading(beam_start, beam_stop, step, VTOL):
    """
    For the WING !
    """
    point_range = np.arange(beam_start, beam_stop + step, step)
    zeros = np.zeros(len(point_range))
    lift_distribution, drag_distribution, torque_distribution, max_th = distribution_from_aero_data(0, beam_stop, step, 'external files/wing_data.txt', point_range)
    Boom_from_centerline = 1.15
    load_factor_vtol = 1.1
    Upwards_load = (load_factor_vtol * const.MTOW) / 2
    Upwards_pointload = point_load(Boom_from_centerline, Upwards_load, point_range)
    # reac_torque = load_distribution(0.2, Boom_from_centerline, step, -1700/step, -1700/step, 'linear', point_range)
    # reac_torque_aileron = load_distribution(0.2, 2.61, step, -22.25 / step, -22.25 / step, 'linear', point_range)
    # elevator_torque = load_distribution(1.15, 2.3, step, -0.5 * 0.118 * 263.429 / step, 0, 'linear', point_range)
    # elevator_torque_2 = load_distribution(0, 1.15, step, 0, -0.5 * 0.118 * 263.429 / step, 'linear', point_range)
    #vertical_tail_moment = point_load(0.88, 3.8 * 109.722 / 2 * 1.15 / step, point_range)
    #torque_from_horizontal_tail = load_distribution(0, 0.88, step, -3.8 * 5 / 2 * 1.15 / step, -3.8 * 5 / 2 * 1.15 / step, 'linear', point_range)
    #reac_force = load_distribution(0, Boom_from_centerline, step, 800, 800, 'linear', point_range)
    if VTOL:
        loading_distribution = Upwards_pointload
        loading_distributionx = 0
    else:
        load_factor_manouvre = 1.5
        loading_distributionz =  load_factor_manouvre * lift_distribution #+ vertical_tail_moment#+ reac_force
        torque_distribution =  3.8*torque_distribution #+ elevator_torque + elevator_torque_2 #+ torque_from_horizontal_tail  + reac_torque + reac_torque_aileron
        loading_distributionx = load_factor_manouvre * drag_distribution
    # ax = plt.axes(projection='3d')
    # ax.plot3D(point_range, loading_distributionx/step, zeros, 'blue')
    # ax.plot3D(point_range, zeros, loading_distributionz/step, 'red')
    # ax.set_xlabel('semi-spanwise location (m)')
    # ax.set_ylabel('internal force (drag) (N)')
    # ax.set_zlabel('internal force (lift) (N)')
    # plt.show()
    plt.plot(point_range, torque_distribution)
    plt.xlabel('semi-spanwise location (m)')
    plt.ylabel('Internal Torque (Nm)')
    plt.show()
    return loading_distributionx, loading_distributionz, torque_distribution, point_range, max_th


loadsx, loadsz, torqueyy, point_range, max_th = combined_loading(0, 1.15, 1/128, False)


def moment_distr_from_load_distr(load_distributionx, load_distributionz, point_range, step):
    momentx_distribution = np.array([])
    momentz_distribution = np.array([])
    zeros = np.zeros(len(point_range))
    # cumulative load distribution
    for i in np.arange(0, len(point_range), 1):
        distances = point_range - point_range[i]
        distances = distances[i:]
        loadsx = load_distributionx[i:]
        loadsz = load_distributionz[i:]
        momentx = np.trapz(loadsx * distances, distances, step)
        momentz = np.trapz(loadsz * distances, distances, step)
                   #+ load_distribution(0, 0.88, step, 3.8 * 109.722 / 2 * 1.15 / step, 3.8 * 109.722 / 2 * 1.15 / step, "linear", point_range)
        #print(momentx, momentz)
        momentx_distribution = np.append(momentx_distribution, momentx/step)
        momentz_distribution = np.append(momentz_distribution, momentz/ step) #+ point_load(0.87, 3.8 * 109.722 / 2 * 1.15, point_range)
    #print(momentx_distribution[0], momentz_distribution[0])
    ax = plt.axes(projection='3d')
    ax.plot3D(point_range, momentx_distribution, zeros, 'blue')
    ax.plot3D(point_range, zeros, momentz_distribution, 'red')
    ax.set_xlabel('semi-spanwise location (m)')
    ax.set_ylabel('internal moment (x) (Nm)')
    ax.set_zlabel('internal moment (z) (Nm)')
    plt.show()
    return momentx_distribution, momentz_distribution, point_range


def Ixxreq(moment_distribution, point_range, maxth, sigmacrit, e_mod, rho, step): #, A):
    Ixxreq = abs(moment_distribution)*max_th/ (sigmacrit) #- 109.722/ 2 / A)
    tskinreq = Ixxreq/ (2*(0.9*max_th/2)**2 * 0.5* max_th/0.12)
    intcomp = moment_distribution/maxth
    req_spar_space = tskinreq/np.sqrt(intcomp/(e_mod*0.9))
    across = tskinreq*0.5*max_th/0.12
    vbox = across*step
    mbox = sum(vbox)*rho
    # plt.plot(point_range, Ixxreq)
    # plt.xlabel('semi-spanwise location (m)')
    # plt.ylabel('required Ixx (m^4)')
    # plt.show()
    # plt.plot(point_range, tskinreq)
    # plt.xlabel('semi-spanwise location (m)')
    # plt.ylabel('required skin thickness (m)')
    # plt.show()
    # plt.plot(point_range, across)
    # plt.xlabel('semi-spanwise location (m)')
    # plt.ylabel('predicted cross-sectional area structural members (m^2)')
    # plt.show()
    # plt.plot(point_range, req_spar_space)
    # plt.xlabel('semi-spanwise location (m)')
    # plt.ylabel('required rib spacing (m)')
    # plt.show()
    print("The maximum required moment of inertia along the wing span of the wing box is:", max(Ixxreq))
    print(max(tskinreq))
    print(mbox)


moment_distributionx, momentdistributionz, point_range = moment_distr_from_load_distr(loadsx, loadsz, point_range, 1/128)
#print(Ixxreq(moment_distributionx, point_range, 0.001, 300*10**6, 71*10**9, 2800, 0.01)) #, 0.000448))


#print("The maximum moment experienced is:", momentdistributionz[0])


#Ixxreq(moment_distribution, point_range, max_th, 300.1*(10 ** 6), 20.1*(10**9),2800, 0.01)

