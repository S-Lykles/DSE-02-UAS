print('Structures <3 !!!')
from DSE import const
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def aero_data_to_numpy(file_name):
    file_dir = Path(__file__).parent
    data = np.loadtxt(file_dir/f"{file_name}")
    data = np.transpose(data)
    y = data[0]
    chord = data[1]
    cl = data[3]
    return y, chord, cl


def find_nearest_point(point, array):
    difference_array = array - point
    index = np.argmin(abs(difference_array))
    value = array[index]
    return value, index


def distribution_from_aero_data(start, stop, step, data_name, point_range):
    y, chord, cl = aero_data_to_numpy(data_name)
    cl_distribution = np.interp(point_range, y, cl)
    chord_distribution = np.interp(point_range, y, chord)
    lift_distribution = cl_distribution * 0.5 * const.rho0 * const.v_cruise**2 * (chord_distribution)

    # Lift only for the start and stop caps
    i_start = find_nearest_point(start, point_range)[1]
    i_stop = find_nearest_point(stop, point_range)[1]

    lift_distribution = lift_distribution[i_start:i_stop + 1]
    empty_array = np.zeros(len(point_range))
    empty_array[i_start:i_stop + 1] = lift_distribution
    lift_distribution = empty_array

    return lift_distribution


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


def combined_loading(beam_start, beam_stop, step):
    point_range = np.arange(beam_start, beam_stop + step, step)

    load1 = load_distribution(0, 4, step, 4, -3, 'linear', point_range)
    load2 = load_distribution(1, 2, step, 400, 400, 'linear', point_range)
    load3 = 0 * distribution_from_aero_data(0.2, beam_stop, step, 'wing_data.txt', point_range)
    load4 = point_load(2, 800, point_range)

    loading_distribution = load2 + load4
    plt.plot(point_range, loading_distribution)
    plt.show()
    return loading_distribution, point_range


loads, point_range = combined_loading(0, 3, 0.01)


def moment_distr_from_load_distr(load_distribution, point_range, step):
    moment_distribution = np.array([])

    # cumulative load distribution
    for i in np.arange(0, len(point_range), 1):
        distances = point_range - point_range[i]
        distances = distances[i:]
        loads = load_distribution[i:]
        print(loads)
        moment = np.trapz(loads * distances, distances, step)
        print(moment)
        moment_distribution = np.append(moment_distribution, moment/step)
    print(moment_distribution[0])
    plt.plot(point_range, moment_distribution)
    plt.show()

moment_distr_from_load_distr(loads, point_range, 0.01)

