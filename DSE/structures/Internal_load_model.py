print('Structures <3 !!!')
import numpy as np
data = np.loadtxt(r"C:\Users\tijme\PycharmProjects\DSE2023\DSE\structures\wing_data")
print(data)

def find_nearest_point(point, array):
    difference_array = array - point
    index = np.argmin(abs(difference_array))
    value = array[index]
    return value, index


def distribution_from_data(start, stop, data_file, point_range):
    pass


def load_distribution(start, stop, step, m1, m2, type, point_range):
    # matching selected points
    start, start_i = find_nearest_point(start, point_range)
    stop, stop_i = find_nearest_point(stop, point_range)

    # creating distribution
    if type == 'linear':
        force_array = np.zeros([len(point_range)])
        a = (m2 - m1)/(stop - start)
        points = np.arange(0, stop - start + step, step)
        distribution = a * points + m1
        force_array[start_i:stop_i + 1] = distribution
        print(force_array)


def point_load(point, magnitude, point_range):
    """
    This code returns a point load, as distributed on a certain data point
    """
    load_lst = np.zeros(len(point_range))
    nearest_point = find_nearest_point(point, point_range)[1]
    load_lst[nearest_point] = magnitude
    print(load_lst)
    return load_lst


def combined_loading(beam_start, beam_stop, step):
    point_range = np.arange(beam_start, beam_stop + step, step)
    load_distribution(0, 4, step, 4, -3, 'linear', point_range)
    load_distribution(2, 5, step, -2, 3, 'linear', point_range)
    point_load(3, 2, point_range)


combined_loading(0, 6, 0.01)


# def moment_distribution(point, magnitude, beam_start, point_range, load_distribution):
#     """
#     This code uses the load distribution as input to create a moment distribution
#     """
#     moment_distr = np.zeros(len(point_range))
#     i = 0
#     for i in range(0, len(point_range)):
#         moment_distr[i] = sum(load_distribution(:i)) / (point_range[i]- beam_start)
#
#         i = i + 1
#
#
#     nearest_point = find_nearest_point(point, point_range)[1]
#     moment_distr[nearest_point] = moment_distr[nearest_point] + magnitude
#
#     print(moment_distr)
#     return moment_distr

