import numpy as np
import matplotlib.pyplot as plt
import math

def calculate_inertia(R, r):
    """
    Calculates moments of inertia of a half-circle shell around axis which are located at the flat surface in the middle

    :param R: Outer radius in m
    :param r: Inner radius in m
    :return Ixx: Moment of inertia around x axis
    :return Iyy: Moment of inertia around y axis
    """
    Ixx = np.pi*(R**4-r**4)/8
    Iyy = Ixx
    return Ixx, Iyy


def calculate_bending_stress(pos_x, pos_y, Mx, My, Ixx, Iyy):
    """
    Calculate bending stress.

    Parameters:
    - pos_x = x coordinates of the postion where the stress is computed
    - pos_y = y coordinates ofthe position where the stress is computed
    - Fz = Force along the z axis
    - Mx = Moment around x axis
    - My = Moment around y axis
    - Ixx = moment of inertia about x axis
    - Iyy = moment of inertia about y axis

    Returns:
    - stress_z = axial stress along z-axis
    """

    stress_z = Mx*Iyy*pos_y/Ixx + My*Ixx*pos_x/Iyy
    return stress_z

def place_booms(N_bottom, N_curve, R):
    angles = np.linspace(np.pi/(N_curve-1), np.pi*(N_curve-2)/(N_curve-1), N_curve-2)
    curve_boom_co_arr = np.column_stack((R * np.cos(angles), R * np.sin(angles)))

    bottom_boom_xco = np.linspace(-R, R, N_bottom)
    bottom_boom_yco = np.zeros(N_bottom)
    bottom_boom_co_arr = np.column_stack((bottom_boom_xco, bottom_boom_yco))

    boom_co_arr = np.vstack((curve_boom_co_arr, bottom_boom_co_arr))
    return boom_co_arr

def boom_coor(N_bottom, N_curve, R):
    boom_co = np.empty((N_curve // 2 + 1, 2))
    dt = np.pi / N_curve
    t = 0
    for i in range(N_curve // 2 + 1):
        x = R * np.cos(t)
        y = R * np.sin(t)
        boom_co[i] = [x, y]
        t = t + dt
    return boom_co

def calculate_shear_stress(Sy, Ixx, coordinates, booms):
    c = -Sy / Ixx
    delta_q = np.array([])
    for i, y in enumerate(coordinates[:, 1]):
        dq = booms[i] * y
        delta_q = np.append(delta_q, [dq*c])
    q_tot = sum(delta_q)
    return delta_q, q_tot

booms = np.array([100, 100, 200])

delta, tot = calculate_shear_stress(100, 10000, boom_coor(0, 5, 1), booms)
print(tot)

def calculate_normal_stress(axial_force, area):
    """
    Calculate normal stress.

    Parameters:
    - axial_force: Axial force (N)
    - area: Cross-sectional area (m^2)

    Returns:
    - Normal stress (Pa)
    """
    return axial_force / area



