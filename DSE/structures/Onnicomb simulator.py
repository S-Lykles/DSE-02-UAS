import numpy as np
import matplotlib.pyplot as plt
import math

def calculate_centroid(R, r, t):
    Cy_halftube = 4/3/np.pi * (R**3-r**3) / (R**2-r**2)
    Cy_bottom = -t/2
    A_halftube = np.pi*(R**2-r**2) / 2
    A_bottom = 2*R*t
    Cy = (Cy_halftube*A_halftube + Cy_bottom*A_bottom) / (A_halftube+A_bottom)
    print(Cy_halftube, Cy_bottom)
    return Cy, Cy_halftube, Cy_bottom, A_halftube, A_bottom


def calculate_inertia(R, r, t):
    """
    Calculates moments of inertia of a half-circle shell around axis which are located at the flat surface in the middle

    :param R: Outer radius in m
    :param r: Inner radius in m
    :return Ixx: Moment of inertia around x axis
    :return Iyy: Moment of inertia around y axis
    """
    Cy, Cy_halftube, Cy_bottom, A_halftube, A_bottom = calculate_centroid(R, r, t)

    Ixx_halftube = np.pi*(R**4-r**4)/8 + A_halftube*(Cy-Cy_halftube)**2
    Ixx_bottom = 2*R*t**3/12 + A_bottom*(Cy-Cy_bottom)**2
    Ixx = Ixx_halftube + Ixx_bottom

    Iyy_halftube = np.pi*(R**4-r**4)/8
    Iyy_bottom = (2*R)**3*t/12
    Iyy = Iyy_halftube + Iyy_bottom

    return Ixx, Iyy


def calculate_bending_stress(boom_co_arr, Mx, My, Ixx, Iyy):
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

    stress_z = Mx*Iyy*boom_co_arr[:,1]/Ixx + My*Ixx*boom_co_arr[:,0]/Iyy
    return stress_z

def place_booms(N_bottom, N_curve, R):
    angles = np.linspace(np.pi/(N_curve-1), np.pi*(N_curve-2)/(N_curve-1), N_curve-2)
    curve_boom_co_arr = np.column_stack((R * np.cos(angles), R * np.sin(angles)))

    bottom_boom_xco = np.linspace(-R, R, N_bottom)
    bottom_boom_yco = np.zeros(N_bottom)
    bottom_boom_co_arr = np.column_stack((bottom_boom_xco, bottom_boom_yco))

    boom_co_arr = np.vstack((curve_boom_co_arr, bottom_boom_co_arr))
    return boom_co_arr[:,0], boom_co_arr[:,1]



R = 1
r = 0.99
t = 0.01
Mx = 50
My = 50


Ixx, Iyy = calculate_inertia(R, r, t)
boom_co_arr = place_booms(50, 150, R)
stress_z = calculate_bending_stress(boom_co_arr, Mx, My, Ixx, Iyy)

# Plot the 3D curve
plt.scatter(boom_co_arr[:,0], boom_co_arr[:,1], c=stress_z, cmap='viridis', label='Axial Stress')

# Add a color bar to show the stress values
cbar = plt.colorbar()
cbar.set_label('Axial Stress')

# Customize the plot
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('3D Curve with Axial Stress')

# Show the plot
plt.show()



# def boom_coor(N_bottom, N_curve, R):
#     boom_co = np.empty((N_curve // 2 + 1, 2))
#     dt = np.pi / N_curve
#     t = 0
#     for i in range(N_curve // 2 + 1):
#         x = R * np.cos(t)
#         y = R * np.sin(t)
#         boom_co[i] = [x, y]
#         t = t + dt
#     return boom_co

def calculate_shear_stress(Sy, Ixx, coordinates, booms):
    c = -Sy / Ixx
    delta_q = np.array([])
    for i, y in enumerate(coordinates[:, 1]):
        dq = booms[i] * y * c
        delta_q = np.append(delta_q, [dq])
    s_tot = sum(delta_q)
    return delta_q, s_tot

# x, y = place_booms(5, 7, 0.4)
#
# plt.figure()
# plt.plot(x, y, marker='x')
# plt.show()

#delta, tot = calculate_shear_stress(100, 10000, boom_coor(0, 5, 1), booms)
def test(sy, t, r, ixx, theta, ybar):
    return -(sy * t * r / ixx) * (-r * np.cos(theta) + r - ybar * theta)

def test2(sy, t, ixx, s, ybar):
    return -(sy * t * ybar * s) / ixx

trange = np.linspace(0, np.pi/2, 100)
srange = np.linspace(0, 0.4, 100)
plt.figure()
plt.plot(np.linspace(0,0.4,100), test(10000, 0.02, 0.4, 0.0001, trange, 0.17))
plt.plot(np.linspace(0,0.4,100), test2(10000, 0.02, 0.0001, srange, 0.17))
plt.show()

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



