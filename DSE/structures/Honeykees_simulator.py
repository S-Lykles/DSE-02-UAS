import numpy as np
import matplotlib.pyplot as plt

def calculate_centroid(R, r, t):
    """
    Calculates the centroid of a half circle
    :param R:
    :param r:
    :param t:
    :return:
    """
    Cy_halftube = 4/3/np.pi * (R**3-r**3) / (R**2-r**2)
    Cy_bottom = -t/2
    A_halftube = np.pi*(R**2-r**2) / 2
    A_bottom = 2*R*t
    Cy = (Cy_halftube*A_halftube + Cy_bottom*A_bottom) / (A_halftube+A_bottom)
    #print(Cy_halftube, Cy_bottom)
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

# def calculate_shear_stress(Sy, Ixx, coordinates, booms):
#     c = -Sy / Ixx
#     delta_q = booms * coordinates[:, 1] * c
#     s_tot = np.sum(delta_q)
#     return delta_q, s_tot

def fuse_shear(fy, ixx, t, r, coord, ybar):
    xco = coord[:, 0]
    yco = coord[:, 1]
    theta = np.arctan2(yco, xco)

    # Calculate circle_shear only for xco between -r and 0
    circle_shear = np.where((yco > 0) & (xco >= -r) & (xco <= 0),
                            -(fy * t * r / ixx) * (-r * np.cos(theta) + r - ybar * theta), 0)

    # Calculate mirror values for circle_shear and concatenate
    circle_shear_mirror = np.where((yco > 0) & (xco <= r) & (xco >= 0),
                                   -(fy * t * r / ixx) * (r * np.cos(theta) + r - ybar * theta), 0)
    circle_shear = np.concatenate((circle_shear, circle_shear_mirror))

    # Calculate straight_shear only for xco between 0 and r
    straight_shear = np.where((yco == 0) & (xco >= 0) & (xco <= r),
                              -(fy * t * ybar * xco) / ixx, 0)

    # Ensure all arrays have the same length by padding with zeros
    max_length = max(len(xco), len(yco), len(circle_shear), len(straight_shear))
    xco = np.pad(xco, (0, max_length - len(xco)))
    yco = np.pad(yco, (0, max_length - len(yco)))
    circle_shear = np.pad(circle_shear, (0, max_length - len(circle_shear)))
    straight_shear = np.pad(straight_shear, (0, max_length - len(straight_shear)))

    return xco, yco, circle_shear, straight_shear

def place_booms(N_bottom, N_curve, R):
    angles = np.linspace(np.pi/(N_curve-1), np.pi*(N_curve-2)/(N_curve-1), N_curve-2)
    curve_boom_co_arr = np.column_stack((R * np.cos(angles), R * np.sin(angles)))

    bottom_boom_xco = np.linspace(-R, R, N_bottom)
    bottom_boom_yco = np.zeros(N_bottom)
    bottom_boom_co_arr = np.column_stack((bottom_boom_xco, bottom_boom_yco))

    boom_co_arr = np.vstack((curve_boom_co_arr, bottom_boom_co_arr))
    return boom_co_arr

# inputs
R = 0.4
r = 0.4 - 0.01
t = 0.01
Mx = 50
My = 50
Sy = -100
s_range = np.linspace(0,0.4,100)
angle_range = np.linspace(0, np.pi/2, 100)

Cy, Cy_halftube, Cy_bottom, A_halftube, A_bottom = calculate_centroid(R, r, t)
Ixx, Iyy = calculate_inertia(R, r, t)
boom_co_arr = place_booms(50, 150, R)

# ##shear
xco, yco, circle_shear, straight_shear = fuse_shear(Sy, Ixx, t, R, boom_co_arr, Cy)

# plt.figure()
# plt.plot(xco, circle_shear)
# plt.plot(xco, straight_shear)
# plt.show()

plt.scatter(xco, yco, c=straight_shear, cmap='viridis')
#plt.scatter(xco, yco, c=circle_shear, cmap='viridis')
plt.colorbar(label='Circle Shear')
plt.xlabel('X Coordinates')
plt.ylabel('Y Coordinates')
plt.title('2D Scatter Plot with Circle Shear Colors')
plt.show()


# ##Bending stress
# stress_z = calculate_bending_stress(boom_co_arr, Mx, My, Ixx, Iyy)
#
# #plot bending
# # Plot the 3D curve
# plt.scatter(boom_co_arr[:,0], boom_co_arr[:,1], c=stress_z, cmap='viridis', label='Axial Stress')
#
# # Add a color bar to show the stress values
# cbar = plt.colorbar()
# cbar.set_label('Axial Stress')
#
# # Customize the plot
# plt.xlabel('X-axis')
# plt.ylabel('Y-axis')
# plt.title('3D Curve with Axial Stress')
#
# # Show the plot
# plt.show()
