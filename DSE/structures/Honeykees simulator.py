import numpy as np
import matplotlib.pyplot as plt

def calculate_centroid(R, r, t):
    """
    Calculates the centroid of a closed half circle shell which represents the upper fuselage (for clarification, ask Onni)
    :param R: Outer radius
    :param r: Inner radius
    :param t: Thickness of the flat bottom plate
    :return Cy: Y coordinate of the structure crossection centroid as measured from the top of the bottom plate
    :return Cy_halftube: Y coordinate of the open tube crossection, measured from the top of the bottom plate
    :return Cy_bottom: Y coordinate of the bottom plate centroid, measured from the top of the bottom plate
    :return A_halftube: area of the open half tube (without bottom plate)
    :return A_bottom: area of the bottom plate
    """
    Cy_halftube = 4/3/np.pi * (R**3-r**3) / (R**2-r**2)
    Cy_bottom = -t/2
    A_halftube = np.pi*(R**2-r**2) / 2
    A_bottom = 2*R*t
    Cy = (Cy_halftube*A_halftube + Cy_bottom*A_bottom) / (A_halftube+A_bottom)
    return Cy, Cy_halftube, Cy_bottom, A_halftube, A_bottom


def calculate_inertia(R, r, t):
    """
    Calculates moments of inertia of a half-circle shell around axis which are located at the flat surface in the middle

    :param R: Outer radius in m
    :param r: Inner radius in m
    :param t: thickness of the bottom plate
    :return Ixx: Moment of inertia around x axis through the centroid
    :return Iyy: Moment of inertia around y axis through the centroid
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
    Calculates the bending stress in each station along the crossection

    :param boom_co_arr: coordinates of the points along the crossection outline (in m)
    :param Mx: moment around the x-axis (in kN)
    :param My: moment around the y-axis (in kN)
    :param Ixx: second moment of area about the x-axis (in m^4)
    :param Iyy: second moment of area about the y-axis (in m^4)
    :return: stress_z: array wth stress at each point along the crossection outline
    """

    stress_z = Mx*Iyy*boom_co_arr[:,1]/Ixx + My*Ixx*boom_co_arr[:,0]/Iyy
    return stress_z

def place_booms(N_bottom, N_curve, R):
    """
    Distributes points equidistantly along the fuselage crossection
    :param N_bottom: number of points on the bottom plate
    :param N_curve: number of points along the curved section
    :param R: outer radius
    :return boom_co_arr: array with the coordinate pairs of each
    """

    angles = np.linspace(np.pi/(N_curve-1), np.pi*(N_curve-2)/(N_curve-1), N_curve-2)
    curve_boom_co_arr = np.column_stack((R * np.cos(angles), R * np.sin(angles)))

    bottom_boom_xco = np.linspace(-R, R, N_bottom)
    bottom_boom_yco = np.zeros(N_bottom)
    bottom_boom_co_arr = np.column_stack((bottom_boom_xco, bottom_boom_yco))

    boom_co_arr = np.vstack((curve_boom_co_arr, bottom_boom_co_arr))
    return boom_co_arr



R = 1
r = 0.99
t = 0.01
Mx = 50
My = 50

Ixx, Iyy = calculate_inertia(R, r, t)
boom_co_arr = place_booms(10, 15, R)
print(boom_co_arr)
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

#make the crossection actually circular
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')

# Show the plot
plt.show()




