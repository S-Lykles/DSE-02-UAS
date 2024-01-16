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

    Ixx_halftube = np.pi*(R**4-r**4)/8 + A_halftube*Cy**2
    Ixx_bottom = 2*R*t**3/12 + A_bottom*(Cy-Cy_bottom)**2
    Ixx = Ixx_halftube + Ixx_bottom

    Iyy_halftube = np.pi*(R**4-r**4)/8
    Iyy_bottom = (2*R)**3*t/12
    Iyy = Iyy_halftube + Iyy_bottom

    return Ixx, Iyy, Cy

def calculate_bending_stress(Ixx, Iyy, Cy, boom_co_arr, Mx, My):
    """
    Calculates the bending stress in each station along the crossection

    :param boom_co_arr: coordinates of the points along the crossection outline (in m)
    :param Mx: moment around the x-axis (in kN)
    :param My: moment around the y-axis (in kN)
    :param Ixx: second moment of area about the x-axis (in m^4)
    :param Iyy: second moment of area about the y-axis (in m^4)
    :return: stress_z: array wth stress at each point along the crossection outline
    """

    stress_z = Mx*Iyy*(boom_co_arr[:,1]-Cy)/Ixx**2 + My*Ixx*boom_co_arr[:,0]/Iyy**2
    return stress_z

def calculate_torsion_stress(boom_co_arr, T, r, t):
    """
    :param boom_co_arr:
    :param T:
    :param r:
    :return:
    """

    stress_shear = T / (np.pi*r**2*t) * np.ones(len(boom_co_arr))
    return stress_shear

def calculate_shear_stress(Ixx, Iyy, Sx, Sy, N_curve, N_bottom, r):

    """
    Calculates the shear stress in the crossection
    :param Ixx: Second moment of aera about x axis
    :param Iyy: Second moment of inertia about y-axis
    :param Sx: Shear force in the x-direction
    :param Sy: Shear force in the y-direction
    :param N_curve: number of points on the boom (this is needed to transform the s-coordinate system to an xy-coordinate system)
    :param N_bottom: number of points on the bottom plate
    :param r: Outer radius
    :return stress_shear: array with the shear stresses which corresponds to the boom placement function
    """

    theta = np.linspace(0, np.pi / 2, int(N_curve/2))
    s2 = np.linspace(0, r, int(N_bottom/2))
    stress12 = t * r * (Sx * r / Iyy * (np.cos(theta) - 1) + Sy / Ixx * (r * np.sin(theta) - Cy * theta))
    stress2 = stress12[-1]
    stress23 = -Sx * t / Iyy * (r * s2 - s2 ** 2 / 2) - Sy * t / Ixx * Cy * s2 + stress2

    shear_flow = np.hstack((stress12[::-1], stress12, stress23, stress23[::-1]))

    return shear_flow

def place_booms(N_bottom, N_curve, R):
    """
    Distributes points equidistantly along the fuselage crossection
    :param N_bottom: number of points on the bottom plate
    :param N_curve: number of points along the curved section
    :param R: outer radius
    :return boom_co_arr: array with the coordinate pairs of each
    """

    #angles = np.linspace(np.pi/(N_curve-1), np.pi*(N_curve-2)/(N_curve-1), N_curve-2) #THIS IS OLD CODE THAT WORKED
    angles = np.linspace(0, np.pi, N_curve)
    curve_boom_co_arr = np.column_stack((R * np.cos(angles), R * np.sin(angles)))

    bottom_boom_xco = np.linspace(-R, R, N_bottom)
    bottom_boom_yco = np.zeros(N_bottom)
    bottom_boom_co_arr = np.column_stack((bottom_boom_xco, bottom_boom_yco))

    boom_co_arr = np.vstack((curve_boom_co_arr, bottom_boom_co_arr))
    return boom_co_arr


if __name__ == '__main__':
    #INPUTS
    R = 0.4 #m
    t = 0.0001 #m
    r = R-t #m
    Mx = 1943 #Nm
    My = 0 #Nm
    Sx = 0 #N
    Sy = 3678 #N
    T = 0 #Nm

    #the N_points have to be even numbers for this to work
    N_points_bottom = 100
    N_points_curve = 500

    Ixx, Iyy, Cy = calculate_inertia(R, r, t)
    boom_co_arr = place_booms(N_points_bottom, N_points_curve, R)

    stress_bending = calculate_bending_stress(Ixx, Iyy, Cy, boom_co_arr, Mx, My) #Pa
    stress_torsion = calculate_torsion_stress(boom_co_arr, T, r, t)
    stress_shear = calculate_shear_stress(Ixx, Iyy, Sx, Sy, N_points_curve, N_points_bottom, R)/t

    plot_stresses = True
    plot_shear = True
    plot_bending = False
    plot_torsion = False
    print('maximum shear stress = ', max(stress_shear/1000000), 'MPa')
    print('maximum bending stress = ', max(stress_bending/100000), 'MPa')

    if plot_stresses == True:
        if plot_shear == True:
            plt.plot([-1, +1], [Cy, Cy], label='Neutral Axis', color='black')
            plt.scatter(0, Cy, marker='x', label='Centroid', color='black')
            plt.legend(loc = 'upper right')

            plt.scatter(boom_co_arr[:, 0], boom_co_arr[:, 1], c=stress_shear, cmap='viridis', label='Shear Stress')
            cbar = plt.colorbar()
            cbar.set_label('Stress [Pa]')

            plt.xlabel('X-coordinate [m]')
            plt.xlim(-0.45, 0.45)
            plt.ylabel('Y-coordinate [m]')
            plt.title('Shear Stress')
            plt.grid()

            # Make x and y-axis use same spacing
            ax = plt.gca()
            ax.set_aspect('equal', adjustable='box')

            plt.show()


        if plot_bending == True:

            plt.scatter(0, Cy, marker='+', label='Centroid')
            plt.legend(loc = 'upper right')

            plt.scatter(boom_co_arr[:,0], boom_co_arr[:,1], c=stress_bending, cmap='viridis', label='Bending Stress')
            cbar = plt.colorbar()
            cbar.set_label('Stress [units TBD]')

            plt.xlabel('X-coordinate [m]')
            plt.ylabel('Y-coordinate [m]')
            plt.title('Bending Stress')
            plt.grid()

            #Make x and y-axis use same spacing
            ax = plt.gca()
            ax.set_aspect('equal', adjustable='box')

            plt.show()

        if plot_torsion == True:
            plt.scatter(0, Cy, marker = '+', label = 'Centroid')
            plt.legend(loc = 'upper right')

            plt.scatter(boom_co_arr[:,0], boom_co_arr[:,1], c=stress_torsion, cmap='viridis', label='Torsion Stress')
            cbar = plt.colorbar()
            cbar.set_label('Stress [Pa]')

            plt.xlabel('X-coordinate [m]')
            plt.ylabel('Y-coordinate [m]')
            plt.title('Torsion Stress')
            plt.grid()

            #Make x and y-axis use same spacing
            ax = plt.gca()
            ax.set_aspect('equal', adjustable='box')

            plt.show()



