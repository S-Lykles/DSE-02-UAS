import numpy as np
import matplotlib.pyplot as plt
from DSE import plot_setting

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
    Cy_halftube = 2*R/np.pi
    Cy_bottom = 0
    A_halftube = np.pi*R*t
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

    Ixx_halftube = np.pi*R**3*t/2 - A_halftube*Cy**2
    Ixx_bottom = A_bottom*Cy**2
    Ixx = Ixx_halftube + Ixx_bottom

    Iyy_halftube = np.pi*R**3*t/2
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

    stress_torsion = T / (np.pi*r**2*t) * np.ones(len(boom_co_arr))
    return stress_torsion

def calculate_shear_flow(Ixx, Iyy, Sx, Sy, N_curve, N_bottom, r):

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

    theta = np.linspace(0, np.pi/2, int(N_curve/2))
    s2 = np.linspace(0, r, int(N_bottom/2))
    shear12 = t * r * (Sy / Ixx * (r * np.sin(theta) - Cy * theta))
    shear2 = shear12[-1]
    shear23 =  - Sy * t / Ixx * Cy * s2 + shear2

    shear_flow = np.hstack((shear12[::-1], shear12, shear23, shear23[::-1]))

    return shear_flow

#THIS FUNCTION DOESN'T WORK PROPERLY, DON'T USE!
def better_shear_flow(Ixx, Sy, N_curve, N_bottom, R, r, Cy):
    theta = np.linspace(0, np.pi/2, int(N_curve/2))
    Qcurve = (R*np.sin(theta/2)/(theta/2)*np.cos(theta/2) - Cy)*(theta**R*t)
    qcurve = Sy/Ixx*Qcurve
    x = np.linspace(R, 0, int(N_bottom/2))
    Qbottom = -Cy*(R-x)*t + (R*np.sin(np.pi/2)/(np.pi/2)*np.cos(np.pi/4)-Cy)*np.pi/2*R*t
    qbottom = 0.5*Sy/Ixx*Qbottom
    print(Qbottom)
    shear_flow = np.hstack((qcurve[::-1], qcurve, qbottom, qbottom[::-1]))

    return shear_flow





if __name__ == '__main__':

    #INPUTS
    N_points_bottom = 100     #even number
    N_points_curve = 500    #even number
    R = 0.125 #m
    t = 0.025 #m
    r = R-t #m
    Mx = 10000 #Nm
    My = 10000 #Nm
    Sx = 0 #N
    Sy = 200000 #N
    T = 0 #Nm


    cross_sec = place_booms(N_points_bottom, N_points_curve, R)
    Ixx, Iyy, Cy = calculate_inertia(R, r, t)

    bending = calculate_bending_stress(Ixx, Iyy, Cy, cross_sec, Mx, My) #Pa
    torsion = calculate_torsion_stress(cross_sec, T, r, t)
    shear_flow = calculate_shear_flow(Ixx, Iyy, Sx, Sy, N_points_curve, N_points_bottom, R)
    shear = shear_flow/t
    shear_flow_2 = better_shear_flow(Ixx, Sy, N_points_curve, N_points_bottom, R, r, Cy)


    plot = True
    if plot == True:
        plt.rcParams.update(plot_setting.report_fast)
        plt.scatter(cross_sec[:, 0], cross_sec[:, 1], c=shear)
        cbar = plt.colorbar()
        plt.plot([-R-R/5, R+R/5], [Cy,Cy], label='neutral axis', color='black')
        plt.scatter(0, Cy, marker='x', color='black', label='centroid')
        cbar.set_label('Stress [Pa]')
        plt.xlabel('x-coordinate')
        plt.xlim(-R-R/5, R+R/5)
        plt.ylabel('y-coordinate')
        plt.grid()
        plt.legend(loc = 'upper right')

        ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        plt.show()
