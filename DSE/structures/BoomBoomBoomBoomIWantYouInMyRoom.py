import numpy as np
import matplotlib.pyplot as plt
from DSE import plot_setting


def create_boom(N, w, h):
    bottom = np.column_stack((np.linspace(-w/2, w/2, int(N/4)), -h/2*np.ones(int(N/4))))
    top = np.column_stack((np.linspace(-w/2, w/2, int(N/4)), h/2*np.ones(int(N/4))))
    left = np.column_stack((-w/2*np.ones(int(N/4)), np.linspace(-h/2, h/2, int(N/4))))
    right = np.column_stack((w/2*np.ones(int(N/4)), np.linspace(-h/2, h/2, int(N/4))))

    boom = np.vstack((bottom, right, top[::-1], left[::-1]))
    return boom

def sec_prop(w, h, t):
    Ixx = 2 * ((w*t**3/12 + w*t*(h/2)**2) + t*h**3/12)
    Iyy = 2 * ((h*t**3/12 + w*t*(w/2)**2) + t*w**3/12)
    return Ixx, Iyy

def stress_bend(Ixx, Iyy, boom_co, Mx, My):
    stress = Mx*Iyy*(boom_co[:,1])/Ixx**2 + My*Ixx*boom_co[:,0]/Iyy**2
    return stress

def shear_flow(Ixx, w, h, t, boom_co, Sy):
    """
    Calculates the absolute shear flow in a cross section. ATTENTION: these are absolute values and forces are only applied in the y-axis direction!
    :param Ixx: moment of inertia about x-axis
    :param w: width of the cross-section
    :param h: height of the cross-section
    :param t: thickness of the cross-section
    :param boom_co: array with the coordinate pairs of each element on the contour of the cross-section (starting in the lower left corner going around counter-clockwise)
    :param Sy: shear force along y-axis through the centroid
    :return stress_shear: array with the shear flow at each point in the crossection
    """
    shear_flow = np.ones(len(boom_co))
    #Bottom plate shear
    x_co = boom_co[:,0][0:int(N/4)]
    shear_flow[0:int(N/4)] = abs(Sy/Ixx * h/2 * 2*abs(x_co)*t /2)
    #Top plate shear
    shear_flow[int(N/2):int(3*N/4)] = shear_flow[0:int(N/4)]
    #Right plate shear
    y_co = boom_co[:,1][int(N/4):int(N/2)]
    shear_flow[int(N/4):int(N/2)] = abs(Sy/Ixx * (2*(h/4 + abs(y_co)/2)*(h/2 - abs(y_co))*t + w*h/2*t) / 2)
    #Left plate shear
    shear_flow[int(3*N/4):int(N)] =  shear_flow[int(N/4):int(N/2)]

    return shear_flow



if __name__ == '__main__':

    #INPUT
    N = 400         #must be divisible by 4
    w = 0.125       #m
    h = 0.125       #m
    t = 0.025       #m
    Mx = 10000      #Nm
    My = 10000      #Nm
    Sx = 00         #N
    Sy = -200000    #N


    boom = create_boom(N, w, h)
    Ixx, Iyy = sec_prop(w, h, t)
    bending = stress_bend(Ixx, Iyy, boom, Mx, My)
    shear = shear_flow(Ixx, w, h, t, boom, Sy)


    plot = True
    if plot == True:
        plt.rcParams.update(plot_setting.report_fast)
        plt.scatter(boom[:,0], boom[:,1], c = shear/t)
        plt.xlabel('x-coordinate')
        plt.ylabel('y-coordinate')
        plt.grid()
        cbar = plt.colorbar()
        cbar.set_label('Stress [Pa]')
        ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        plt.show()
