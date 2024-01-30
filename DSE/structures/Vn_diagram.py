import numpy as np
from DSE.plot_setting import *
import matplotlib.pyplot as plt
from DSE.power_curves.rotor_tool import rotor_sizing_tool, P_profile_drag, P_induced
from DSE import const
from sympy import symbols, Eq, solve

def lift(rho, v, s, cl):
    return 0.5 * rho * v**2 * s * cl

v_range = np.linspace(0, 50, 100)
pos_lim = 3.8
neg_lim = -1.52
cl_max = 1.47
W = 160*9.81
s = 3.5
b = 6
rho0 = 1.225
vcr = 42
vne = 1.25 * vcr
amax = 15
cla = np.deg2rad(0.099) #per dgree
ucr = 15.2
ud = 7.6

def intersection(rho, s, cl, pos_lim, neg_lim):
    x = symbols('x')

    L = 0.5 * rho * x**2 * s * cl
    n1 = L / W
    n2 = -L / W

    equation = Eq(n1, pos_lim)
    equation1 = Eq(n2, neg_lim)

    intersection_points = solve(equation, x)
    intersection_points1 = solve(equation1, x)

    return intersection_points[1], intersection_points1[1]

def delta_gust(rho, s, cl, w, u, v):
    return ((0.5 * rho * cl * u * v ) / (w / s))

gcr = delta_gust(rho0, s, cl_max, W, ucr, vcr)
gd = delta_gust(rho0, s, cl_max, W, ud, vne)

gust_limits = np.array([[0,1], [vcr, 1+gcr], [vne, 1+gd], [vne, 1-gd], [vcr, 1-gcr], [0,1]])

va, va_n = intersection(rho0, s, cl_max, pos_lim, neg_lim)

L = lift(rho0, v_range, s, cl_max)
n1 = L / W
n2 = -L / W

mask = n1 <= pos_lim
mask2 = n2 >= neg_lim

vn_limits = np.array([[va_n,neg_lim], [vcr, neg_lim], [vne, 0], [vne, pos_lim], [va, pos_lim]])

plt.figure(dpi=300)
plt.title('V-n diagram', fontsize='x-large')
plt.plot(v_range[mask], n1[mask], color='green')
plt.plot(v_range[mask2], n2[mask2], color='green')
plt.plot(gust_limits[:,0], gust_limits[:, 1], linestyle='-', label='Gust Envelope')
plt.plot(vn_limits[:, 0], vn_limits[:, 1], linestyle='-', color='green', label='Vn Diagram')
#plt.plot(vcr, 0, marker='x', label='Vcr')
#plt.plot(va, 0, marker='x', color='red', label='Va')
plt.xlabel('V [m/s]')
plt.ylabel('n [-]')
plt.legend()
plt.savefig('Vn_diagram.png')
plt.show()

