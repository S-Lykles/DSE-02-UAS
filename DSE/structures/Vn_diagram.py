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
W = 180*9.81
s = 3.5
b = 6
rho0 = 1.225
vne = 50
vcr = 42

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

va, va_n = intersection(rho0, s, cl_max, pos_lim, neg_lim)

L = lift(rho0, v_range, s, cl_max)
n1 = L / W
n2 = -L / W

mask = n1 <= pos_lim
mask2 = n2 >= neg_lim

y_b = [neg_lim, neg_lim, 0, pos_lim, pos_lim]
x_b = [va_n, vcr, vne, vne, va]

plt.figure()
plt.title('V-n diagram', fontsize='x-large')
plt.plot(v_range[mask], n1[mask], color='green')
plt.plot(v_range[mask2], n2[mask2], color='green')
plt.plot(x_b, y_b, linestyle='-', color='green')
#plt.plot(vne, 0, marker='x', label='Vne')
plt.plot(vcr, 0, marker='x', label='Vcr')
plt.plot(va, 0, marker='x', color='red', label='Va')
plt.xlabel('V [m/s]')
plt.ylabel('n [-]')
plt.legend()
plt.show()
