from .inputs import *
import numpy as np
from numpy import pi
'''
def generate_Preq_ac(W, rho, S, AR, e, Cd0, eff_prop, t_start, t_end, step):
    v = np.linspace(t_start, t_end, step)

    Preq_tot = 0.5*rho*S*v**3*Cd0 + W**2 / (0.5*rho*S*v*AR*pi*e)
    Preq_shaft = Preq_tot / eff_prop
    Preq_lst = Preq_shaft

    return Preq_lst, v
'''

def generate_Preq_ac(W, S, rho, CD, CL, eff_prop):
    v = np.sqrt(2*W / (rho*S*CL))
    D = W * (CD/ CL)
    Preq_tot = D * v
    Preq_shaft = Preq_tot / eff_prop
    Preq_lst = Preq_shaft

    return Preq_lst, v


def find_optimum_range_and_endurance_speed(Preq_array, v_array):
    # Find for endurance
    index_Preq_min = np.argmin(Preq_array)
    v_max_endurance = v_array[index_Preq_min]

    # Find for range
    slope_array = np.divide(Preq_array, v_array)
    index_slope_min = np.argmin(slope_array)
    v_max_range = v_array[index_slope_min]
    return v_max_endurance, v_max_range