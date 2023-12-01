from inputs import *
import numpy as np
from numpy import pi

def generate_Preq_ac(W, rho, S, AR, e, Cd0, eff_prop, t_start, t_end, step):
    v = np.linspace(t_start, t_end, step)

    Preq_tot = 0.5*rho*S*v**3*Cd0 + W**2 / (0.5*rho*S*v*AR*pi*e)
    Preq_shaft = Preq_tot / eff_prop
    Preq_lst = Preq_shaft

    return Preq_lst, v