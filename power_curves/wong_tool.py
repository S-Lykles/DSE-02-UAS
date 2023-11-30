from inputs import *
import numpy as np
from numpy import pi

def generate_Preq_ac(W, rho, S, AR, e, Cd0, t_start, t_end, step):
    v = np.arange(t_start, t_end, step)

    Preq_lst = []
    for i in v:
        C_L = 2*W / (rho*S*i**2)
        C_D = Cd0 + C_L**2/(pi*AR*e)
        Preq = 0.5*rho*S*i**3*Cd0+ W**2 / (0.5*rho*S*i*AR*pi*e)
        Preq_lst.append(Preq)

    return Preq_lst, v