from DSE import const
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def compute_twist(T_C, G, l, t_theta_root, P_r_root, A_root, dt_dl, dP_r_dl, dA_dl, l_end):
    if l => l_end:
         P_r = 0
         A = 0
         t_theta = 0
    else:
        P_r = P_r_root - dP_r_dl * l
        A = A_root - dA_dl * l
        t_theta = t_theta_root - dt_dl * l

    theta =