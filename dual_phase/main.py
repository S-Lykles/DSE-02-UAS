import numpy as np
import matplotlib.pyplot as plt
import const as c
import power_curves.inputs as i
from power_curves.mass_frac import fuel_weight
from aero.cl_cd import dragpolar

s = 0.5
b = np.arange(1,10+s,s)
c = np.arange(1,10+c,c)

CL = np.zeros((len(b), len(c)))
CD = np.zeros((len(b), len(c)))

for i in range()
CL, CD = dragpolar(6, 3.7)
mass_f = fuel_weight(CL, CD, SFC, S, which='payload', v_cruise=const.v_cruise*1.1)


