import numpy as np
import matplotlib.pyplot as plt
import const
from .compound_helicopter import dragpolar_heli
from .cl_cd import dragpolar
# from .titl_wing import dragpolar_tilt


b, S = 6, 3.763

CL_dp, CD_dp = dragpolar(b, S, 0, 1)
CL_heli, CD_heli = dragpolar_heli(b, S, 0, 1)

plt.plot(CD_dp, CL_dp, label='Dual phase')
plt.plot(CD_heli, CL_heli, label='Compound helicopter')
plt.legend()
plt.show()