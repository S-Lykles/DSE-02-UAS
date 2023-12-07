import numpy as np

g0 = 9.81  # [m/s^2]
rho0 = 1.225  # [kg/m^3]

deg2rad = np.pi / 180
rad2deg = 1 / deg2rad

HP2W = 745.699872
W2HP = 1 / HP2W

ft2m = 0.3048
m2ft = 1 / ft2m
lbsft2Nm = 1.3558


kts2ms = 0.514444
ms2kts = 1 / kts2ms

# SFC


# Disk loading


# Requirements
MTOW = 160 * g0  # [N]
Payload = 50 * g0  # [N]
R_cruise = 185e3  # [m]
T_cruise = 90 * 60  # [s]
T_loiter_pay = 20 * 60  # [s]
T_loiter_end = 10 * 3600  # [s]
v_cruise = R_cruise / T_cruise  # [m/s]