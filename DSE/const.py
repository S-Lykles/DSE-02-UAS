import numpy as np


xlemac = 3.2
g0 = 9.81  # [m/s^2]
rho0 = 1.225  # [kg/m^3]
T0 = 288.15  # [K]
R = 287.058  # [J/kg/K]
gamma = 1.4  # [-]
P0 = 101325  # [Pa]
alpha = -0.0065  # [K/m]
a = np.sqrt(gamma * R * T0)  # [m/s]

V_min = 42 # [m/s]
def m2rho(m):
    """Return density of air at given altitude in kg/m^3"""
    return rho0 * (1 - alpha * m / T0)**(g0 / alpha / R - 1)


HPtoWatt = 745.6998720065

deg2rad = np.pi / 180
rad2deg = 1 / deg2rad

HP2W = 745.699872
W2HP = 1 / HP2W

ft2m = 0.3048
m2ft = 1 / ft2m
lbsft2Nm = 1.3558


kts2ms = 0.514444
ms2kts = 1 / kts2ms

WH2J = 3600


# SFC
J_KG2WH_KG = 1 / 3600
WH_KG2J_KG = 1 / J_KG2WH_KG

# Disk loading

# Requirements
total_mass = 160
MTOW = total_mass * g0  # [N]
Payload = 50 * g0  # [N]
R_cruise = 185e3  # [m]
T_cruise = 3600 * 1.25  # [s]
T_loiter_pay = 20 * 60  # [s]
T_loiter_end = 10 * 3600  # [s]
v_cruise = R_cruise / T_cruise  # [m/s]
h_loiter = 500  # [m]

P_pay_pay = 400  # [W] Power required for payload during payload phase
P_pay_end = 800  # [W] Power required for payload during endurance phase


P_aux = 1000  # [W] Power required for auxiliaries

