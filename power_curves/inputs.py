import numpy as np
from numpy import pi
from aero.cl_cd import *

#CONSTANTS
g0 = 9.81

#AIRCRAFT PARAMETERS
M_gross = 160
W	= M_gross * g0 #Given in req's #Input

V_max =	50 #Minimun required cruise speed * 1.3
psi_deg	= 20 #Arbitrarily chosen
psi_rad	= psi_deg * (np.pi / 180)
SFC = 0.5 #Given in ppt heli design
C_T_sig = 0.11 #Check this value later heli design ppt 18
k = 1.1
eff_prop = 0.7 #PLACEHOLDER
vc = 2

#Rotor design parameters
DL = 230
N = 4
A_eq = 5 * 0.0929

#ENVIRONMENTAL PARAMETERS
rho	= 1.2 #Arbitrarily chosen
V_g	= 9.2 #Given in req's

#AERODYNAMIC PARAMETERS
S = 3.763 #input
b = 6 #input

AR = S**2 / b #PLACEHOLDER
e = 0.77 #PLACEHOLDER
Cd0 = 0.02345 #PLACEHOLDER
Cl_alpha_rot = 5.73 #Given in ppt heli design

CL, CD = dragpolar(b, S, 0.2, 0.8) #From cl_cd.py

#conversion factors
HPtoWatt = 745.699872
