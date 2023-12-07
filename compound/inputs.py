import numpy as np
import const

#CONSTANTS
g0 = 9.81

#AIRCRAFT PARAMETERS
M_gross = 230
W	= M_gross * g0 #Given in req's #Input

v_max =	50 #Minimun required cruise speed * 1.3
psi_deg	= 20 #Arbitrarily chosen
psi_rad	= psi_deg * (np.pi / 180)
SFC = 0.5 #Given in ppt heli design
C_T_sig = 0.11 #Check this value later heli design ppt 18
k = 1.1
eff_prop = 0.7 #PLACEHOLDER
vc = 2

#Rotor design parameters
DL = 230
N = 1
A_eq = 5 * 0.0929

#ENVIRONMENTAL PARAMETERS
rho	= 1.2 #Arbitrarily chosen
V_g	= 9.2 #Given in req's

#AERODYNAMIC PARAMETERS
S = 2.8 #input
b = 5 #input
cl_max = 1.47

AR = S**2 / b #PLACEHOLDER
e = 0.77 #PLACEHOLDER
Cd0 = 0.027 #PLACEHOLDER
Cl_alpha_rot = 5.73 #Given in ppt heli design


#conversion factors
HPtoWatt = 745.699872
