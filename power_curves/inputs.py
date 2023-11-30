import numpy as np
from numpy import pi, cos

#AIRCRAFT PARAMETERS
W	= 1569.6
V_max =	60 #Minimun required cruise speed * 1.3
psi_deg	= 20 #Arbitrarily chosen
psi_rad	= psi_deg * 180 / pi
SFC = 0.5 #Given in ppt heli design
C_T_sig = 0.11 #Check this value later heli design ppt 18
k = 1.1


#ENVIRONMENTAL PARAMETERS
rho	= 1.2 #Arbitrarily chosen
V_g	= 9.2 #Given in req's


#AERODYNAMIC PARAMETERS
S = 2 #PLACEHOLDER
AR = 8 #PLACEHOLDER
e = 0.8 #PLACEHOLDER
Cd0 = 0.01 #PLACEHOLDER
Cl_alpha_rot = 5.73 #Given in ppt heli design
