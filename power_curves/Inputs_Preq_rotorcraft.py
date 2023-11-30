import numpy as np
from numpy import pi, cos

DL	= 200 #Changes for each configuration
W	= 1569.6
V_max =	60 #Minimun required cruise speed * 1.3
rho	= 1.2 #Arbitrarily chosen
psi_deg	= 20 #Arbitrarily chosen
psi_rad	= psi_deg * 180 / pi
V_g	= 9.2 #Given in req's
Cl_alpha = 5.73 #Given in ppt heli design
SFC = 0.5 #Given in ppt heli design
C_T_sig = 0.11 #Check this value later heli design ppt 18
A_eq = 0.09
k = 1.1
N = 4

#rotor sizing
R               = np.sqrt(W/N/DL/pi)
V_tip           = 150*(2*R)**0.171
D_v             = 0.04*W
k_dl            = 1 + D_v/W
omega           = V_tip/R
Vne             = 1.1* V_max
mu_Vne          = 1.1*Vne/(omega*R)
Advance_ratio   = V_max / V_tip

#level flight
T_level         = W*k_dl
C_T_level       = T_level/ (rho*pi*R**2*omega**2*R**2)
sig_level       = C_T_level/C_T_sig

#turning flight
n_z             = 1 / cos(psi_rad)
T_turn          = W * k_dl * n_z
C_T_turn        = T_turn / (rho * pi * R**2 * (omega*R)**2)
sig_turn        = C_T_turn / C_T_sig


# T_gust = n_z*k_dl*
# C_T_gust = T_gust/ (Vne*pi*R**2*omega**2*R**2)

sig_max = max(sig_level, sig_turn)