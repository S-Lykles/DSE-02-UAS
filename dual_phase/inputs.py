import const
import numpy as np

"""Aircraft Parameters"""
eta = 0.75  # The propeller efficiency (Needs to be checked!)


"""Design Parameters"""
DL = 230  # Disk loading
N = 4   # Number of rotors

b = 6  # Wing span
S = 3.763  # Wing area
# N_bl = 5 Number of blades (Is this not calculated in a function??)

# T = from ISA
# R = 287
# gamma = 1.4
# c_s = np.sqrt(gamma*R*T)
# M_max = (V_max + V_tip)/c_s
# if M_max<0.95

# c =
#
# alpha_m =
# rho =
# sigma =
# omega =
# V =
# D_v =
# k =
# v_i =
# A_eq =

# psi_deg	= 20 #Arbitrarily chosen
# psi_rad	= np.radians(psi_deg)
# SFC = 0.5 #Given in ppt heli design
# C_T_sig = 0.11 #Check this value later heli design ppt 18
# k = 1.1
# eff_prop = 0.7 #PLACEHOLDER
# vc = 2
#
#
# #ENVIRONMENTAL PARAMETERS
# rho	= 1.2 #Arbitrarily chosen
# V_g	= 9.2 #Given in req's
#
# #AERODYNAMIC PARAMETERS
# S = 3.763 #input
# b = 6 #input
#
# AR = S**2 / b #PLACEHOLDER
# e = 0.77 #PLACEHOLDER
# Cd0 = 0.02345 #PLACEHOLDER
# Cl_alpha_rot = 5.73 #Given in ppt heli design
#
# # CL, CD = dragpolar(b, S, 0.2, 0.8) #From cl_cd.py
#
# #conversion factors
# HPtoWatt = 745.699872