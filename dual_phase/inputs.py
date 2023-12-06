import const #AIRCRAFT PARAMETERS
import numpy as np
n = 0.1

M_gross = 160
W = M_gross * const.g0 #Given in req's #Input

V_max =	50 #PLACEHOLDER, Minimun required cruise speed * 1.3
V = np.arange(0,V_max+n, n)

DL
R = np.sqrt(W / (DL * np.pi))
D = R * 2
V_tip = 140 * D ** 0.171
# T = from ISA
# R = 287
# gamma = 1.4
# c_s = np.sqrt(gamma*R*T)
# M_max = (V_max + V_tip)/c_s
# if M_max<0.95

N_bl = 5
c =

alpha_m =
rho =
sigma =
omega =
V =
D_v =
k =
v_i =
A_eq =

# psi_deg	= 20 #Arbitrarily chosen
# psi_rad	= np.radians(psi_deg)
# SFC = 0.5 #Given in ppt heli design
# C_T_sig = 0.11 #Check this value later heli design ppt 18
# k = 1.1
# eff_prop = 0.7 #PLACEHOLDER
# vc = 2
#
# #Rotor design parameters
# DL = 230
# N = 4
# A_eq = 5 * 0.0929
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