from DSE.aero.aero_constants import *
from DSE.Locations import *
from math import *
import numpy as np

CL = CL_cruise
alpha0 = alpha_0012[14]*pi/180 # initial equilibrium condition at alpha = 2 deg
CD_0 = cd_23012[14] 
CL0 = cl_23012[14]
gamma0 = 2*pi/180
alpha = alpha_23012[16]*pi/180 # cruise at alpha = 3 deg 
CD = cd_23012[16]
CT = 0.01
V = 43
T = 288.15 - 0.0065 * 500
M0 = V/(sqrt(1.4*287.15*T))
CDM = 0 # No effect on cd due to increasing Mach number
CLM = 0 # No effect on cl due to increasing Mach number
CLM_w = 0 # No effect on cl due to increasing Mach number
CLM_h = 0 # No effect on cl due to increasing Mach number

l_fr, l_aft, l_acw,l_h,h_p,h_acw,h_h,z_h, Xcg = locations()

Xac_w = l_acw
S_w = S
St =  S_w + S_h
c = c_bar

# Longitudinal aerodynamic forces
CX = CL * alpha - CD + CT
CZ = -CL - CD * alpha

# Velocity stability derivatives
CXu = M0**2 / (1 - M0**2) *CL0 * alpha0 - 3 * CD_0 - 3 * CL0 * atan(gamma0) - M0 * CDM
CZu = -M0**2 / (1 - M0**2)  *CL0 - M0 * CDM * alpha0
print("cx",CX)
print("cz",CZ)
print("cxu",CXu)
print("czu",CZu)
Cmu = M0 * (CLM_w * (Xcg - Xac_w) * S_w / St/ c - CLM_h * (Xac_h - Xcg) * S_h / St/ c * V_h**2 * V**(-2)) + CTu_w * Zcg_w * S_w / St/ c - CTu_h * Zcg_h * S_h / St/ c * V_h**2 * V**(-2)



# Angle of Attack stability derivatives
CXalpha = CL_alpha * alpha0 + CL0 - CDalpha + CTalpha
CZalpha = -CLalpha - CDalpha * alpha0 - CD0
# Cmalpha = 

# Pitch rate stbility derivatives
CXq = 0
CZq = CL_alpha_wing * (Xcg - Xac_w) * S_w / c / St- Cl_alpha_h * (Xac_h - Xcg) * S_h * V_h**2 / c / St/ V / V
Cmq = -(CL_alpha_wing * S_w * (Xcg - Xac_w)**2 / St/ c / c + Cl_alpha_h * (Xac_h -  Xcg)**2 * S_h * V_h**2 / c / c / St/ V /V)

# Angle of attack rate stability derivatives
CXalphadott = 0
CZalphadott = - Cl_alpha_h * S_h * V_h**2 * deps_dalpha * lw / St/ V / V / c
Cmalphadott = - Cl_alpha_h * S_h * V_h**2 * deps_dalpha * lw * (Xac_h - Xcg)/ St/ V / V / c / c

# Elevator deflection derivatives
CXde = 0
# CZde = 
# Cmde = 