from DSE.aero.aero_constants import *
# from DSE.aero.xfoil_plots import *
from DSE.Locations import *
from math import *
import numpy as np
from DSE.const import deg2rad, rad2deg

CL = CL_cruise
# initial equilibrium condition at 2 deg before cruise angle
alpha_initial_eq = 2 * deg2rad
CD_initial_eq = 0.0077
CL_initial_eq = 1.2

gamma0 = 2*pi/180

# Values taken during cruise
alpha = alpha_23012_wing[index_cruise_wing]*pi/180 # cruise at alpha where cl.\/cd is max
CD = cd_23012_wing[index_cruise_wing]
CT = 0.01
V = 42
V_h = V # assume no interference of freestream air for horizontal tail
T = 288.15 - 0.0065 * 500
M0 = V/(sqrt(1.4*287.15*T))
CDM = CD_initial_eq * M0 / (1-M0**2) # No effect on cd due to increasing Mach number
CLM = CL_initial_eq * M0 / (1-M0**2) # No effect on cl due to increasing Mach number
CLM_w = CLM # No effect on cl due to increasing Mach number
CLM_h = CLM * S_h / S # No effect on cl due to increasing Mach number

CL_alpha = CL_alpha_wing
CD_alpha = CD_alpha_wing
CTalpha = 0 #  assumption form marilenas paper
CT0 = CD_initial_eq + CL_initial_eq * atan(gamma0) # citation 14 from marilenas control derivatives paper
CTu_w = 0
CTu_h = CTu_w * S_h / S # assumption using ratio of volumes

l_fr, l_aft, l_acw,l_h,h_p,h_acw,h_h,z_h,X_lemac, Xcg, Zac, Zh = locations()

Xac_w = l_acw
Xac_h = l_h
Zcg_w = Zac
Zcg_h = Zh
lw = l_h - l_acw
S_w = S
St =  S_w
c = c_bar

deps_dalpha = 0 # Assumption


# Longitudinal aerodynamic forces
CX = CL * alpha - CD + CT
CZ = -CL - CD * alpha

# Velocity stability derivatives
CXu = M0**2 / (1 - M0**2) * CL_initial_eq * alpha_initial_eq - 3 * CD_initial_eq - 3 * CL_initial_eq * tan(gamma0) - M0 * CDM
CZu = -M0**2 / (1 - M0**2)  * CL_initial_eq - M0 * CDM * alpha_initial_eq
Cmu = M0 * (CLM_w * (Xcg - Xac_w) * S_w / St/ c - CLM_h * (Xac_h - Xcg) * S_h / St/ c * V_h**2 * V**(-2)) + CTu_w * Zcg_w * S_w / St/ c - CTu_h * Zcg_h * S_h / St/ c * V_h**2 * V**(-2)



# Angle of Attack stability derivatives
CXalpha = CL_alpha * alpha_initial_eq + CL_initial_eq - CD_alpha + CTalpha
CZalpha = -CL_alpha - CD_alpha * alpha_initial_eq - CD_initial_eq
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

A = np.array([[CXu, CZu, Cmu],
     [CXalpha, CZalpha, 0],
     [CXq, CZq, Cmq],
     [CXalphadott, CZalphadott, 0]])
print("""
          CXu, CZu, Cmu,
          CXalpha, CZalpha, 0,
          CXq, CZq, Cmq,
          CXalphadott, CZalphadott, 0
""",A)
