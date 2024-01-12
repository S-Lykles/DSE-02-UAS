from DSE.aero.aero_constants import *
from math import *
import numpy as np

# longitudinal aerodynamic forces
CX = CL * alpha - CD + CT
CZ = -CL - CD * alpha

# Velocity stability derivatives
CXu = M0**2 / (1 - M0**2) *CL0 * alpha0 - 3 * CD_0 - 3 * CL0 * atan(gamma0) - M0 * CDm
CZu = -M0**2 / (1 - M0**2)  *CL0 - M0 * CDm * alpha0
Cmu = M0 * (CLm_w * (Xcg - Xac_w) * S_w / S / c - CLm_h * (Xac_h - Xcg) * S_h / S / c * V_h**2 * V**(-2)) + CTu_w * Zcg_w * S_w / S / c - CTu_h * Zcg_h * S_h / S / c * V_h**2 * V**(-2)

# Angle of Attack stability derivatives
CXalpha = CL_alpha * alpha0 + CL0 - CDalpha + CTalpha
CZalpha = -CLalpha _ CDalpha * alpha0 - CD0

# Pitch rate stbility derivatives
CZq = CL_alpha_wing * (Xcg - Xac_w) * S_w / c / S - Cl_alpha_h * (Xac_h - Xcg) * S_h * V_h**2 / c / S / V / V
Cmq = -(CL_alpha_wing * S_w * (Xcg - Xac_w)**2 / S / c / c + Cl_alpha_h * (Xac_h -  Xcg)**2 * S_h * V_h**2 / c / c / S / V /V)

# Angle of attack rate stability derivatives
CZalphadott = - Cl_alpha_h * S_h 