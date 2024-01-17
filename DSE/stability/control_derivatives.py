import numpy as np
from DSE import const
from DSE.Locations import locations
from DSE.aero import aero_constants
from DSE.stability.tail_sizing import horizontal_tail_sizing
V = 42 # placeholder
V_h = V
d_dt = 99999 # placeholder time step

rho = const.rho0
S = aero_constants.S
Sh = horizontal_tail_sizing()[0] # placeholder horizontal tail surface
b = aero_constants.b
c_bar = aero_constants.c_bar
m = const.total_mass
Cd = aero_constants.CD_cruise[0] # placeholder, input from aerodyamics
CL = aero_constants.CL_cruise[0]
CL_h = aero_constants.Cl_cruise_h
CL0 = aero_constants.CL_0
sweep_ang_25_c = aero_constants.sweep_ang_25_c_rad
CL_alpha_cruise = 9999 # placeholder, input from aerodynamics CL_alpaha at CL cruise.
CL_alpha_CL_0 = 9999 # placeholder, input from aerodynamics CL_alpha at CL=0
Theta_0 = 9999 # placeholder
T = 288.15 - 0.0065 * 500
M0 = V/(np.sqrt(1.4*287.15*T))
CDM = Cd * M0 / (1-M0**2)
CL_alpha_w = aero_constants.CL_alpha_wing
Cd_alpha = aero_constants.CD_alpha_wing
Cd0_w = aero_constants.CD0_wing
Cr_w = aero_constants.c_root
taper_w = aero_constants.taper
Cl_h = aero_constants.Cl_cruise_h
Cl_alpha_h = aero_constants.Cl_alpha_h
Cd0_h = 9999 # placeholder, input from aerodynamics
AR_h = 6.8 # import from horizonal
bh = 2.3 # import from horizontal
Cl_alpha_v = aero_constants.Cl_alpha_v
AR_v = 1.9
sweep_v = 22
eta_v = 0.95    # assumption
M =0.12 # base
beta = np.sqrt(1-M**2)
eta = 0.95

CD_alpha_w = aero_constants.CL_alpha_wing * 2 * CL / (np.pi * b*b/S*e)
Ixx = -9999 # placeholder, input from structures
Iyy = -9999 # placeholder, input from structures
Ixz = -9999 # placeholder, input from structures
# Kxz = -9999
Jxy = Ixz/(m*b**2)
Ky_2 = Iyy/(m*c_bar**2)
mu_c = m/(rho*S*c_bar)
mu_b = m/(rho*S*b)
Kx_2 = Ixx/(m*b**2)
Dc = c_bar/V * d_dt
Db = b/V * d_dt
T1 = -9999 # placeholder, input from propulsion
T2 = -9999 # placeholder, input from propulsion
T3 = -9999 # placeholder, input from propulsion
T4 = -9999 # placeholder, input from propulsion
Tp = -9999 # placeholder, input from propulsion
Lw = -9999 # placeholder, input from propulsion
Lh = -9999 # placeholder, input from propulsion
q_rad= -9999 # placeholder, input for control

l_fr, l_aft, l_acw,l_h,h_p,h_acw,h_h,z_h,X_lemac, Xcg, Zac, Zh = locations()

vtol=False
if vtol:
    CX0 = 0
    CZ0 = -1*(T1+T2+T3+T4)/(0.5*rho*S*V**2)
    CXu = 0
    CZu = -9999
    CMu = -9999
    CXalpha = 0
    CZalpha = 0
    Cmalpha = 0
    CXalphadott = 0
    CZalphadott = 0
    Cmalphadott = 0
    CZq = (T1+T2+T3+T4)*np.sin(q_rad)
    Cnr = -9999
    Cmq = -9999
    CYr = 0
    CYp = -1.87  # Ref(lit) : K.W. Booth. Effect of horizontal-tail chord on the calculated subsonic span loads and stability derivatives of isolated unswept tail assemblies in sideslip and steady roll. Technical report, NASA Memo 4-1-59 L, 1959.
    Clr = -9999
    CXu = 0
    CZu = 0
    Cmu = 0
    CXq = 0

else:
    CX0 = Tp / (0.5*rho*S*V**2) - Cd
    CZ0 = -CL_w - CL_h*(S_h/S) * (V_h/V**2)
    CXalpha = - Cd_alpha
    CZalpha = - aero_constants.CL_alpha_wing - aero_constants.Cl_alpha_h 
    Cmalpha = aero_constants.CL_alpha_wing * l_acw / aero_constants.c_bar- aero_constants.Cl_alpha_h * l_h / aero_constants.c_bar- CD_alpha_w * Zac / aero_constants.c_bar  
    CXalphadott = 0
    CZalphadott = - Cl_alpha_h * (V_h/V)**2 * deps_dalpha * aero_constants.S_h * l_h / S / aero_constants.c_bar
    Cmalphadott = - Cl_alpha_h * (V_h/V)**2 * deps_dalpha * aero_constants.S_h * l_h**2 / S / aero_constants.c_bar/ aero_constants.c_bar
    CZq = (Lw+Lh)*np.sin(q_rad)
    Cnr = -9999
    Cmq = -9999
    CYp = -2*  8/(np.pi*3) *eta_v * (bv*Sv/(b*S))* (Cl_alpha_v * AR_v) / (2 + np.sqrt(4 + (((AR_v * beta) / eta) ** 2) * (((np.tan(sweep_v * const.deg2rad)) ** 2 / beta ** 2) + 1)))
    CYr = -9999
    Clp = -1* (((CL_alpha_w + Cd0_w)*Cr_w*b)/(24*S) * (1+3*taper_w)) - (( (( (Cl_alpha_h*AR_h)/(2+np.sqrt(4+(AR_h*beta/eta)**2))) + Cd0_h))/6)  # Radians
    Clr = -9999
    Cnp = -lv / b * CYp - 1 / 8 * (CL + CL_h * Sh / S * bh / b)
    CXu = -3 * CD0 - 3 * CL0 * tan(Theta_0) - M0 * CDM # Caughey, D. A., Introduction to Aircraft Stability and Control Course Notes for AE5070, 2011
    CZu = -M0**2 / (1 - M0**2)  * (CL + CL_h * (Sh/S))
    CMu = (2/c_bar) * (CL * l_acw - Cl_h * l_h - Cd0_w * Zac + C_t * Z_m) * ((2 * Z_m)/(V * c_bar))
    CXq = 0

A = np.array(
    [[CXu, CZu, Cmu],
     [CXalpha, CZalpha, Cmalpha],
     [CXq, CZq, Cmq],
     [CXalphadott, CZalphadott, Cmalphadott]])
print("""
          CXu, CZu, Cmu,
          CXalpha, CZalpha, 0,
          CXq, CZq, Cmq,
          CXalphadott, CZalphadott, 0
""",A)



# Db = -9999
# Dc = -9999