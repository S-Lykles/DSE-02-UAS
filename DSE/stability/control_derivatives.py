import numpy as np
from DSE import const
from DSE.Locations import locations
from DSE.aero import aero_constants

V = -9999 # placeholder
d_dt = 99999 # placeholder time step

rho = const.rho0
S = aero_constants.S
b = aero_constants.b
c_bar = aero_constants.c_bar
m = const.total_mass
Cd = 9999 # placeholder, input from aerodyamics
CL = aero_constants.CL_cruise
CL_h = aero_constants.Cl_cruise_h
sweep_ang_25_c = aero_constants.sweep_ang_25_c_rad
CL_alpha_cruise = 9999 # placeholder, input from aerodynamics CL_alpaha at CL cruise.
CL_alpha_CL_0 = 9999 # placeholder, input from aerodynamics CL_alpha at CL=0
CL_alpha_w = aero_constants.CL_alpha_wing
Cd0_w = aero_constants.CD0_wing
Cr_w = aero_constants.c_root
taper_w = aero_constants.taper
Cl_alpha_h = aero_constants.Cl_alpha_h
Cd0_h = 9999 # placeholder, input from aerodynamics
AR_h = 6.8 # import from horizonal
Sh = 0.538 # import from horizontal
bh = 2.3 # import from horizontal
Cl_alpha_v = aero_constants.Cl_alpha_v
AR_v = 1.9
sweep_v = 22
eta_v = 0.95    # assumption
M =0.12 # base
beta = np.sqrt(1-M**2)
eta = 0.95
bv = 0.848956720089143
Sv = 0.379330269781324
lv = 1.2354894391032434


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

T1 = -9999 # placeholder, input from propulsion
T2 = -9999 # placeholder, input from propulsion
T3 = -9999 # placeholder, input from propulsion
T4 = -9999 # placeholder, input from propulsion
Tp = -9999 # placeholder, input from propulsion
Lw = -9999 # placeholder, input from propulsion
Lh = -9999 # placeholder, input from propulsion
q_rad= -9999 # placeholder, input for control

l_fr, l_aft, l_acw,l_h,h_p,h_acw,h_h,z_h,X_lemac, Xcg, Zac, Zh = locations()

vtol=True
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

else:
    CX0 = Tp / (0.5*rho*S*V**2) - Cd
    CZ0 = -CL - CL_h*(Sh/S) * (Vh/V**2)
    CXu = -9999
    CZu = -9999
    CMu = -9999
    CXalpha = - CD_alpha
    CZalpha = - aero_constants.CL_alpha_wing - aero_constants.Cl_alpha_h 
    Cmalpha = aero_constants.CL_alpha_wing * l_acw / aero_constants.c - aero_constants.Cl_alpha_h * l_h / aero_constants.c - CD_alpha_w * Zac / aero_constants.c   
    CXalphadott = 0
    CZalphadott = - CL_alphadott_w - CL_alphadott_h
    Cmalphadott = - CL_alphadott_w * l_acw / c- CL_alphadott_h * l_h / aero_constants.c - CDalphadott * Zac / c
    CZq = (Lw+Lh)*np.sin(q_rad)
    Cnr = -9999
    Cmq = -9999
    CYp = -2*  8/(np.pi*3) *eta_v * (bv*Sv/(b*S))* (Cl_alpha_v * AR_v) / (2 + np.sqrt(4 + (((AR_v * beta) / eta) ** 2) * (((np.tan(sweep_v * const.deg2rad)) ** 2 / beta ** 2) + 1)))
    CYr = -9999
    Clp = -1* (((CL_alpha_w + Cd0_w)*Cr_w*b)/(24*S) * (1+3*taper_w)) - (( (( (Cl_alpha_h*AR_h)/(2+np.sqrt(4+(AR_h*beta/eta)**2))) + Cd0_h))/6)  # Radians
    Clr = -9999
    Cnp = -lv / b *CYp - 1/8*( CL + CL_h*Sh/S*bh/b)



# Db = -9999
# Dc = -9999
