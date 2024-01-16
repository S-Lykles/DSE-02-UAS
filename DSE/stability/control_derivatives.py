import numpy as np
from DSE import const
from DSE.Locations import locations
from DSE.aero import aero_constants

V = -9999 # placeholder
V_h = V
d_dt = 99999 # placeholder time step

rho = const.rho0
S = aero_constants.S
b = aero_constants.b
c_bar = aero_constants.c_bar
m = const.total_mass
Cd = 9999 # placeholder, input from aerodyamics
CD0 = 9999 # placeholder, input from aerodynamics, just normal Cd0 if initial equilibrium calculated for at 0 angle of attack
CL = aero_constants.CL_cruise
CL_h = aero_constants.Cl_cruise_h
CL0 = aero_constants.CL_0
sweep_ang_25_c = aero_constants.sweep_ang_25_c_rad
CL_alpha_cruise = 9999 # placeholder, input from aerodynamics CL_alpaha at CL cruise.
CL_alpha_CL_0 = 9999 # placeholder, input from aerodynamics CL_alpha at CL=0
Theta_0 = 9999 # placeholder
V = 42
T = 288.15 - 0.0065 * 500
M0 = V/(sqrt(1.4*287.15*T))
CDM = Cd * M0 / (1-M0**2)
CL_alpha_w = aero_constants.CL_alpha_wing
Cd0_w = aero_constants.CD0_wing
Cr_w = aero_constants.c_root
taper_w = aero_constants.taper
Cl_alpha_h = aero_constants.Cl_alpha_h
Cd0_h = 9999 # placeholder, input from aerodynamics
AR_h = 6.8 # import from horizonal
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
    CXu = 0
    CZu = 0
    Cmu = 0
    CXq = 0

else:
    CX0 = Tp / (0.5*rho*S*V**2) - Cd
    CZ0 = -CL - CL_h*(Sh/S) * (Vh/V**2)
    CXu = -3 * CD0 - 3 * CL0 * tan(Theta_0) - M0 * CDM # Caughey, D. A., Introduction to Aircraft Stability and Control Course Notes for AE5070, 2011
    CZu = -9999
    CMu = -9999
    CXalpha = - CD_alpha
    CZalpha = - aero_constants.CL_alpha_wing - aero_constants.Cl_alpha_h 
    Cmalpha = aero_constants.CL_alpha_wing * l_acw / aero_constants.c - aero_constants.Cl_alpha_h * l_h / aero_constants.c - CD_alpha_w * Zac / aero_constants.c   
    CXalphadott = 0
    CZalphadott = - CN_h_alpha * (V_h/V)**2 * deps_dalpha * aero_constants.S_h * l_h / S / aero_constants.c
    Cmalphadott = - CN_h_alpha * (V_h/V)**2 * deps_dalpha * aero_constants.S_h * l_h**2 / S / aero_constants.c / aero_constants.c
    CZq = (Lw+Lh)*np.sin(q_rad)
    Cnr = -9999
    Cmq = -9999
    CYp = -1.87  # Ref(lit) : K.W. Booth. Effect of horizontal-tail chord on the calculated subsonic span loads and stability derivatives of isolated unswept tail assemblies in sideslip and steady roll. Technical report, NASA Memo 4-1-59 L, 1959.
    CYr = -9999
    Clp = -1* (((CL_alpha_w + Cd0_w)*Cr_w*b)/(24*S) * (1+3*taper_w)) - (( (( (Cl_alpha_h*AR_h)/(2+np.sqrt(4+(AR_h*beta/eta)**2))) + Cd0_h))/6)
    Clr = -9999
    CXu = 0
    CZu = 0
    Cmu = 0
    CXq = 0

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



# Db = -9999
# Dc = -9999
