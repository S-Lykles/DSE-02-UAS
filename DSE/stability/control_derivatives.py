import numpy as np
from DSE import const
from DSE.Locations import locations
from DSE.aero import aero_constants
from DSE.stability.tail_sizing import horizontal_tail_sizing

de_da = horizontal_tail_sizing()[4]

V = 42 # placeholder
Vh = V
d_dt = 99999 # placeholder time step
T = 288.15 - 0.0065 * 500


##   !!! Imported Values !!!
##  !! still undifend and not placed !!
Theta_0 = 9999 # placeholder
CL_alpha_cruise = 9999 # placeholder, input from aerodynamics CL_alpaha at CL cruise.
CL_alpha_CL_0 = 9999 # placeholder, input from aerodynamics CL_alpha at CL=0
Cd0_h = 9999 # placeholder, input from aerodynamics


# Base
rho = const.rho0
m = const.total_mass
M =0.12

# Propulsion
C_t = 9999 # placeholder, input from propulsion
T1 = -9999 # placeholder, input from propulsion
T2 = -9999 # placeholder, input from propulsion
T3 = -9999 # placeholder, input from propulsion
T4 = -9999 # placeholder, input from propulsion
Tp = -9999 # placeholder, input from propulsion
Lw = -9999 # placeholder, input from propulsion

# Wing properties
S = aero_constants.S
b = aero_constants.b
c_bar = aero_constants.c_bar
Cd = aero_constants.CD_cruise[0] # placeholder, input from aerodyamics
CL_w = aero_constants.CL_cruise[0]
CL_w = aero_constants.CL_cruise[0]
CL0 = aero_constants.CL_0
sweep_ang_25_c = aero_constants.sweep_ang_25_c_rad
CL_alpha_w = aero_constants.CL_alpha_wing
Cd_alpha = aero_constants.CD_alpha_wing
Cd0_w = aero_constants.CD0_wing
Cr_w = aero_constants.c_root
taper_w = aero_constants.taper
e = aero_constants.e

# Tail properties
    # Horizontal
Sh = horizontal_tail_sizing()[0] # placeholder horizontal tail surface
AR_h = 6.8 # import from horizonal
Sh = 0.538 # import from horizontal
bh = 2.3 # import from horizontal
CL_h = aero_constants.Cl_cruise_h
Cl_alpha_h = aero_constants.Cl_alpha_h
eta = 0.95
de_da = horizontal_tail_sizing()[4]
V_h = Vh

    # Vertical
bv = 0.848956720089143
Sv = 0.379330269781324
lv = 1.2354894391032434
AR_v = 1.9
sweep_v = 22
eta_v = 0.95    # assumption
Cl_alpha_v = aero_constants.Cl_alpha_v

# Initial Calucaltions
M0 = V/(np.sqrt(1.4*287.15*T))
CDM = Cd * M0 / (1-M0**2)
beta = np.sqrt(1-M**2)


CD_alpha_w = aero_constants.CL_alpha_wing * 2 * CL_w / (np.pi * b*b/S*e)
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
Lh = 0.8*Lw
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
    CZq = -9999
    Cmq = -9999

    CYr = -9999
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
    CZ0 = -CL_w - CL_h*(Sh/S) * (V_h/V**2)
    CXalpha = - Cd_alpha
    CZalpha = - aero_constants.CL_alpha_wing - aero_constants.Cl_alpha_h
    Cmalpha = aero_constants.CL_alpha_wing * l_acw / aero_constants.c_bar- aero_constants.Cl_alpha_h * l_h / aero_constants.c_bar- CD_alpha_w * Zac / aero_constants.c_bar
    CXalphadott = 0
    CZalphadott = - Cl_alpha_h * (V_h/V)**2 * de_da * aero_constants.S_h * l_h / S / aero_constants.c_bar
    Cmalphadott = - Cl_alpha_h * (V_h/V)**2 * de_da * aero_constants.S_h * l_h**2 / S / aero_constants.c_bar/ aero_constants.c_bar

    CZq = -CL_alpha_w - CL_alpha_h*l_h*Sh/(c_bar*S)*(Vh/V)**2
    Cmq = CL_alpha_w * l_acw**2/c_bar**2 - CL_alpha_h

    CYr_v1 = 2*(Vv/V)**2*Sv*lv/(S*b)*CL_alpha_v1
    CYr_v2 = 2*(Vv/V)**2*Sv*lv/(S*b)*CL_alpha_v2
    CYr = CYr_v1+CYr_v2 + 2*CY_alpha_p*(Vp/V)**2*Sv*l_p/(S*b)

    Clr = CL_w+CL_h*Sh*bh/(S*b)*(Vh/V)**2 - zv/b*(CYr_v1+CYr_v2)
    Cnr = lv/b*CYr_v1+lv/b*CYr_v2

    CYp = -2*  8/(np.pi*3) *eta_v**2 * (bv*Sv/(b*S))* (Cl_alpha_v * AR_v) / (2 + np.sqrt(4 + (((AR_v * beta) / eta) ** 2) * (((np.tan(sweep_v * const.deg2rad)) ** 2 / beta ** 2) + 1)))
    Clp = -1* (((CL_alpha_w + Cd0_w)*Cr_w*b)/(24*S) * (1+3*taper_w)) - (( (( (Cl_alpha_h*AR_h)/(2+np.sqrt(4+(AR_h*beta/eta)**2))) + Cd0_h))/6)  # Radians
    Cnp = -lv / b * CYp - 1 / 8 * (CL_w + CL_h * Sh / S * bh / b)
    CXu = -3 * Cd0_w - 3 * CL0 * np.tan(Theta_0) - M0 * CDM # Caughey, D. A., Introduction to Aircraft Stability and Control Course Notes for AE5070, 2011
    CZu = -M0**2 / (1 - M0**2)  * (CL_w + CL_h * (Sh/S))
    CMu = (2/c_bar) * (CL_w * l_acw - CL_h * l_h - Cd0_w * Zac + C_t * h_p) * ((2 * h_p)/(V * c_bar))
    CXq = 0

A = np.array(
    [[CXu, CZu, Cmu],
     [CXalpha, CZalpha, Cmalpha],
     [CXq, CZq, Cmq],
     [CXalphadott, CZalphadott, Cmalphadott]])
print("""
          CXu, CZu, CMu,
          CXalpha, CZalpha, 0,
          CXq, CZq, Cmq,
          CXalphadott, CZalphadott, 0
""",A)



# Db = -9999
# Dc = -9999