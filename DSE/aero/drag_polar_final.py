from aero_constants import *
from DSE import const

#atmospheric constants
rho = const.m2rho(500)
mu = 1.81*10**(-5)

#wing constants
C_L = np.linspace(0,1.5,1000)
AR = b**2 / S
S_ref = S
e = 1.78*(1 - 0.045 * AR **0.68) - 0.64
M = const.V_min / (const.R*const.T0*const.gamma)**0.5
xc_m_wing = 0.298 # maximum thickness location on the wing chord
Lambda_m_wing = np.arctan(np.tan(sweep_ang_25_c_rad)-(4 / AR) * ((xc_m_wing - 0.25)*((1-taper)/(1+taper))))

#empennage constants
tc_h = 0.12
tc_v = 0.12
xc_m_h = 0.3
xc_m_v = 0.3
taper_h = 1 # placeholder! change when stab&control have sizing done
taper_v = 1 # placeholder! change when stab&control have sizing done
sweep_ang_50_c_rad_h = 0 # placeholder! change when stab&control have sizing done
sweep_ang_50_c_rad_v = 22 * const.deg2rad # placeholder! change when stab&control have sizing done
Lambda_m_h = np.arctan(np.tan(sweep_ang_50_c_rad_h)-(4 / AR) * ((xc_m_h - 0.5)*((1-taper_h)/(1+taper_h))))
Lambda_m_v = np.arctan(np.tan(sweep_ang_50_c_rad_v)-(4 / AR) * ((xc_m_v - 0.5)*((1-taper_v)/(1+taper_v))))



#fuselage constants
l_fus = 2.0
d_fus = 0.8

def Reynolds(rho, V, l, mu):
    return rho*V*l/mu

def C_f_laminar(R):
    return 1.328/R**0.5

def C_f_turbulent(R, M):
    return 0.455 / ((np.log(R))**2.58*(1+0.144*M**2)**0.65)

def FF_wing(xc_m, tc, M, Lambda_m):
    FF_wing = (1 + 0.6*tc/xc_m + 100*tc**4)*(1.34 * M**0.18 * np.cos(Lambda_m)**0.28)
    return FF_wing

def FF_fuselage(f):
    return 1 + 60 / f**3 + f / 400

def FF_nacelle(f):
    return 1 + (0.35/f)

def CD0(Cf_list, FF_list, Q_list, Swet_list, Sref, CD_misc, CD_LP):
    sm = 0
    for Cf, FF, Q, Swet in Cf_list, FF_list, Q_list, Swet_list:
        sm += Cf*FF*Q*Swet
    CD0 = sm / Sref + CD_misc + CD_LP
    return CD0

def drag_polar(CL, CL_minD, AR, e, CD0):
    CDi = (CL-CL_minD)**2 / (pi * AR * e)
    return CD0 + CDi

FF_wing = FF_wing(xc_m_wing, tc, M, Lambda_m_wing)
FF_h = FF_wing(xc_m_h, tc_h, M, Lambda_m_h)
FF_v = FF_wing(xc_m_v, tc_v, M, Lambda_m_v)
FF_fuselage = FF_fuselage(l_fus/d_fus)
FF_list = (FF_wing, FF_h, FF_v, FF_fuselage)

Q_wing = 1.0
Q_fus = 1.0
Q_h = 1.08
Q_v = 1.08
Q_list = [Q_wing, Q_fus, Q_h, Q_v]

