from aero_constants import *
from DSE import const
import itertools
import matplotlib.pyplot as plt

#atmospheric constants
rho = const.m2rho(500)
mu = 1.81*10**(-5)
M = const.V_min / (const.R*const.T0*const.gamma)**0.5

#wing constants
C_L = np.linspace(0,1.5,1000)
AR = b**2 / S
S_ref = S
S_wet_wing = S*2
xc_m_wing = 0.298 # maximum thickness location on the wing chord
Lambda_m_wing = np.arctan(np.tan(sweep_ang_25_c_rad)-(4 / AR) * ((xc_m_wing - 0.25)*((1-taper)/(1+taper))))

#empennage constants
tc_h = 0.12
tc_v = 0.12
xc_m_h = 0.3
xc_m_v = 0.3
taper_h = 1 # placeholder! change when stab&control have sizing done
taper_v = 1 # placeholder! change when stab&control have sizing done
S_h = 0.539
S_wet_h = S_h * 2
c_root_h = 0.234
S_v = 0.37933 # for a single vertical tail
S_wet_v = S_v * 4
c_root_v = 0.4468
c_bar_h = c_root_h * (2/3) * ((1+taper_h+taper_h**2)/(1+taper_h))
c_bar_v = c_root_v * (2/3) * ((1+taper_v+taper_v**2)/(1+taper_v))
sweep_ang_50_c_rad_h = 0 # placeholder! change when stab&control have sizing done
sweep_ang_50_c_rad_v = 22 * const.deg2rad # placeholder! change when stab&control have sizing done
Lambda_m_h = np.arctan(np.tan(sweep_ang_50_c_rad_h)-(4 / AR) * ((xc_m_h - 0.5)*((1-taper_h)/(1+taper_h))))
Lambda_m_v = np.arctan(np.tan(sweep_ang_50_c_rad_v)-(4 / AR) * ((xc_m_v - 0.5)*((1-taper_v)/(1+taper_v))))

#winglet constants
xc_m_winglet = xc_m_wing
S_wet_winglet = S_winglet * 2
AR_winglet = b_winglet**2 / (S_winglet / 2)
Lambda_m_winglet = np.arctan(np.tan(sweep_ang_25_c_rad_winglet)-(4 / AR) * ((xc_m_winglet - 0.25)*((1-taper_winglet)/(1+taper_winglet))))
c_bar_winglet = c_root_winglet * (2/3) * ((1+taper_winglet+taper_winglet**2)/(1+taper_winglet))

#fuselage constants
l_fus = 2.0
d_fus = 0.8
S_wet_fus = pi * d_fus * 2 * l_fus

#other drag components
CD_misc_prop = 0.02
CD_misc_lg = 0.01
CD_misc = CD_misc_prop + CD_misc_lg
frac_CD_LP = 0.075

frac_lam_list = [0.5, 0, 0.5, 0.5, 0.5]

def Reynolds(rho, V, l, mu):
    return rho*V*l/mu


def Reynolds_per_component():
    R_wing = Reynolds(rho, const.V_min, c_bar, mu)
    R_fus = Reynolds(rho, const.V_min, l_fus, mu)
    R_h = Reynolds(rho, const.V_min, c_bar_h, mu)
    R_v = Reynolds(rho, const.V_min, c_bar_v, mu)
    R_winglet = Reynolds(rho, const.V_min, c_bar_winglet, mu)
    return [R_wing, R_fus, R_h, R_v, R_winglet]


def C_f_laminar(R):
    return 1.328/(R**0.5)

def C_f_turbulent(R, M):
    return 0.455 / (((np.log10(R))**2.58)*(1+0.144*M**2)**0.65)

def C_f_per_component():
    Reynold_list = Reynolds_per_component()
    C_f_list = []
    for (R, frac_lam) in itertools.zip_longest(Reynold_list, frac_lam_list):
        C_f = C_f_laminar(R) * frac_lam + C_f_turbulent(R, M) * (1-frac_lam)
        C_f_list.append(C_f)
    return C_f_list

def FF_wing(xc_m, tc, M, Lambda_m):
    FF_wing = (1 + 0.6*tc/xc_m + 100*tc**4)*(1.34 * M**0.18 * np.cos(Lambda_m)**0.28)
    return FF_wing

def FF_fuselage(f):
    return 1 + 60 / f**3 + f / 400

def FF_nacelle(f):
    return 1 + (0.35/f)

def CD0(Cf_list, FF_list, Q_list, Swet_list, Sref, CD_misc, frac_CD_LP):
    sm = 0
    CD_list = []
    for Cf, FF, Q, Swet in itertools.zip_longest(Cf_list, FF_list, Q_list, Swet_list):
        CD_list.append(Cf*FF*Q*Swet/Sref)
        sm += Cf*FF*Q*Swet
    print('CD', CD_list)
    CD0 = sm / Sref + CD_misc + ((frac_CD_LP)/(1-frac_CD_LP))* (sm / Sref + CD_misc)
    return CD0

def drag_polar(CL, AR, e, CD0):
    CDi = CL**2 / (pi * AR * e)
    return CD0 + CDi



FF_w = FF_wing(xc_m_wing, tc, M, Lambda_m_wing)
FF_h = FF_wing(xc_m_h, tc_h, M, Lambda_m_h)
FF_v = FF_wing(xc_m_v, tc_v, M, Lambda_m_v)
FF_fuselage = FF_fuselage(l_fus/d_fus)
FF_winglet = FF_wing(xc_m_winglet, tc_winglet, M, Lambda_m_winglet)


Q_wing = 1.0
Q_fus = 1.0
Q_h = 1.08
Q_v = 1.08
Q_winglet = 1.08 #PLACEHOLDER

FF_list = [FF_w, FF_fuselage, FF_h, FF_v,  FF_winglet]
Q_list = [Q_wing, Q_fus, Q_h, Q_v, Q_winglet]
S_wet_list = [S_wet_wing, S_wet_fus, S_wet_h, S_wet_v, S_wet_winglet]
C_f_list = C_f_per_component()
print(S_ref, CD_misc, frac_CD_LP)
print('C_f', C_f_list)
print('FF', FF_list)
print('Q', Q_list)
print('S_wet', S_wet_list)

CD0 = CD0(C_f_list, FF_list, Q_list, S_wet_list, S_ref, CD_misc, frac_CD_LP)
CD = drag_polar(C_L, AR, e, CD0)


