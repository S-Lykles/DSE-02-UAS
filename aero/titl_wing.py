import numpy as np
from matplotlib import pyplot as plt
from math import pi, sqrt, log10

# Wetted area from statistics
def s_wet(MTOM):
    Sw = 0.262*MTOM**0.745
    return Sw

#Calculate cd0 of the fuselage
def cd0_fuselage(h,v,c,Sf):
    # ISA Calculations
    T0 = 288.15
    p0 = 101325
    T = T0 - 0.0065*h
    a = sqrt(1.4*287.05*T)                                  # speed of sound
    p = (T/T0)**(-9.81/(-0.0065*287.05))*p0
    rho = p/287.05/T

    # Reynolds and Mach number
    mu = 1.81*10**(-5)
    Re = rho*v*c/mu
    M= v/a

    # cd0 for fuselage
    cf = 0.455/((log10(Re))**2.58*(1+0.144*M**2)**0.65)     #friction coefficient from literature
    Sw = s_wet(160)                                         # Fuselage area seen from above
    cd0_fus = Sw/Sf*cf
    return cd0_fus,cf,Re, rho, T, p, M

# Calculate the cl and cd, updated with N_number of engines, evenly spaced, d_eng diameter of engine, assuming engine close to optimal slenderness ratio, fully turbulent flow over engine casing
def dragpolar_tilt_wing(b,S,h,v,c,Sf,d_eng,N_eng,Lambda,c_root):
    Sw = s_wet(160)                                         # Full confguration drag estimation of short‑to‑medium range fxed‑wing UAVs and its impact on initial sizing optimization
    cf_e = 0.01                                             # Lit. from erwin
    A = b**2/S
    e = 1.78*(1-0.045*A**0.68)-0.64                         # Preliminary Design Method and Prototype Testing of a Novel Rotors Retractable Hybrid VTOL UAV
    cd0_clean = e * Sw/S * cf_e                              # Lit. from erwin
    # cd0 of engine bodies
    cd0_eng = N_eng * 0.006                                 # Assuming fully turbulent flow around streamlined body, Hoerner drag
    # Calculating the area of the wing taken up by the engines
    S_eng = 0
    delta_eng = 1 / (N_eng/2)
    for i in range(N_eng/2):
        S_eng = S_eng + d_eng * (Lambda * c_root + i * delta_eng * (1-Lambda) * c_root)
    cd0_wing = cd0_clean + cd0_eng + ((S_eng/Sw) * cd0_clean)
    cd0_fus,cf,Re,rho, T, p, M = cd0_fuselage(h,v,c,Sf)
    cl = np.linspace(-0.4,1.5,150)
    cd_tilt_wing = cd0_wing + cl**2/(pi*A*e) + cd0_fus
    return cl,cd_tilt_wing