import numpy as np
from matplotlib import pyplot as plt
from math import pi, sqrt, log10
import const

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

# Calculate the cl and cd 
def dragpolar_dual(b,S,h=500,v=const.v_cruise,c=None,Sf=2,CL_start=0.,CL_end=1.2,CL_step=1000):
    S_ref = 3.763
    if c == None:
        c = S/b 
        c = 0.5                            
    Sw = s_wet(160)                                         # Full confguration drag estimation of short‑to‑medium range fxed‑wing UAVs and its impact on initial sizing optimization
    cf_e = 0.01                                             # Lit. from erwin
    A = b**2/S             
    e = 1.78*(1-0.045*A**0.68)-0.64                         # Preliminary Design Method and Prototype Testing of a Novel Rotors Retractable Hybrid VTOL UAV 
    cd0_wing = e * Sw/S * cf_e                              # Lit. from erwin
    cd0_fus,cf,Re,rho, T, p, M = cd0_fuselage(h,v,c,Sf)
    cl = np.linspace(CL_start,CL_end,CL_step)
    cd_dual = 2*cd0_wing + cl**2/(pi*A*e) + cd0_fus       # (1/0.34) due to presence propellors -> Aerodynamic performance of aircraft wings with stationary vertical lift propellers
    return cl,cd_dual




