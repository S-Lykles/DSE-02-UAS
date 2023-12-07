import numpy as np
from matplotlib import pyplot as plt
from math import pi, sqrt, log10
from cl_cd import s_wet, cd0_fuselage


# Calculate the cl and cd, updated with N_number of engines, evenly spaced, d_eng diameter of engine, assuming engine close to optimal slenderness ratio, fully turbulent flow over engine casing
def dragpolar_tilt_wing(b,S,h,v,c,Sf,d_eng,N_eng,Lambda,CL_start=0.,CL_end=1.2,CL_step=1000):
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
    c_root = b*Lambda/(1+Lambda)
    for i in range(N_eng//2):
        S_eng = S_eng + d_eng * (Lambda * c_root + i * delta_eng * (1-Lambda) * c_root)
    cd_prop = (S_eng/Sw) * cd0_clean
    cd0_tiltwing = cd0_clean + cd0_eng + cd_prop
    cd0_fus,cf,Re,rho, T, p, M = cd0_fuselage(h,v,c,Sf)

    cl = np.linspace(-0.4,2,150)
    cd_tilt_wing = cd0_tiltwing + cl**2/(pi*A*e) + cd0_fus
    # print("cd0 tilt wing is",cd0_tiltwing)
    # print("cd0 fus",cd0_fus)
    return cl,cd_tilt_wing, cd_prop

