
import numpy as np
from matplotlib import pyplot as plt
from math import pi, sqrt, log10

def s_wet(MTOM):
    Sw = 0.262*MTOM**0.745
    return Sw


def cd0_fuselage(h,v):
    # ISA Calculations
    T0 = 288.15
    p0 = 101325
    T = T0 - 0.0065*h
    a = sqrt(1.4*8.314*T)
    p = (T/T0)**(-9.81/(-0.0065*287.05))*p0
    rho = p/287.05/T

    # Reynolds and Mach number
    c = 0.6
    mu = 1.81*10**(-5)
    Re = rho*v*c/mu
    M= v/a

    # cd0 for fuselage
    cf = 0.455/((log10(Re))**2.58*(1+0.144*M**2)**0.65)
    Sw = s_wet(160)
    Sf = 2
    cd0_fus = Sw/Sf*cf
    return cd0_fus,cf,Re, rho, T, p

def dragpolar(b,S):                                  
    Sw = s_wet(160)                      # Full confguration drag estimation of short‑to‑medium range fxed‑wing UAVs and its impact on initial sizing optimization
    cf_e = 0.01                                 # Lit. from erwin
    A = b**2/S             
    e = 1.78*(1-0.045*A**0.68)-0.64             # Preliminary Design Method and Prototype Testing of a Novel Rotors Retractable Hybrid VTOL UAV 
    cd0_wing = e * Sw/S * cf_e                     # Lit. from erwin
    cd0_fus,cf,Re,rho, T, p = cd0_fuselage(500,43)
    cl = np.linspace(0,1.5,100)
    cd = (cd0_wing + cl**2/(pi*A*e)) * (1/0.34) + cd0_fus     # (1/0.34) due to presence propellors -> Aerodynamic performance of aircraft wings with stationary vertical lift propellers
    return cl,cd

def plotdragpolar(b,S):
    cl,cd = dragpolar(b,S)
    plt.plot(cd,cl)
    plt.xlabel("cd")
    plt.ylabel("cl")
    plt.title("cd vs cl")
    # plt.xlim(0,0.5)
    # plt.ylim(0,1.5)
    plt.show()

if __name__ == "__main__":
    b = 6                                       # maximum width according to REQ
    S = 3.763                                 # found from statistics
    plotdragpolar(b,S)
    h = 500
    v = 43
    cd0_fus,cf,Re, rho, T, p = cd0_fuselage(h,v)
    print("T is",T)
    print("p is",p)
    print("rho is",rho)
    print(cd0_fus)







