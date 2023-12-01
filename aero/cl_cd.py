
import numpy as np
from matplotlib import pyplot as plt
from math import pi, sqrt

def dragpolar(b,S):
    MTOW = 160                                  # REQ
    Sw = 0.262*MTOW**0.745                      # Full confguration drag estimation of short‑to‑medium range fxed‑wing UAVs and its impact on initial sizing optimization
    cf_e = 0.01                                 # Lit. from erwin
    A = b**2/S             
    e = 1.78*(1-0.045*A**0.68)-0.64             # Preliminary Design Method and Prototype Testing of a Novel Rotors Retractable Hybrid VTOL UAV 
    cd0 = e * Sw/S * cf_e                       # Lit. from erwin
    cl = np.linspace(0,1.5,100)
    cd = (cd0 + cl**2/(pi*A*e)) * (1/0.34)      # (1/0.34) due to presence propellors -> Aerodynamic performance of aircraft wings with stationary vertical lift propellers
    return cl,cd

def plotdragpolar(b,S):
    cl,cd = dragpolar(b,S)
    plt.plot(cd,cl)
    plt.xlabel("cd")
    plt.ylabel("cl")
    plt.title("cd vs cl")
    plt.xlim(0,0.5)
    plt.ylim(0,1.5)
    plt.show()

if __name__ == "__main__":
    b = 6                                       # maximum width according to REQ
    S = 3.763                                   # found from statistics
    plotdragpolar(b,S)





