
import numpy as np
from matplotlib import pyplot as plt

def dragpolar(b,S):
    cf_e = 0.01
    A = b**2/S             # Lit. from erwin
    e = 1.78*(1-0.045*A**0.68)-0.64                 #
    Sw = 0.262*160**0.745   # 
    cd0 = e * Sw/S * cf_e * 2
    cl = np.linspace(-0.2,1.2,101)
    cd = cd0 + cl**2/(3.14*A*e) 
    return cl,cd,e

def plotdragpolar(b,S):
    cl,cd = dragpolar(b,S)
    plt.plot(cd,cl)
    plt.xlabel("cd")
    plt.ylabel("cl")
    plt.show()


if __name__ == "__main__":
    plotdragpolar(6,3.763)





