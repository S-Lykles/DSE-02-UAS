import numpy as np
from pathlib import Path
from DSE import const
from matplotlib import pyplot as plt

file_dir = Path(__file__).parent
data = np.genfromtxt(file_dir / r"clarky.csv", delimiter=",", skip_header=1)
alpha_clarky = data[:, 0]
cl_clarky = data[:, 1]
Cl_func_clarky = lambda a: np.interp(a, alpha_clarky, cl_clarky)


# r = np.linspace(0, 0.5, 1000)
# omega = (5000*2*np.pi) / 60  # 5000 rpm

def chord_dist(r, a=0.03, b=-0.06):
    return a + b*r

def twist_dist(r, a=14, b=-22):
    return a + b*r

def dT_vi(r, vi, rho=const.rho0):
    return 4 * np.pi * rho * vi**2 * r

def dT_bem(r, omega, chord, twist, vi, Vc=0, rho=const.rho0, Cl_func=Cl_func_clarky):
    alpha = twist - np.rad2deg((Vc + vi) / (omega * r))
    return 0.5 * rho * (omega * r)**2 * chord * Cl_func(alpha)

def vi_bem(r, dT, rho=const.rho0):
    return np.sqrt(dT / (4 * np.pi * rho * r))

def dT():
    r = np.linspace(0.05, 0.5, 1000)
    omega = (5000*2*np.pi) / 60  # 5000 rpm
    # chord = chord_dist(r)
    chord = 0.03
    twist = 3/r
    vi = np.zeros_like(r)

    dT1 = dT_bem(r, omega, chord, twist, vi)
    vi = vi_bem(r, dT1)
    dT2 = dT_bem(r, omega, chord, twist, vi)

    while np.linalg.norm(dT1 - dT2) > 1e-6:
        dT1 = dT2
        vi = vi_bem(r, dT1)
        dT2 = dT_bem(r, omega, chord, twist, vi)
        print(f'norm = {np.linalg.norm(dT1 - dT2)}')

    plt.plot(r, dT1)
    plt.show()

    plt.plot(r, vi)
    plt.show()
        
    


if __name__ == "__main__":
    dT()