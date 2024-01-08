import numpy as np
from pathlib import Path
from DSE import const
from matplotlib import pyplot as plt

file_dir = Path(__file__).parent
data = np.genfromtxt(file_dir / r"clarky.csv", delimiter=",", skip_header=1)
alpha_clarky = data[:, 0]
cl_clarky = data[:, 1]
cd_clarky = data[:, 2]
Cl_func_clarky = lambda a: np.interp(a, alpha_clarky, cl_clarky, left=0, right=0)
Cd_func_clarky = lambda a: np.interp(a, alpha_clarky, cd_clarky)

# r = np.linspace(0, 0.5, 1000)
# omega = (5000*2*np.pi) / 60  # 5000 rpm

def chord_dist(r, a=0.03*2, b=-0.02):
    return a + b*r

def solidity(R, r, chord, N=2):
    return N * np.trapz(chord, r) / (np.pi * R**2)

def twist_dist(r, a=14, b=-22):
    return a + b*r

def dT_vi(r, vi, rho=const.rho0):
    return 4 * np.pi * rho * vi**2 * r

def dT_dr_bem(r, omega, chord, twist, vi, N=2, Vc=0, rho=const.rho0, Cl_func=Cl_func_clarky):
    alpha = twist - np.rad2deg((Vc + vi) / (omega * r))
    return 0.5 * rho * (omega * r)**2 * chord * Cl_func(alpha) * N

def vi_bem(r, dT_dr, rho=const.rho0):
    return np.sqrt(dT_dr / (4 * np.pi * rho * r))

def solve_dT_dr(r, omega, chord, twist, N=2, Cl_func=Cl_func_clarky, max_iter=100, tol=1e-6):
    """twist in deg"""
    vi = np.zeros_like(r)

    for _ in range(max_iter):
        dT1 = dT_dr_bem(r, omega, chord, twist, vi, N=N, Cl_func=Cl_func)
        vi = vi_bem(r, dT1)
        dT2 = dT_dr_bem(r, omega, chord, twist, vi, N=N, Cl_func=Cl_func)
        if np.linalg.norm(dT1 - dT2)/len(r) < tol:
            break
        dT1 = dT2

    return dT2

def dD_bem(r, dT, omega, chord, twist, N=2, Vc=0, rho=const.rho0, Cd_func=Cd_func_clarky):
    vi = vi_bem(r, dT)
    alpha = twist - np.rad2deg((Vc + vi) / (omega * r))
    dDp = 0.5 * const.rho0 * (omega*r)**2 * chord * Cd_func(alpha) * N
    return (dT * (Vc + vi) / (omega * r) + dDp)


if __name__ == "__main__":
    rho = const.rho0
    R = 0.5
    r = np.linspace(0.1, R, 1000)
    N = 2
    omega = (5000*2*np.pi) / 60  # 5000 rpm
    Vtip = omega * R
    print(f'Tip speed: {Vtip:.2f} m/s')
    A = np.pi * R**2
    chord = chord_dist(r)
    print(f'Solidity: {solidity(R, r, chord, N=N):.2f}')
    twist = twist_dist(r)
    # twist = 3 / r
    dT = solve_dT_dr(r, omega, chord, twist, N=N)
    vi = vi_bem(r, dT)
    dD = dD_bem(r, dT, omega, chord, twist, N=N)
    plt.plot(r, dT/dD)
    plt.show()
    dQ = dD * r

    plt.plot(r, dT)
    # to account for tip losses we integrate only to 0.97 R
    T = np.trapz(np.where(r<=0.97*R, dT, 0), r)
    T_ideal = np.trapz(dT, r)
    print(f"tip loss % {100*(T_ideal-T)/T_ideal:.2f}")
    print(f"thrust: {T:.2f} N")
    print(f"Disk loading: {T / (np.pi * R**2):.2f} N/m^2")
    plt.show()

    plt.plot(r, vi)
    plt.show()

    Q = np.trapz(dQ, r)
    P = Q * omega
    print(f"Torque: {Q:.2f} Nm")
    print(f"Power: {P:.2f} W")
    print(f"Power loading {P /T:.2f} W/N")
    
    CT = T / (rho * A * Vtip**2)
    CP = P / (rho * A * Vtip**3)
    M = CT / CP * np.sqrt(CT/2)
    print(f"CT: {CT:.2f}")
    print(f"figure of merit: {M:.2f}")


    # plt.plot(r, twist - np.rad2deg((vi) / (omega * r)))
    # plt.show()
        
    