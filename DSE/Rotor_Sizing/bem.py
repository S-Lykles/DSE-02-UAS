import numpy as np
from pathlib import Path
from DSE import const
from matplotlib import pyplot as plt
from scipy.optimize import minimize
from DSE.Rotor_Sizing.profile import chord_dist, twist_dist, P_chord_bounds, P0_chord, P0_twist, P_twist_bounds
from DSE.Rotor_Sizing.airfoil import Cl_func_clarky, Cd_func_clarky, Cp_func_clarky

# file_dir = Path(__file__).parent
# data = np.genfromtxt(file_dir / r"clarky.csv", delimiter=",", skip_header=1)
# alpha_clarky = data[:, 0]
# cl_clarky = data[:, 1]
# cd_clarky = data[:, 2]
# Cl_func_clarky = lambda a: np.interp(a, alpha_clarky, cl_clarky, right=0)
# Cd_func_clarky = lambda a: np.interp(a, alpha_clarky, cd_clarky)
Cl_func_clarky0 = lambda M, a: Cl_func_clarky(M, a)
Cd_func_clarky0 = lambda M, a: Cd_func_clarky(M, a)
# Cl_func_clarky = lambda M, a: 2 * np.pi * np.deg2rad(a)
# Cd_func_clarky = lambda M, a: 0.01 * np.pi * np.deg2rad(a)**2

def solidity(R, r, chord, N=2):
    return N * np.trapz(chord, r) / (np.pi * R**2)


def dT_vi(r, vi, rho=const.rho0):
    return 4 * np.pi * rho * vi**2 * r


def vi_bem(r, dT_dr, N, rho=const.rho0):
    return np.sqrt(dT_dr * N / (4 * np.pi * rho * r))


def dT_dr_bem(r, omega, chord, twist, vi, Vc=0, rho=const.rho0, Cl_func=Cl_func_clarky, Cd_func=Cd_func_clarky):
    theta = np.arctan((Vc + vi) / (omega * r))
    alpha = twist - theta
    M = (omega * r) / const.a
    # Cpfrac = (1 + (const.gamma-1)/2*M**2) / (1 + (const.gamma-1)/2)
    # Cp_crit = 2 / const.gamma/M**2 * (Cpfrac**(const.gamma/(const.gamma-1))-1)
    # cp_min = Cp_func_clarky(M, np.rad2deg(alpha))
    # if np.any(cp_min < Cp_crit):
    #     raise ValueError("Cp_min < Cp_crit")
    dL = 0.5 * rho * (omega*r)**2 * chord * Cl_func(M, np.rad2deg(alpha))
    dD = 0.5 * rho * (omega*r)**2 * chord * Cd_func(M, np.rad2deg(alpha))
    dT = dL * np.cos(theta) - dD * np.sin(theta)
    dQ = (dL * np.sin(theta) + dD * np.cos(theta)) * r
    return dT, dQ


def solve_dT_dr(r, omega, chord, twist, N, Cl_func=Cl_func_clarky, Cd_func=Cd_func_clarky, max_iter=100, tol=1e-6):
    vi1 = np.zeros_like(r)

    for _ in range(max_iter):
        dT1, dQ = dT_dr_bem(r, omega, chord, twist, vi1, Cl_func=Cl_func, Cd_func=Cd_func)
        vi2 = vi_bem(r, dT1, N)
        vi = (vi1*2 + vi2) / 3
        dT2, dQ = dT_dr_bem(r, omega, chord, twist, vi, Cl_func=Cl_func, Cd_func=Cd_func)
        if np.linalg.norm(dT1 - dT2)/len(r) < tol:
            break
        vi1 = vi
        dT1 = dT2
    else:
        raise ValueError("Failed to converge")
    return dT2, dQ


def solve_omega(r, omega0, chord, twist, N, T_req=const.MTOW/4, tip_frac=0.97, max_iter=100, tol=1e-6):
    if hasattr(solve_omega, "omega0"):
        omega = solve_omega.omega0
    else:
        omega = omega0
    for _ in range(max_iter):
        dT, _ = solve_dT_dr(r, omega, chord, twist, N)
        T = np.trapz(np.where(r<=tip_frac*R, dT, 0), r) * N
        # T = np.trapz(dT, r) * 0.9 * N
        if np.abs((T - T_req)/T_req) < tol:
            break
        omega = omega * (1 + T_req/T) / 2
    else:
        raise ValueError("Failed to converge")
    solve_omega.omega0 = omega
    return omega


def figure_of_merit(x, N=2, T_req=const.MTOW/4, rpm0=4000, tip_frac=0.97, max_iter=1000, tol=1e-4, optim=False):
    R = x[0]
    P_twist = x[1:1+len(P0_twist)]
    P_chord = x[1+len(P0_twist):]
    r = np.linspace(0.1, R, 100)
    omega = (rpm0*2*np.pi) / 60
    A = np.pi * R**2
    Vtip = omega * R
    chord = chord_dist(r, P_chord)
    twist = twist_dist(r, P_twist)

    omega = solve_omega(r, omega, chord, twist, N, T_req=T_req, tip_frac=tip_frac, max_iter=max_iter, tol=tol)
    dT, dQ = solve_dT_dr(r, omega, chord, twist, N)
    T = np.trapz(np.where(r<=tip_frac*R, dT, 0), r) * N
    # T = np.trapz(dT, r) * 0.9 * N
    Q = np.trapz(dQ, r) * N
    P = Q * omega
    Ct = T / (const.rho0 * A * Vtip**2)
    Cp = P / (const.rho0 * A * Vtip**3)
    fom = Ct / Cp * np.sqrt(Ct/2)
    return fom


if __name__ == "__main__":
    rho = const.rho0
    R = 0.5
    N = 4
    x0 = np.array([0.5, *P0_twist, *P0_chord])

    print(const.a * 0.68 / 0.5 * 60 / (2 * np.pi))

    bounds = [(0.5, 0.5), *P_twist_bounds, *P_chord_bounds]
    func = lambda x: -figure_of_merit(x, N=N, optim=True)
    res = minimize(func, x0, bounds=bounds, method='Nelder-Mead', options={'maxiter':10000})
    print(res)
    R = res.x[0]
    P_twist = res.x[1:1+len(P0_twist)]
    P_chord = res.x[1+len(P0_twist):]

    r = np.linspace(0.1, R, 100)
    chord = chord_dist(r, P_chord)
    twist = twist_dist(r, P_twist)

    omega = solve_omega(r, 4000*2*np.pi/60, chord, twist, N=N, T_req=const.MTOW/4, tip_frac=0.97, max_iter=1000)
    dT, dQ = solve_dT_dr(r, omega, chord, twist, N)
    vi = vi_bem(r, dT, N)
    alpha = np.arctan(vi / (omega * r))

    S = np.trapz(chord, r)
    A = np.pi * R**2
    Vtip = omega * R
    print(f'Aspect ratio: {R**2/S:.3f}')
    print(f"Mtip = {Vtip/const.a:.3f}")
    print(f"rpm: {omega*60/(2*np.pi):.3f}")
    T = np.trapz(np.where(r<=0.97*R, dT, 0), r) * N
    # T = np.trapz(dT, r) * 0.9 * N
    print(f'Disk loading [N/m^2]: {T/A:.3f}')

    plt.plot(r, vi)
    plt.xlabel('r [m]')
    plt.ylabel('vi [m/s]')
    plt.ylim(bottom=0)
    plt.show()

    plt.plot(r, np.rad2deg(alpha))
    plt.xlabel('r [m]')
    plt.ylabel('alpha [deg]')
    plt.ylim(bottom=0)
    plt.show()
    Q = np.trapz(dQ, r) * N
    P = Q * omega
    Ct = T / (const.rho0 * A * Vtip**2)
    Cp = P / (const.rho0 * A * Vtip**3)
    fom = Ct / Cp * np.sqrt(Ct/2)
    print(f'Ct: {Ct:.5f}') 
    print(f"figure of merit: {fom:.3f}")
    print(f"R: {R:.3f}")

    r = np.linspace(0., R, 100)
    chord = chord_dist(r, P_chord)
    print(f'Solidity: {solidity(R, r, chord, N):.5f}')

    P_bounds_upper, P_bounds_lower = zip(*P_chord_bounds)
    P_bounds_upper = np.array(P_bounds_upper) * R  # type: ignore
    P_bounds_lower = np.array(P_bounds_lower) * R  # type: ignore
    xp = np.linspace(0, R, len(P0_chord))
    plt.plot(r, chord)
    plt.scatter(xp, P_bounds_lower)
    plt.scatter(xp, P_bounds_upper)
    plt.scatter(xp, P_chord*R)

    plt.xlabel('r [m]')
    plt.ylabel('chord [m]')
    plt.ylim(bottom=0)
    # set aspect ratio to same scale
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()

    P_bounds_upper, P_bounds_lower = zip(*P_twist_bounds)
    P_bounds_upper = np.array(P_bounds_upper)  # type: ignore
    P_bounds_lower = np.array(P_bounds_lower)  # type: ignore
    xp = np.linspace(0, R, len(P0_twist))
    twist = twist_dist(r, P_twist)
    plt.plot(r, np.rad2deg(twist))
    plt.scatter(xp, P_bounds_lower)
    plt.scatter(xp, P_bounds_upper)
    plt.scatter(xp, P_twist)

    plt.xlabel('r [m]')
    plt.ylabel('twist [deg]')
    plt.ylim(bottom=0)
    plt.show()
