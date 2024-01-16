import numpy as np
from pathlib import Path
from DSE import const
from matplotlib import pyplot as plt
from scipy.optimize import minimize
from DSE.Rotor_Sizing.profile import chord_dist, twist_dist, P_chord_bounds, P0_chord, P0_twist, P_twist_bounds
from DSE.Rotor_Sizing.airfoil import Cl_func_clarky, Cd_func_clarky, Cp_func_clarky
import pickle

file_dir = Path(__file__).parent
Cl_func_clarky0 = lambda M, a: Cl_func_clarky(M, a)
Cd_func_clarky0 = lambda M, a: Cd_func_clarky(M, a)

def solidity(R, r, chord, N):
    return N * np.trapz(chord, r) / (np.pi * R**2)

def integrate_dT_dr(r, dT, N):
    R = r[-1]
    return np.trapz(np.where((r>=0.25*R) & (r<=0.95*R), dT, 0), r) * N

def integrate_dT_dA(r, dT, N):
    R = np.max(r)
    dT = np.mean(dT, axis=0)
    return np.trapz(np.where((r>=0.25*R) & (r<=0.95*R), dT, 0), r) * N

# def dT_vi(r, vi, rho=const.rho0):
#     return 4 * np.pi * rho * vi**2 * r

def vi_bem(r, dT_dr, N, Vc, rho=const.rho0):
    return np.sqrt(dT_dr * N / (4 * np.pi * rho * r)) * 1.15 # 1.15 is for tip and 3d effects

def dT_vi_2d(r, Vi, N, U, Vc, rho=const.rho0):
    """return dT/dA"""
    # Vi *= 1.15
    V = np.sqrt(U**2 + (Vc+Vi)**2)
    V = Vi
    return 2 * (2*np.pi*r) * rho * V * (Vi + Vc) / N / 1.15**2


def dT_dr_bem(r, omega, chord, twist, vi, Vc=0, rho=const.rho0, Cl_func=Cl_func_clarky0, Cd_func=Cd_func_clarky0, which='alpha'):
    phi = np.arctan((Vc + vi) / (omega * r))
    if which == 'alpha':
        alpha = twist
    elif which == 'theta':
        alpha = twist - phi
    else: 
        raise ValueError(f"which must be 'alpha' or 'theta', not {which}")
    V = np.sqrt((omega * r)**2 + (Vc + vi)**2)
    M = V / const.a
    dL = 0.5 * rho * V**2 * chord * Cl_func(M, np.rad2deg(alpha))
    # print(Cl_func(M, np.rad2deg(alpha)))
    dD = 0.5 * rho * V**2 * chord * Cd_func(M, np.rad2deg(alpha))
    dT = dL * np.cos(phi) - dD * np.sin(phi)
    dQ = (dL * np.sin(phi) + dD * np.cos(phi)) * r
    return dT, dQ

def dT_dA_bem(r, psi, omega, chord, theta, vi, U, Vc, rho=const.rho0, Cl_func=Cl_func_clarky0, Cd_func=Cd_func_clarky0):
    U_tot = omega * r + U * np.cos(psi)
    phi = np.arctan((Vc + vi) / U_tot)
    alpha = theta - phi
    V = np.sqrt((U_tot**2 + (Vc + vi)**2))
    M = V / const.a
    shape = M.shape
    dL = 0.5 * rho * V**2 * chord * Cl_func(M.flatten(), np.rad2deg(alpha).flatten()).reshape(shape)
    # print(Cl_func(M.flatten(), np.rad2deg(alpha).flatten()).reshape(shape)[0,:])
    dD = 0.5 * rho * V**2 * chord * Cd_func(M.flatten(), np.rad2deg(alpha).flatten()).reshape(shape)
    dT = dL * np.cos(phi) - dD * np.sin(phi)
    dQ = (dL * np.sin(phi) + dD * np.cos(phi)) * r
    return dT, dQ


def solve_dT_dr(r, omega, chord, twist, N, Cl_func=Cl_func_clarky0, Cd_func=Cd_func_clarky0, max_iter=100, tol=1e-6, which='alpha'):
    # print(sum([hash(np.sum(x)) for x in [r, omega, chord, twist, N]]))
    vi1 = np.zeros_like(r)

    for _ in range(max_iter):
        dT1, dQ = dT_dr_bem(r, omega, chord, twist, vi1, Cl_func=Cl_func, Cd_func=Cd_func, which=which)
        vi2 = vi_bem(r, dT1, N, 0)
        vi = (vi1*2 + vi2) / 3
        dT2, dQ = dT_dr_bem(r, omega, chord, twist, vi, Cl_func=Cl_func, Cd_func=Cd_func, which=which)
        if np.linalg.norm(dT1 - dT2)/len(r) < tol:
            break
        vi1 = vi
        dT1 = dT2
    else:
        raise ValueError("Failed to converge")
    return dT2, dQ

def solve_dT_dr_2d(r, omega, chord, theta, N, U, Vc, Cl_func=Cl_func_clarky0, Cd_func=Cd_func_clarky0, max_iter=500, tol=1e-8):
    psi = np.linspace(0, 2*np.pi, len(r)//2)
    r, psi = np.meshgrid(r, psi)
    vi1 = np.full_like(r, 0)
    vi2 = np.full_like(r, 30)
    # vi2 = vi_bem(r, 2*(2*np.pi*r)*dT_dA_bem(r, psi, omega, chord, theta, vi1, U, Vc), N, Vc)

    for _ in range(max_iter):
        dT1, dQ = dT_dA_bem(r, psi, omega, chord, theta, vi1, U, Vc, Cl_func=Cl_func, Cd_func=Cd_func)
        dT2, dQ = dT_dA_bem(r, psi, omega, chord, theta, vi2, U, Vc, Cl_func=Cl_func, Cd_func=Cd_func)
        dT1v = dT_vi_2d(r, vi1, N, U, Vc)
        dT2v = dT_vi_2d(r, vi2, N, U, Vc)
        y1 = dT1 - dT1v
        y2 = dT2 - dT2v
        assert (y1 * y2 < 1e-3).all()  # assert root is between r1 and r2
        vi = (vi1 + vi2) / 2
        dT, _ = dT_dA_bem(r, psi, omega, chord, theta, vi, U, Vc, Cl_func=Cl_func, Cd_func=Cd_func)
        dTv = dT_vi_2d(r, vi, N, U, Vc)
        ym = dT - dTv
        if np.linalg.norm(ym)/len(r) < tol:
            break
        vi1 = np.where(ym > 0, vi, vi1)
        vi2 = np.where(ym < 0, vi, vi2)
    else:
        raise ValueError("Failed to converge")
    return dT2, dQ, vi

def solve_c_mean(r, omega, c_mean0, chord, alpha, N, T_req=1.1*const.MTOW/4, tip_frac=0.95, max_iter=100, tol=1e-6):
    if hasattr(solve_c_mean, "c_mean0"):
        c_mean = solve_c_mean.c_mean0
    else:
        c_mean = c_mean0
    for _ in range(max_iter):
        c = np.mean(chord)
        dT, _ = solve_dT_dr(r, omega, c_mean/c*chord, alpha, N)
        T = integrate_dT_dr(r, dT, N)
        if np.abs((T - T_req)/T_req) < tol:
            break
        c_mean = c_mean * (1 + T_req/T) / 2
    else:
        raise ValueError("Failed to converge")
    solve_c_mean.c_mean0 = c_mean
    return c_mean


def figure_of_merit(x, omega, N=2, T_req=1.1*const.MTOW/4, tip_frac=0.95, max_iter=1000, tol=1e-4, optim=False, c_mean0=0.05):
    R = x[0]
    P_alpha = x[1:1+len(P0_twist)]
    P_chord = x[1+len(P0_twist):]
    r = np.linspace(0.07, R, 100)
    A = np.pi * R**2
    Vtip = omega * R
    chord = chord_dist(r, P_chord)
    alpha = twist_dist(r, P_alpha)

    c_mean = solve_c_mean(r, omega, c_mean0, chord, alpha, N, T_req=T_req, tip_frac=tip_frac, max_iter=max_iter, tol=tol)
    c = np.mean(chord)
    dT, dQ = solve_dT_dr(r, omega, c_mean/c*chord, alpha, N)
    T = integrate_dT_dr(r, dT, N)
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
    c_mean0 = 0.05
    omega = 4000 * 2 * np.pi / 60
    A = np.pi * R**2

    # print(const.a * 0.68 / 0.5 * 60 / (2 * np.pi))

    bounds = [(0.5, 0.5), *P_twist_bounds, *P_chord_bounds]
    func = lambda x: -figure_of_merit(x, omega, N=N, optim=True)
    res = minimize(func, x0, bounds=bounds, method='Nelder-Mead', options={'maxiter':10000})
    print(res)
    R = res.x[0]
    P_twist = res.x[1:1+len(P0_twist)]
    P_chord = res.x[1+len(P0_twist):]

    r = np.linspace(0.07, R, 100)
    chord = chord_dist(r, P_chord)
    alpha = twist_dist(r, P_twist)

    c_mean = solve_c_mean(r, 4000*2*np.pi/60, c_mean0, chord, alpha, N=N, T_req=1.1*const.MTOW/4, max_iter=1000)
    fac = c_mean / np.mean(chord)
    chord *= fac
    print(np.sum(chord))
    print(np.sum(alpha))
    dT, dQ = solve_dT_dr(r, omega, chord, alpha, N)
    vi_test = vi_bem(r, dT, N, 0)
    dT_test, _ = dT_dr_bem(r, omega, chord, alpha, vi=vi_test, which='alpha')
    print(sum([hash(np.sum(x)) for x in [r, omega, chord, alpha, N]]))
    plt.plot(r, dT)
    plt.show()
    vi = vi_bem(r, dT, N, 0)
    phi = np.arctan(vi / (omega * r))
    theta = phi + alpha
    plt.plot(r, np.rad2deg(phi), label='phi')
    plt.plot(r, np.rad2deg(alpha), label='alpha')
    plt.plot(r, np.rad2deg(phi+alpha), label='theta')
    plt.xlabel('r [m]')
    plt.ylabel('angle [deg]')
    plt.legend()
    plt.show()

    params = {'R': R, 'N': N, 'omega': omega, 'P_twist': P_twist, 'P_chord': P_chord, 'r': r, 'c_mean': c_mean, 'theta': theta, 'debug_dT': dT, 'chord': chord}
    pickle.dump(params, open(file_dir/'rotor_params.pkl', 'wb'))
    
    Vtip = omega * R
    print(f"Mtip = {Vtip/const.a:.3f}")
    print(f"rpm: {omega*60/(2*np.pi):.3f}")
    T = integrate_dT_dr(r, dT, N)
    print(f"T = {T:.3f} N")
    print(f'Tip and hub loss: {1 - T/(np.trapz(dT, r)*N):.3f}')
    print(f'Disk loading [N/m^2]: {T/A:.3f}')

    plt.plot(r, vi)
    plt.xlabel('r [m]')
    plt.ylabel('vi [m/s]')
    plt.ylim(bottom=0)
    plt.show()

    Q = np.trapz(dQ, r) * N
    P = Q * omega
    print(f'Torque: {Q:.3f} Nm')
    print(f'Power: {P:.3f} W')
    Ct = T / (const.rho0 * A * Vtip**2)
    Cp = P / (const.rho0 * A * Vtip**3)
    fom = Ct / Cp * np.sqrt(Ct/2)
    print(f'Ct: {Ct:.5f}') 
    print(f"figure of merit: {fom:.3f}")
    print(f"R: {R:.3f}")

    M = (omega * r) / const.a
    print(f'mean cl: {Cl_func_clarky0(M, np.rad2deg(alpha)).mean()}')


    r = np.linspace(0., R, 100)
    chord = chord_dist(r, P_chord) * fac
    print(f'Solidity: {solidity(R, r, chord, N):.5f}')
    S = np.trapz(chord, r)
    print(f'Aspect ratio: {R**2/S:.3f}')

    xp = np.linspace(0, R, len(P0_chord))
    plt.plot(r, chord)
    plt.plot(r, chord*c_mean/np.mean(chord))
    plt.scatter(xp, list(zip(*P_chord_bounds))[0])
    plt.scatter(xp, list(zip(*P_chord_bounds))[1])
    plt.scatter(xp, P_chord*R)

    plt.xlabel('r [m]')
    plt.ylabel('chord [m]')
    plt.ylim(bottom=0)
    # set aspect ratio to same scale
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()


    xp = np.linspace(0, R, len(P0_twist))
    alpha = twist_dist(r, P_twist)
    plt.plot(r, np.rad2deg(alpha))
    plt.scatter(xp, list(zip(*P_twist_bounds))[0])
    plt.scatter(xp, list(zip(*P_twist_bounds))[1])
    plt.scatter(xp, P_twist)

    plt.xlabel('r [m]')
    plt.ylabel('alpha [deg]')
    plt.ylim(bottom=0)
    plt.show()

    
    r = np.linspace(0.07, R, 100)
    chord = chord_dist(r, P_chord)
    alpha = twist_dist(r, P_twist)
    M = (omega * r) / const.a
    import matplotlib.colors as mcolors
    color_below_zero = 'pink'  # Color for values below 0
    color_above_zero = 'white'   # Color for values above 0
    epsilon = 1e-5
    cmap = mcolors.LinearSegmentedColormap.from_list('custom_colormap', 
                                                    [(0, color_below_zero), 
                                                    (0.5 - epsilon, color_below_zero), 
                                                    (0.5 + epsilon, color_above_zero), 
                                                    (1, color_above_zero)])
    fig, ax = plt.subplots()
    plt.xlim(0.1, 0.7)
    plt.ylim(-3, 15)

    plt.plot(M, np.rad2deg(alpha))

    alpha_range = np.linspace(-2, 15, 100)
    M_grid = np.linspace(0.1, 0.7, 50)
    M_grid, alpha_range = np.meshgrid(M_grid, alpha_range)
    Cpfrac = (1 + (const.gamma-1)/2*M_grid**2) / (1 + (const.gamma-1)/2)
    Cp_crit = 2 / const.gamma/M_grid**2 * (Cpfrac**(const.gamma/(const.gamma-1))-1)
    # plt.plot(M_grid, Cp_crit)
    cl = Cl_func_clarky(M_grid, alpha_range)
    cd = Cd_func_clarky(M_grid, alpha_range)
    cp = Cp_func_clarky(M_grid, alpha_range)
    cl_cd = cl / cd
    # ax.contourf(M_grid, alpha, cl, levels=100)
    # plt.scatter(*list(zip(*alpha_max)))
    CS = ax.contourf(M_grid, alpha_range, cp-Cp_crit, cmap=cmap, levels=1)
    CS = ax.contour(M_grid, alpha_range, cl/cd, colors='k', linestyles='-')
    ax.clabel(CS, inline=1, fontsize=10)
    ax.set_xlabel("M")
    ax.set_ylabel(r"$\alpha$")
    ax.set_title("Clark Y $C_L / C_D$")
    plt.show()
