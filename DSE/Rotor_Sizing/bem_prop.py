import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares, brentq, minimize, root, dual_annealing
from DSE import const
from DSE.Rotor_Sizing.airfoil import Cl_func_clarky, Cd_func_clarky, Cp_func_clarky
from DSE.Rotor_Sizing.profile import chord_dist, twist_dist
from pathlib import Path
import pandas as pd

def dT_be(a, b, V_inf, r, omega, chord, twist, rho=const.rho0, Cl_func=Cl_func_clarky, Cd_func=Cd_func_clarky):
    """Calculate thrust and torque distribution using blade element theory"""
    V0 = V_inf * (1 + a)
    V2 = omega * r * (1 - b)
    V = np.sqrt(V0**2 + V2**2)
    phi = np.arctan(V0/V2)
    alpha = twist - phi
    M = V / const.a
    dL = 0.5 * rho * V**2 * chord * Cl_func(M, np.rad2deg(alpha))
    dD = 0.5 * rho * V**2 * chord * Cd_func(M, np.rad2deg(alpha))
    dT = dL * np.cos(phi) - dD * np.sin(phi)
    dQ = (dL * np.sin(phi) + dD * np.cos(phi)) * r
    return dT, dQ


def Prandtl_tip_loss(a, b, V_inf, r, R, omega):
    V0 = V_inf * (1 + a)
    V2 = omega * r * (1 - b)
    phi = np.arctan(V0/V2)
    f = N / 2 * (R - r) / (r * np.sin(phi))
    F = 2 / np.pi * np.arccos(np.exp(-f))
    # return 1
    return F

def dT_momentum(a, V_inf, r, N, F, rho=const.rho0):
    """Calculate thrust distribution using momentum theory"""
    return F * 4 * np.pi * r * V_inf**2 * (1 + a) * a / N


def dQ_momentum(a, b, V_inf, r, omega, N, F, rho=const.rho0):
    """Calculate torque distribution using momentum theory"""
    return F * 4 * np.pi * r**3 * V_inf * (1 + a) * b * omega / N
    

def solve_bem(V_inf, r_blade, omega, chord, twist, N, rho=const.rho0):
    R = np.max(r_blade)
    a, b = np.zeros_like(r_blade), np.zeros_like(r_blade)
    def F(x, r, c, t):
        a, b = x
        dT, dQ = dT_be(a, b, V_inf, r, omega, c, t, rho=rho)
        F = Prandtl_tip_loss(a, b, V_inf, r, R, omega)
        # F = 1
        return (dT - dT_momentum(a, V_inf, r, N, F, rho=rho))**2 + (dQ - dQ_momentum(a, b, V_inf, r, omega, N, F, rho=rho))**2

    eps = 1e-3
    for i, (r, c, t) in enumerate(zip(r_blade[:-1], chord[:-1], twist[:-1])):
        if i == 0:
            x0 = np.array([0, 0])
        else:
            x0 = np.array([a[i-1], b[i-1]])
        res = minimize(F, x0=x0, args=(r, c, t), bounds=([-0.5, 1.5], [-1+eps, 0.7]))
        d = F(res.x, r, c, t)
        if res.fun > 1e-2:
            a = b = np.linspace(-1+1e-4, 1-1e-4, 100)
            a, b = np.meshgrid(a, b)
            Fprandtl = Prandtl_tip_loss(a, b, V_inf, r, R, omega)
            dT, dQ = dT_be(a, b, V_inf, r, omega, c, t)
                # exclude last element becasue force will always be zero there
            dTm = dT_momentum(a, V_inf, r, N, Fprandtl)
            dQm = dQ_momentum(a, b, V_inf, r, omega, N, Fprandtl)
            ddT = dT - dTm
            ddQ = dQ - dQm
            C1 = plt.contour(a, b, ddT, colors='k')
            plt.clabel(C1, inline=True, fontsize=10)
            C2 = plt.contour(a, b, ddQ, colors='r')
            plt.clabel(C2, inline=True, fontsize=10)
            C3 = plt.contour(a, b, ddT**2+ddQ**2, colors='b', levels=[0, 1e4, 2e4, 3e4, 4e4, 1e5, 1e6, 2e6, 3e6, 1e7])
            plt.clabel(C3, inline=True, fontsize=10)
            plt.scatter(res.x[0], res.x[1])
            plt.show()
            F(res.x, r, c, t)
            raise ValueError("BEM did not converge")
        a[i], b[i] = res.x
    a[-1], b[-1] = -1 + eps, 1 - eps
    dT, dQ = dT_be(a, b, V_inf, r_blade, omega, chord, twist, rho=rho)
    return dT, dQ, a, b


def solve_pitch(T, V_inf, r, omega, chord, twist, N, rho=const.rho0):
    def F(x):
        dT, _, _, _ = solve_bem(V_inf, r, omega, chord, twist + x, N, rho=rho)
        return np.trapz(dT, r) * N - T
    try:
        pitch, res = brentq(F, 0, np.pi/2, full_output=True)
    except ValueError:
        dT, _, _, _ = solve_bem(V_inf, r, omega, chord, twist+np.pi/4, N, rho=rho)
        plt.plot(r, dT)

        plt.show()
    if not res.converged:
        raise ValueError(res)
    return pitch
        

def prop_eff(T, V, r, omega, chord, twist, N, rho=const.rho0, Cl_func=Cl_func_clarky, Cd_func=Cd_func_clarky):
    pitch = solve_pitch(T, V, r, omega, chord, twist, N, rho=rho)
    dT, dQ, _, _ = solve_bem(V, r, omega, chord, twist + pitch, N, rho=rho)
    P = np.trapz(dQ, r) * N * omega
    return T * V / P


def prop_optim(T_cruise, V_cruise, T_loiter, V_loiter, r, omega, N, nchord, ntwist):
    def F(x):
        Pchord = x[:nchord]
        Ptwist = x[nchord:]
        chord = chord_dist(r, Pchord)
        twist = twist_dist(r, Ptwist)
        eff_cruise = prop_eff(T_cruise, V_cruise, r, omega, chord, twist, N)
        eff_loiter = prop_eff(T_loiter, V_loiter, r, omega, chord, twist, N)
        return -(eff_cruise + eff_loiter)

    bounds = [(0.1, 0.7)] * nchord + [(0,45)] * ntwist
    bounds[-1] = (0,0)
    x0 = np.append(np.linspace(0.2, 0.1, nchord), np.linspace(45, 0, ntwist))
    ret = minimize(F, x0, bounds=bounds, options={'maxiter': 10})
    print(ret)
    
    return ret.x


if __name__ == "__main__":
    import matplotlib.colors as mcolors
    color_below_zero = 'pink'  # Color for values below 0
    color_above_zero = 'white'   # Color for values above 0
    epsilon = 1e-5
    cmap = mcolors.LinearSegmentedColormap.from_list('custom_colormap', 
                                                    [(0, color_below_zero), 
                                                    (0.5 - epsilon, color_below_zero), 
                                                    (0.5 + epsilon, color_above_zero), 
                                                    (1, color_above_zero)])
    N = 2
    r = np.linspace(0.07, 0.52, 20)
    RPM = 3000
    omega = RPM * 2 * np.pi / 60

    nchord = 5
    ntwist = 5
    T_cruise = 270
    T_loiter = 190
    V_cruise = 41.2
    V_loiter = 23.9
    # x = prop_optim(270, 41.2, 190, 23.9, r, omega, N, nchord, ntwist)
    Pchord = [0.3, 0.2, 0.1, 0.09, 0.1]
    Ptwist = [40, 28, 10, 0]
    chord = chord_dist(r, Pchord)
    twist = twist_dist(r, Ptwist)
    cruise_eff = prop_eff(T_cruise, V_cruise, r, omega, chord, twist, N)
    loiter_eff = prop_eff(T_loiter, V_loiter, r, omega, chord, twist, N)
    pitch_cruise = solve_pitch(T_cruise, V_cruise, r, omega, chord, twist, N)
    pitch_loiter = solve_pitch(T_loiter, V_loiter, r, omega, chord, twist, N)
    print(cruise_eff, loiter_eff)

    for pitch, V_inf in [(pitch_cruise, V_cruise), (pitch_loiter, V_loiter)]:
        dT, dQ, a, b = solve_bem(V_inf, r, omega, chord, twist + pitch, N)
        plt.plot(r, dT)
        plt.xlabel('r [m]')
        plt.ylabel('dT [N]')
        plt.show()

        V0 = V_inf * (1 + a)
        V2 = omega * r * (1 - b)
        V = np.sqrt(V0**2 + V2**2)
        phi = np.arctan(V0/V2)
        alpha = twist + pitch - phi
        M = V / const.a
        fig, ax = plt.subplots()
        plt.xlim(0.1, 0.7)
        plt.ylim(-3, 15)
        ax.plot(M[:-1], np.rad2deg(alpha[:-1]))
        # print(np.rad2deg(alpha))
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

        plt.plot(r, dT/r)
        plt.show()

    # Pchord = x[:nchord]
    r = np.linspace(0, 0.51, 100)
    chord = chord_dist(r, Pchord)
    twist = twist_dist(r, Ptwist)
    # Ptwist = x[nchord:]
    fig, ax = plt.subplots(2, 1)
    ax[0].plot(r, chord)
    ax[0].set_ylabel('chord [m]')
    ax[0].minorticks_on()
    ax[0].grid(which='major', color='#DDDDDD', linewidth=0.8)
    ax[0].grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
    ax[0].set_ylim(bottom=0)
    ax[1].plot(r, np.rad2deg(twist))
    ax[1].set_ylabel('twist [deg]')
    ax[1].set_ylim(bottom=0)
    ax[1].set_xlabel('r [m]')
    ax[1].minorticks_on()
    ax[1].grid(which='major', color='#DDDDDD', linewidth=0.8)
    ax[1].grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
    plt.savefig(Path(__file__).parent/'propeller.pdf', bbox_inches='tight')
    plt.show()
    data = pd.DataFrame({'r': r, 'chord': chord, 'twist': np.rad2deg(twist)})
    data.to_csv('propeller.csv', index=False)
    # print(prop_eff(200, V_inf, r, omega, chord, twist, N))
    # print(np.rad2deg(pitch))