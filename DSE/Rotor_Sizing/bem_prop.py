import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares, brentq, minimize
from DSE import const
from DSE.Rotor_Sizing.airfoil import Cl_func_clarky, Cd_func_clarky, Cp_func_clarky
from DSE.Rotor_Sizing.profile import chord_dist, twist_dist


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


def Prandtl_tip_loss(a, b, V_inf, r, omega):
    R = np.max(r)
    V0 = V_inf * (1 + a)
    V2 = omega * r * (1 - b)
    phi = np.arctan(V0/V2)
    f = N / 2 * (R - r) / (r * np.sin(phi))
    F = 2 / np.pi * np.arccos(np.exp(-f))
    return 1
    # return F

def dT_momentum(a, V_inf, r, N, F, rho=const.rho0):
    """Calculate thrust distribution using momentum theory"""
    return F * 4 * np.pi * r * V_inf**2 * (1 + a) * a / N


def dQ_momentum(a, b, V_inf, r, omega, N, F, rho=const.rho0):
    """Calculate torque distribution using momentum theory"""
    return F * 4 * np.pi * r**3 * V_inf * (1 + a) * b * omega / N
    

def solve_bem(V_inf, r, omega, chord, twist, N, rho=const.rho0):
    a, b = np.zeros_like(r), np.zeros_like(r)
    def F(x):
        a, b = x[:len(r)], x[len(r):]
        dT, dQ = dT_be(a, b, V_inf, r, omega, chord, twist, rho=rho)
        F = Prandtl_tip_loss(a, b, V_inf, r, omega)
        # exclude last element becasue force will always be zero there
        return np.append((dT - dT_momentum(a, V_inf, r, N, F, rho=rho)), (dQ - dQ_momentum(a, b, V_inf, r, omega, N, F, rho=rho)))
    x0 = np.append(a, b)

    lb, rb = -np.ones(r.shape[0]*2)+1e-5, np.ones(r.shape[0]*2)-1e-5
    # assert (F(lb)*F(rb) < 0).all(), "Must be at least one root in the interval"
    res = least_squares(F, x0, bounds=(lb, rb), ftol=1e-6, xtol=1e-6, gtol=1e-6)
    if not res.success or np.abs(F(res.x)).max() > 1e-3:
        raise ValueError(res.message)
    a, b = res.x[:len(r)], res.x[len(r):]
    # a[-1], b[-1] = -1, 1-1e-4
    dT, dQ = dT_be(a, b, V_inf, r, omega, chord, twist, rho=rho)
    return dT, dQ


def solve_pitch(T, V_inf, r, omega, chord, twist, N, rho=const.rho0):
    def F(x):
        dT, _ = solve_bem(V_inf, r, omega, chord, twist + x, N, rho=rho)
        return np.trapz(dT, r) * N - T
    pitch, res = brentq(F, 0, np.pi/2, full_output=True)
    if not res.converged:
        raise ValueError(res)
    return pitch
        

def prop_eff(T, V, r, omega, chord, twist, N, rho=const.rho0, Cl_func=Cl_func_clarky, Cd_func=Cd_func_clarky):
    pitch = solve_pitch(T, V, r, omega, chord, twist, N, rho=rho)
    dT, dQ = solve_bem(V, r, omega, chord, twist + pitch, N, rho=rho)
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
        return eff_cruise + eff_loiter

    bounds = [(0, None)] * nchord + [(0, None)] * ntwist
    bounds[-1] = (0,0)
    x0 = np.append(np.linspace(0.2, 0.1, nchord), np.linspace(45, 0, ntwist))
    ret = minimize(F, x0, bounds=bounds, options={'maxiter': 10})
    print(ret)
    
    return ret.x


if __name__ == "__main__":
    N = 3
    r = np.linspace(0.1, 0.52, 15)
    RPM = 3000
    omega = RPM * 2 * np.pi / 60
    chord = 0.1 * np.ones_like(r)
    twist = np.deg2rad(3) / (r * 2)

    V_inf = 41.2

    # a = b = np.linspace(-1+1e-4, 1-1e-4, 100)
    # a, b = np.meshgrid(a, b)
    # R = 0.4
    # c = 0.354
    # t = 0.6057
    # dT, dQ = dT_be(a, b, V_inf, R, omega, c, t)
    #     # exclude last element becasue force will always be zero there
    # dTm = dT_momentum(a, V_inf, R, N, 1)
    # dQm = dQ_momentum(a, b, V_inf, R, omega, N, 1)
    # ddT = dT - dTm
    # ddQ = dQ - dQm
    # C1 = plt.contour(a, b, ddT, colors='k')
    # plt.clabel(C1, inline=True, fontsize=10)
    # C2 = plt.contour(a, b, ddQ, colors='r')
    # plt.clabel(C2, inline=True, fontsize=10)
    # plt.show()
    # dT_be(a, b, V_inf, r, omega, chord, twist)

    # dT, dQ = solve_bem(15, r, omega, chord, twist, N)
    # plt.plot(r, dT)
    # plt.show()
    # print(np.trapz(dT, r)*N)

    # pitch = solve_pitch(200, V_inf, r, omega, chord, twist, N).linspace(-1+1e-4, 1-1e-4, 100)
    # a, b = np.meshgrid(a, b)
    # R = 0.4
    # c = 0.354
    # t = 0.6057
    # dT, dQ = dT_be(a, b, V_inf, R, omega, c, t)
    #     # exclude last element becasue force will always be zero there
    # dTm = dT_momentum(a, V_inf, R, N, 1)
    # dQm = dQ_momentum(a, b, V_inf, R, omega, N, 1)
    # ddT = dT - dTm
    # ddQ = dQ - dQm
    # C1 = plt.contour(a, b, ddT, colors='k')
    # plt.clabel(C1, inline=True, fontsize=10)
    # C2 = plt.contour(a, b, ddQ, colors='r')
    # plt.clabel(C2, inline=True, fontsize=10)
    # plt.show()
    # dT_be(a, b, V_inf, r, omega, chord, twist
    # dT, dQ = solve_bem(V_inf, r, omega, chord, twist+pitch, N)
    # plt.plot(r, dT)
    # plt.show()
    nchord = 5
    ntwist = 5
    x = prop_optim(270, 41.2, 190, 23.9, r, omega, N, nchord, ntwist)
    Pchord = x[:nchord]
    Ptwist = x[nchord:]
    chord = chord_dist(r, Pchord)
    twist = twist_dist(r, Ptwist)
    fig, ax = plt.subplots(2, 1)
    ax[0].plot(r, chord)
    ax[1].plot(r, np.rad2deg(twist))
    plt.show()
    # print(prop_eff(200, V_inf, r, omega, chord, twist, N))
    # print(np.rad2deg(pitch))