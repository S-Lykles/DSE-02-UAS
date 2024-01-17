import numpy as np
from pathlib import Path
from DSE import const
from matplotlib import pyplot as plt
from DSE.Rotor_Sizing.profile import chord_dist, twist_dist, P_chord_bounds, P0_chord, P0_twist, P_twist_bounds
from DSE.Rotor_Sizing.airfoil import Cl_func_clarky, Cd_func_clarky, Cp_func_clarky



def dT_vi_2d(r, Vi, N, U, Vc, rho=const.rho0, max_iter=200, tol=1e-3):
    """return dT/dA"""
    F = lambda x: 0.5*(x**2 + 2*x*Vc) - (Vi+Vc)*(np.sqrt(U**2 + (Vc+x)**2)-np.sqrt(U**2+Vc**2))
    dFdx = lambda x: 0.5*(2*x + 2*Vc) - (Vi+Vc)*(Vc+x)/np.sqrt(U**2 + (Vc+x)**2)
    x0 = Vi*2 + 1
    for _ in range(max_iter):
        x1 = x0 - F(x0)/dFdx(x0)
        if np.linalg.norm(x1-x0)/len(x0) < tol:
            break
        x0 = x1
    else:
        raise ValueError("Failed to converge")
    dV = np.sqrt(U**2 + (Vc+x0)**2)-np.sqrt(U**2+Vc**2)
    return (2*np.pi*r) * rho * (Vi+Vc) * dV / N / 1.15**2


def dT_be_2d(r, psi, omega, chord, theta, vi, U, Vc, rho=const.rho0, Cl_func=Cl_func_clarky, Cd_func=Cd_func_clarky):
    U_tot = omega * r + U * np.cos(psi)
    phi = np.arctan((Vc + vi) / U_tot)
    alpha = theta - phi
    V = np.sqrt((U_tot**2 + (Vc + vi)**2))
    M = V / const.a
    shape = M.shape
    dL = 0.5 * rho * V**2 * chord * Cl_func(M.flatten(), np.rad2deg(alpha).flatten()).reshape(shape)
    dD = 0.5 * rho * V**2 * chord * Cd_func(M.flatten(), np.rad2deg(alpha).flatten()).reshape(shape)
    dT = dL * np.cos(phi) - dD * np.sin(phi)
    dQ = (dL * np.sin(phi) + dD * np.cos(phi)) * r
    return dT, dQ


def integrate_dT_2d(r, dT, N):
    R = np.max(r)
    dT = np.mean(dT, axis=0)
    return np.trapz(np.where((r>=0.25*R) & (r<=0.95*R), dT, 0), r) * N


def solve_dT_2d(r, omega, chord, theta, N, U, Vc, Cl_func=Cl_func_clarky, Cd_func=Cd_func_clarky, max_iter=100, tol=1e-4):
    psi = np.linspace(0, 2*np.pi, len(r)//2)
    r, psi = np.meshgrid(r, psi)
    vi1 = np.full_like(r, -10)
    vi2 = np.full_like(r, 35)
    F = lambda vi: dT_be_2d(r, psi, omega, chord, theta, vi, U, Vc, Cl_func=Cl_func, Cd_func=Cd_func)[0] - dT_vi_2d(r, vi, N, U, Vc)
    y1, y2 = F(vi1), F(vi2)

    for i in range(max_iter):
        assert (y1 * y2 < 1e-3).all()  # assert root is between r1 and r2
        # dT1, dQ = dT_be_2d(r, psi, omega, chord, theta, vi1, U, Vc, Cl_func=Cl_func, Cd_func=Cd_func)
        # dT2, dQ = dT_be_2d(r, psi, omega, chord, theta, vi2, U, Vc, Cl_func=Cl_func, Cd_func=Cd_func)
        # dT1v = dT_vi_2d(r, vi1, N, U, Vc)
        # dT2v = dT_vi_2d(r, vi2, N, U, Vc)
        # y1 = dT1 - dT1v
        # y2 = dT2 - dT2v
        vi = (vi1 + vi2) / 2
        ym = F(vi)
        # dT, _ = dT_be_2d(r, psi, omega, chord, theta, vi, U, Vc, Cl_func=Cl_func, Cd_func=Cd_func)
        # dTv = dT_vi_2d(r, vi, N, U, Vc)
        # ym = dT - dTv
        if np.linalg.norm(ym)/len(r) < tol:
            break
        vi1 = np.where(ym > 0, vi, vi1)
        vi2 = np.where(ym < 0, vi, vi2)
    else:
        raise ValueError("Failed to converge")
    dT2, dQ = dT_be_2d(r, psi, omega, chord, theta, vi, U, Vc, Cl_func=Cl_func, Cd_func=Cd_func)
    return dT2, dQ, vi