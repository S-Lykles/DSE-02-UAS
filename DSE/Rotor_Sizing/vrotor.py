import numpy as np
from DSE import const
from DSE.Rotor_Sizing.profile import chord_dist, twist_dist
import pickle
from DSE.Rotor_Sizing.bem import solve_dT, integrate_dT, vi_bem, dT_be, dT_vi
from DSE.Rotor_Sizing.bem_2d import dT_be_2d, solve_dT_2d, integrate_dT_2d
from DSE.Rotor_Sizing.airfoil import Cp_func_clarky
from pathlib import Path
from matplotlib import pyplot as plt
from scipy.optimize import brentq
import pandas as pd
from DSE.plot_setting import report_tex, set_size, textwidth

plt.rcParams.update(report_tex)


file_dir = Path(__file__).parent
rotor_params = pickle.load(open("DSE/Rotor_Sizing/rotor_params.pkl", "rb"))
N = rotor_params["N"]
OMEGA = rotor_params["omega"]
P_twist = rotor_params["P_twist"]
debug_dT = rotor_params["debug_dT"]
r_ROTOR = rotor_params["r"]
R = np.max(r_ROTOR)
# alpha = twist_dist(r, P_twist)
THETA = rotor_params["theta"]
CHORD = rotor_params["chord"]

df = pd.DataFrame({"r [M]": r_ROTOR, "theta [deg]": np.rad2deg(THETA), "chord [M]": CHORD})
df.to_csv(file_dir/"rotor_params.csv", index=False)

def rotor_perf(throttle, U=0, Vc=0):
    """Calculate rotor performance at a given throttle setting.
    Parameters
    ----------
    throttle : float
        Throttle setting [0-1.1]
    U : float, optional
        Freestream velocity [m/s], by default 0
    Vc : float, optional
        Vertical climb velocity [m/s], by default 0

    Returns
    -------
    T : float
        Thrust [N]
    P : float
        Power [W]
    Q : float
        Torque [Nm]
    """
    assert 0 <= throttle <= 1.1, "Throttle must be between 0 and 1.1, not in %"
    dT, dQ, _ = solve_dT_2d(r_ROTOR, OMEGA*throttle, CHORD, THETA, N, U, Vc)
    T = integrate_dT_2d(r_ROTOR, dT, N)
    Q = np.trapz(np.mean(dQ, axis=0), r_ROTOR) * N
    P = Q * OMEGA
    return T, P, Q


def rotor_perf_T(T, v, alpha, max_iter=100, tol=1e-4):
    U, Vc = np.cos(alpha) * v, -np.sin(alpha) * v
    T_func = lambda throttle: rotor_perf(throttle, U, Vc)[0] - T
    x0, x1 = 0.1, 1.1
    try:
        x, ret = brentq(T_func, x0, x1, xtol=1e-3, rtol=1e-3, full_output=True)
    except ValueError:
        return np.inf
    # print(ret.function_calls)
    return rotor_perf(x, U, Vc)[1]


def plot_throttle_curve(U=0, Vc=0):
    """
    Plot rotor performance at different throttle settings.
    Parameters
    ----------
    U : float, optional
        Freestream velocity [m/s], by default 0
    Vc : float, optional
        Vertical climb velocity [m/s], by default 0
    """
    fig, ax = plt.subplots(2, 1)
    throttle = np.linspace(0.5, 1.05, 10)
    
    for u in [0, 5, 10]:
        T, P, Q = np.array([rotor_perf(t, u, Vc) for t in throttle]).T
        o = OMEGA * throttle * 60 / (2*np.pi)
        A = np.pi * R**2
        Vtip = R * o
        ax[0].plot(o, T, label=f"$U={u}$ m/s")
        ax[0].set_ylabel("Thrust [N]")
        ax[0].minorticks_on()
        ax[0].grid(which='major', color='#DDDDDD', linewidth=0.8)
        ax[0].grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
        ax[1].plot(o, P)
        ax[1].set_ylabel("Power [W]")
        ax[1].set_xlabel("Rotor speed [RPM]")
        ax[1].minorticks_on()
        ax[1].grid(which='major', color='#DDDDDD', linewidth=0.8)
        ax[1].grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
        # ax[2].plot(throttle*o, Q)
    # ax[2].set_ylabel("Torque [Nm]")
    # [ax.legend() for ax in ax]
    ax[0].legend()
    plt.tight_layout()
    plt.savefig(file_dir/"throttle_curve.pdf")
    plt.show()


def interp_power():
    """Interpolate power curve"""
    throttle = np.linspace(0.5, 1.05, 1000)
    U = np.linspace(0, 25, 10)
    Vc = np.linspace(-6, 6, 10)
    # omega = OMEGA * throttle
    T, P, Q = np.array([rotor_perf(t) for t in throttle]).T


def plot_rotor_2d(r, omega, chord, theta, N, U, Vc, z='inflow', highlight_crit=True):
    dT, dQ, vi = solve_dT_2d(r, omega, chord, theta, N, U, Vc)
    psi = np.linspace(0, 2*np.pi, len(r)//2)
    r_grid, psi_grid = np.meshgrid(r, psi)
    M = np.sqrt((r_grid * omega + U*np.cos(psi_grid))**2 + (via + Vc)**2) / const.a
    alpha = np.rad2deg(theta) - np.rad2deg(np.arctan((Vc + via) / (omega * r_grid + U*np.cos(psi_grid))))

    if z == 'inflow':
        z = vi
    elif z == 'dT_dr':
        z = dT
    elif z == 'alpha':
        z = alpha
    elif z == 'M':
        z = M
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    C = ax.contour(psi, r, z.T, levels=5, colors='k')
    ax.clabel(C, inline=True, fontsize=10)
    ax.set_rticks([])
    ax.grid(False)

    if highlight_crit:
        import matplotlib.colors as mcolors
        color_below_zero = 'pink'  # Color for values below 0
        color_above_zero = 'white'   # Color for values above 0
        epsilon = 1e-5
        cmap = mcolors.LinearSegmentedColormap.from_list('custom_colormap', 
                                                        [(0, color_below_zero), 
                                                        (0.5 - epsilon, color_below_zero), 
                                                        (0.5 + epsilon, color_above_zero), 
                                                        (1, color_above_zero)])
        shape = M.shape
        Cp = Cp_func_clarky(M.flatten(), alpha.flatten()).reshape(shape)
        
        Cpfrac = (1 + (const.gamma-1)/2*M**2) / (1 + (const.gamma-1)/2)
        Cp_crit = 2 / const.gamma/M**2 * (Cpfrac**(const.gamma/(const.gamma-1))-1)
        dCp = Cp - Cp_crit
        C2 = ax.contourf(psi, r, dCp.T, levels=1, cmap=cmap, vmin=-.5, vmax=.5)

    plt.show()


# def solve_T()


if __name__ == "__main__":

    U = 10
    Vc = 0
    
    # fig, ax = plt.subplots(2, 1, figsize=set_size(textwidth, 0.5*textwidth))
    # ax[0].plot(r_ROTOR, CHORD)
    # ax[0].set_xlabel("r [m]")
    # ax[0].set_ylabel("chord [m]")
    # # ax[0].set_ylim(0,
    # ax[0].minorticks_on()
    # ax[0].grid(which='major', color='#DDDDDD', linewidth=0.8)
    # ax[0].grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
    # ax[0].set_aspect('equal', adjustable='box')
    # ax[1].plot(r_ROTOR, np.rad2deg(THETA))
    # ax[1].set_xlabel("r [m]")
    # ax[1].set_ylabel("$\\theta$ [deg]")
    # ax[1].set_ylim(bottom=0)
    # ax[1].minorticks_on()
    # ax[1].grid(which='major', color='#DDDDDD', linewidth=0.8)
    # ax[1].grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
    
    # plt.tight_layout()
    # plt.savefig(file_dir/"rotor_params.pdf")
    # plt.show()
    
    OMEGA = OMEGA * 1
    # dT, dQ = solve_dT_dr(r_ROTOR, OMEGA, CHORD, THETA, N, which='theta')
    # T = integrate_dT_dr(r_ROTOR, dT, N)
    # print(f"T = {T:.3f} N")
    dTA, dQa, via = solve_dT_2d(r_ROTOR, OMEGA, CHORD, THETA, N, U, Vc)
    # dTA_new = dT_vi(r_ROTOR, via, N, U, Vc)
    Ta = integrate_dT_2d(r_ROTOR, dTA, N)

    Q = np.trapz(np.mean(dQa, axis=0), r_ROTOR) * N
    print(f"Q = {Q:.3f} Nm")
    print(f"P = {Q*OMEGA:.3f} W")

    print(f"Ta = {Ta:.3f} N")

    # plot_rotor_2d(r_ROTOR, OMEGA, CHORD, THETA, N, U, Vc, z='dT_dr', highlight_crit=True)
    plot_throttle_curve()
    


