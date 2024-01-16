import numpy as np
from DSE import const
from DSE.Rotor_Sizing.profile import chord_dist, twist_dist
import pickle
from DSE.Rotor_Sizing.bem import solve_dT_dr, solve_dT_dr_2d, integrate_dT_dr, integrate_dT_dA, vi_bem, dT_dr_bem, dT_dA_bem
from DSE.Rotor_Sizing.airfoil import Cp_func_clarky
from pathlib import Path
from matplotlib import pyplot as plt


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
    dT, dQ, _ = solve_dT_dr_2d(r_ROTOR, OMEGA*throttle, CHORD, THETA, N, U, Vc)
    T = integrate_dT_dA(r_ROTOR, dT, N)
    Q = np.trapz(np.mean(dQ, axis=0), r_ROTOR) * N
    P = Q * OMEGA
    return T, P, Q


def plot_throttle_curve(U=0, Vc=0):
    """Plot rotor performance at different throttle settings.
    Parameters
    ----------
    U : float, optional
        Freestream velocity [m/s], by default 0
    Vc : float, optional
        Vertical climb velocity [m/s], by default 0
    """
    throttle = np.linspace(0.2, 1.1, 20)
    T, P, Q = np.array([rotor_perf(t, U, Vc) for t in throttle]).T
    fig, ax = plt.subplots(2, 1)
    o = OMEGA * 60 / 2 / np.pi
    ax[0].plot(throttle*o, T)
    ax[0].set_ylabel("Thrust [N]")
    ax[1].plot(throttle*o, P)
    ax[1].set_ylabel("Power [W]")
    # ax[2].plot(throttle*o, Q)
    # ax[2].set_ylabel("Torque [Nm]")
    # [ax.legend() for ax in ax]
    plt.show()


def plot_rotor_2d(r, omega, chord, theta, N, U, Vc, z='inflow', highlight_crit=True):
    dT, dQ, vi = solve_dT_dr_2d(r, omega, chord, theta, N, U, Vc)
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


if __name__ == "__main__":

    U = 5
    Vc = 1
    
    # OMEGA = OMEGA * 1.01
    dT, dQ = solve_dT_dr(r_ROTOR, OMEGA, CHORD, THETA, N, which='theta')
    T = integrate_dT_dr(r_ROTOR, dT, N)
    print(f"T = {T:.3f} N")
    dTA, dQa, via = solve_dT_dr_2d(r_ROTOR, OMEGA, CHORD, THETA, N, U, Vc)
    Ta = integrate_dT_dA(r_ROTOR, dTA, N)
    Q = np.trapz(dQ, r_ROTOR) * N
    print(f"Q = {Q:.3f} Nm")
    print(f"P = {Q*OMEGA:.3f} W")

    print(f"Ta = {Ta:.3f} N")

    plot_rotor_2d(r_ROTOR, OMEGA, CHORD, THETA, N, U, Vc, z='inflow', highlight_crit=True)
    plot_throttle_curve()