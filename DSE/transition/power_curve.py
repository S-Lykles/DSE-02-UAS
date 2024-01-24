import numpy as np
import matplotlib.pyplot as plt
from DSE import const
from DSE.Rotor_Sizing.vrotor import rotor_perf_T
from DSE.plot_setting import report_tex, set_size
# from DSE.aero.Old.cl_cd import dragpolar_dual
# from DSE.Rotor_sizing.prop import prop_perf


aero_func = lambda alpha: (1.3*np.pi*alpha, 0.01 + 1e-3*(1.3*np.pi*alpha)**2)  # Cl, Cd
R_prop = 0.5
A_prop = np.pi * R_prop**2
def prop_perf(v, T):
    # if T < 0: return np.full_like(v, np.inf)
    return np.where(T>0, 0.5 * T * v * (np.sqrt(T * 2 / (A_prop * v**2 * const.rho0) + 1) + 1), np.inf)
# prop_perf = lambda v, T: 0.5 * T v

# v = np.linspace(0.1, 5, 100)
# T = 500
# plt.plot(v, prop_perf(v, T))
# plt.show()


if __name__ == "__main__":
    S= 10
    v_range = np.linspace(0.1, 25, 20)
    alpha = np.deg2rad(np.linspace(-1, 15, 30))
    alpha_opt = np.zeros_like(v_range)
    P = np.zeros((v_range.shape[0], alpha.shape[0]))
    Nr = 4

    # Tr = rotor thrust
    # Tp = propeller thrust
    for i, v in enumerate(v_range):
        CL, CD = aero_func(alpha)
        L = 0.5 * const.rho0 * v**2 * S * CL
        D = 0.5 * const.rho0 * v**2 * S * CD
        
        Fx = -D
        Fy = -const.MTOW + L
        F_req = np.sqrt(Fx**2 + Fy**2)
        angle_F_req = np.arctan2(-Fy, -Fx)
        Tr_req = F_req * np.cos(np.pi/2 + alpha - angle_F_req) / Nr
        Tp_req = F_req * np.sin(np.pi/2 + alpha - angle_F_req)

        Pp = prop_perf(v, Tp_req)
        Pr = [rotor_perf_T(Tr, v, a)  for a, Tr in zip(alpha, Tr_req)]
        plt.show()

        P[i, :] = Pp + Pr
        print(f'{(i+1)/len(v_range)*100:.1f}%')

    i = np.argmin(P, axis=1)
    Pmin = np.min(P, axis=1)
    alpha_opt = np.broadcast_to(alpha, P.shape)[range(P.shape[0]), i]

    plt.plot(v_range, alpha_opt.flatten())
    plt.show()
    plt.plot(v_range, Pmin)
    plt.show()



