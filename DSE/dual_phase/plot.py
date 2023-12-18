import numpy as np
import matplotlib.pyplot as plt
from DSE.power_curves.plot_fuel import plot_fuel, plot_cd0
from DSE.power_curves.plot_power import plot_power_curves
from DSE.aero.cl_cd import dragpolar_dual


if __name__ == '__main__':
    k_dl = 1.01
    n = 10
    CL_max = 1.53
    for which in ['endurance', 'payload']:
        plot_fuel(np.linspace(2,6.3,n), np.linspace(1.5,5,n), lambda b, S: dragpolar_dual(b, S,CL_start=0.01),which=which,name=None)
    # # plot_fuel(np.linspace(2,6.3,n), np.linspace(1.5,5,n), lambda b, S: dragpolar_dual(b, S,CL_start=0.01),which='payload')
    # plot_cd0(np.linspace(2,6.3,n), np.linspace(1.5,5,n), lambda b, S: S*dragpolar_dual(b, S, CL_start=0.01)[1][0])
    CD0 = dragpolar_dual(6, 3.76, CL_start=0.0)[1][0]
    # plot_power_curves([400, 500, 600], [(6, 3.76), (5, 3.5), (4, 3.2), (3, 2.7)], 1, lambda b, S: dragpolar_dual(b, S, CL_start=0.01, CL_end=CL_max), CD0, S_design=3.76, k_dl=k_dl, name='dual_phase_pc.svg')

    # for b,s in [(6, 3.763), (5, 3.5), (4, 3.2), (3, 2.7)]:
    #     CL, CD = dragpolar_dual(b, s, CL_start=0.01)
        # plt.plot(CL, CD*s, label=f'b={b}, S={s}')

    # plt.ylim(bottom=0)
    # plt.legend()
    # plt.show()
