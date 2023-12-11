import numpy as np
import matplotlib.pyplot as plt
from power_curves.plot_fuel import plot_fuel, plot_cd0
from power_curves.plot_power import plot_power_curves
from aero.cl_cd import dragpolar_dual


if __name__ == '__main__':
    n = 10
    plot_fuel(np.linspace(2,6.3,n), np.linspace(1.5,5,n), lambda b, S: dragpolar_dual(b, S,CL_start=0.01),which='endurance')
    # plot_fuel(np.linspace(2,6.3,n), np.linspace(1.5,5,n), lambda b, S: dragpolar_dual(b, S,CL_start=0.01),which='payload')
    # plot_cd0(np.linspace(2,6.3,n), np.linspace(1.5,5,n), lambda b, S: S*dragpolar_dual(b, S, CL_start=0.01)[1][0])
    plot_power_curves([200, 300, 400, 500, 600], [(6, 3.763), (5, 3.5), (4, 3.2), (3, 2.7)], 1, lambda b, S: dragpolar_dual(b, S, CL_start=0.01), 0.01, 3.763)

    for b,s in [(6, 3.763), (5, 3.5), (4, 3.2), (3, 2.7)]:
        CL, CD = dragpolar_dual(b, s, CL_start=0.01)
        # plt.plot(CL, CD*s, label=f'b={b}, S={s}')

    # plt.ylim(bottom=0)
    # plt.legend()
    # plt.show()
