import numpy as np
import matplotlib.pyplot as plt
from power_curves.plot_fuel import plot_fuel, plot_cd0
from power_curves.plot_power import plot_power_curves
from aero.compound_helicopter import dragpolar_comp


if __name__ == '__main__':
    n = 100
    b_design = 6
    S_design = 3.763
    CD0 = dragpolar_comp(b_design, S_design, CL_start=0.0)[1][0]
    print(CD0)
    # plot_fuel(np.linspace(2,6.3,n), np.linspace(1.5,5,n), lambda b, S: dragpolar_comp(b, S, CL_start=0.01))

    plot_cd0(np.linspace(2,6.3,n), np.linspace(1.5,5,n), lambda b, S: dragpolar_comp(b, S, CL_start=0., CL_step=1)[1][0])

    # plot_power_curves([100, 160, 220, 280], [(6, 3.763), (5, 3.5), (4, 3.2), (3, 2.7)], 1, lambda b, S: dragpolar_comp(b, S, CL_start=0.01), CD0, S_design)


    # for b,s in [(6, 3.763), (5, 3.5), (4, 3.2), (3, 2.7)]:
        # CL, CD = dragpolar_comp(b, s, CL_start=0.01)
        # plt.plot(CL, CD*s, label=f'b={b}, S={s}')

    # plt.legend()
    # plt.show()

    
