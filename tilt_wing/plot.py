import numpy as np
from power_curves.plot_fuel import plot_fuel, plot_cd0
from power_curves.plot_power import plot_power_curves
from aero.titl_wing import dragpolar_tilt_wing
import const


if __name__ == '__main__':
    h = 500
    v = const.v_cruise
    c = 0.5
    S_f = 2
    d_eng = 0.5
    N_eng = 4
    Lambda = 0.45

    CL_max = 1.53

    n = 100
    for which in ['endurance', 'payload']:
        plot_fuel(np.linspace(2,6.3,n), np.linspace(1.5,5,n), lambda b, S: dragpolar_tilt_wing(b, S, h, v, c, S_f, d_eng, N_eng, Lambda, CL_start=0.01, CL_end=CL_max)[:-1], b_design=5, which=which,name='tilt_wing_mf_'+which+".svg")
    CD0 = dragpolar_tilt_wing(5, 3.76, h, v, c, S_f, d_eng, N_eng, Lambda, CL_start=0.0)[1][0]
    # plot_power_curves([400, 500, 600], [(5, 3.76), (4.5, 3.5), (4, 3.2), (3, 2.7)], 1, lambda b, S: dragpolar_tilt_wing(b, S, h, v, c, S_f, d_eng, N_eng, Lambda, CL_start=0.01, CL_end=CL_max)[:-1], CD0, 3.76, name='tilt_wing_pc.svg')
    