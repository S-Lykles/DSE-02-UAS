import pytest
import numpy as np
import matplotlib.pyplot as plt
from DSE import const
from DSE.performance.payload_range import *
from DSE.aero.cl_cd import dragpolar_dual
from DSE import plot_setting


# @pytest.mark.skip()
@pytest.mark.parametrize("b, S", [(6, 3.763), (5, 3.763), (5, 3), (4, 2)])
@pytest.mark.parametrize("W1", [np.linspace(const.MTOW*0.6, const.MTOW*1.1, 10), const.MTOW, const.MTOW*0.9, const.MTOW*1.1])
def test_analytical_loiter(b, S, W1):
    """If auxiliary and payload power is zero, the analytical solution should be the same as the numerical solution"""

    eta = 0.75
    SFC = 320 / (1e3 * 1e3 * 3600)
    CL, CD = dragpolar_dual(b, S, CL_start=0.1, CL_end=1.2, CL_step=int(1e5))

    P_extra = 0
    if type(W1) == float:
        W1 = [W1]
    t = [10*3600]
    W = Wf_loiter(W1, t, SFC, eta, CL, CD, S, P_extra)
    W2 = W[:, -1]
    assert np.all(W2 < W1)
    idx = np.argmax(CL**3/CD**2)
    E = loiter_eq(W1, W2, S, CL[idx], CD[idx], eta, SFC)
    assert np.allclose(E, t[-1])


# @pytest.mark.skip()
@pytest.mark.parametrize("b, S", [(6, 3.763), (5, 3.763), (5, 3), (4, 2)])
@pytest.mark.parametrize("W1", [np.linspace(const.MTOW*0.6, const.MTOW*1.1, 10), const.MTOW, const.MTOW*0.9, const.MTOW*1.1])
def test_analytical_range(b, S, W1):
    """If auxiliary and payload power is zero, the analytical solution should be the same as the numerical solution"""

    eta = 0.75
    SFC = 320 / (1e3 * 1e3 * 3600)
    CL, CD = dragpolar_dual(b, S, CL_start=0.1, CL_end=1.2, CL_step=int(1e5))

    P_extra = 0
    R = np.linspace(0, const.R_cruise, 100)
    if type(W1) == float:
        W1 = [W1]
    W = Wf_range(W1, R, SFC, eta, CL, CD, S, P_extra, v_cruise=0)  # disable cruise speed requirement
    W2 = W[:, -1]
    assert np.all(W2 < W1)
    idx = np.argmax(CL/CD)
    R_breguet = range_eq(W1, W2, CL[idx], CD[idx], eta, SFC)
    assert np.allclose(R_breguet, R[-1])


def run_Wf_loiter():
    eta = 0.75
    S = 3.763
    b = 6
    SFC = 320 / (1e3 * 1e3 * 3600)
    CL, CD = dragpolar_dual(b, S, CL_start=0.1, CL_end=1.2)
    P_extra = 1000 + 800
    W1 = const.MTOW
    t, W = Wf_loiter(W1, SFC, eta, CL, CD, S, P_extra)
    # derivative of W wrt t
    dtdW = np.gradient(W, t)

    plt.plot(t/3600, dtdW)
    plt.xlabel('Time [h]')
    plt.ylabel('dW/dt [N/s]')
    plt.show()

    plt.plot(t/3600, W/const.g0)
    plt.xlabel('Time [h]')
    plt.ylabel('Weight [kg]')
    plt.show()


def run_Wf_cruise():
    eta = 0.75
    S = 3.763
    b = 6
    SFC = 320 / (1e3 * 1e3 * 3600)
    CL, CD = dragpolar_dual(b, S, CL_start=0.1, CL_end=1.2, CL_step=int(1e5))
    P_extra = 1000 + 800
    W1 = const.MTOW
    R, W = Wf_range(W1, SFC, eta, CL, CD, S, P_extra)
    # derivative of W wrt R
    dRdW = np.gradient(W, R)

    plt.plot(R/1e3, dRdW)
    plt.xlabel('Range [km]')
    plt.ylabel('dW/dR [N/m]')
    plt.show()

    plt.plot(R/1e3, W/const.g0)
    plt.xlabel('Range [km]')
    plt.ylabel('Weight [kg]')
    plt.show()


def run_Wf_payload():
    eta = 0.75
    S = 3.763
    b = 6
    SFC = 320 / (1e3 * 1e3 * 3600)
    CL, CD = dragpolar_dual(b, S, CL_start=0.1, CL_end=1.2, CL_step=int(1e4))
    P_max = 40e3
    P_aux = 1000
    
    R, W = Wf_payload_mission(CL, CD, S, eta, SFC, P_max, P_aux, drop_payload=50)

    plt.plot(R/1e3, W/const.g0)
    plt.xlabel('Range [km]')
    plt.ylabel('Weight [kg]')
    plt.show()

def run_Wf_endurance():
    eta = 0.75
    S = 3.763
    b = 6
    SFC = 320 / (1e3 * 1e3 * 3600)
    CL, CD = dragpolar_dual(b, S, CL_start=0.1, CL_end=1.2, CL_step=int(1e4))
    P_max = 40e3
    P_aux = 1000
    OEW = 100 * const.g0
    Payload = 20 * const.g0
    Wf_max = const.MTOW - OEW - Payload

    R, t, W = T_endurance_mission(Wf_max, CL, CD, S, eta, SFC, P_max, P_aux)
    ix_req = np.argmin(np.abs(R-const.R_cruise))
    iy_req = np.argmin(np.abs(t-const.T_loiter_end))
    Wf_req = W[iy_req, ix_req]
    print(f'Fuel used for endurance req: {Wf_req/const.g0:.3f}')

     # W is 2d array with R on columns and t on rows
    plt.figure(figsize=plot_setting.set_size(plot_setting.slidewidth/2, plot_setting.slideheight/2))
    plt.grid(False)
    plt.contourf(R/1e3, t/3600, W/const.g0, 200)
    c = plt.colorbar()
    c.set_label('Fuel Mass [kg]')
    plt.contour(R/1e3, t/3600, W/const.g0, levels=[Wf_max/const.g0], colors='black')
    plt.scatter([const.R_cruise/1e3], [const.T_loiter_end/3600], marker='x', color='black', label='VFS requirement')
    plt.legend()
    plt.xlabel('Range [km]')
    plt.ylabel('Time [h]')
    plt.title('Fuel used as a function of range and loiter time')
    plt.show()

def run_payload_range_diagram():
    eta = 0.75
    S = 3.763
    b = 6
    SFC = 320 / (1e3 * 1e3 * 3600)
    CL, CD = dragpolar_dual(b, S, CL_start=0.1, CL_end=1.2, CL_step=int(1e4))
    P_max = 40e3
    P_aux = 1000
    Wf_max = 20 * const.g0
    Payload_max = 70 * const.g0
    OEW = 100 * const.g0
    payload_range_diagram(OEW, Wf_max, Payload_max, CL, CD, S, eta, SFC, P_max, P_aux)



if __name__ == '__main__':
    # run_Wf_loiter()
    # run_Wf_cruise()
    # run_Wf_payload()
    plt.rcParams.update(plot_setting.report_tex)
    run_Wf_endurance()
    run_payload_range_diagram()
    # pytest.main(["-s"])
