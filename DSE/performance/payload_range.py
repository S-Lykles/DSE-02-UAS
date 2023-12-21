import numpy as np
from scipy.integrate import solve_ivp
from numpy.typing import NDArray
import matplotlib.pyplot as plt
from DSE import plot_setting
from DSE import const


def range_eq(W1, W2, CL: float, CD: float, eta: float, SFC: float):
    """Brequet range equation"""
    SFC = SFC * const.g0 # [kg/W/s] -> [N/W/s]
    return eta * CL / (SFC * CD) * np.log(W1/W2)


def loiter_eq(W1, W2, S, CL, CD, eta, SFC, rho=const.m2rho(500)):
    """Endurance equation similar to Breguet range equation"""
    SFC = SFC * const.g0 # [kg/W/s] -> [N/W/s]
    return eta * CL**(3/2) * np.sqrt(2*rho*S) / (SFC * CD) * (1/np.sqrt(W2) - 1/np.sqrt(W1))


# Take of phase
def Wf_take_off(Pmax, SFC, t=4*60):
    """
    Calculates fuel weight during take-off phase in N
    """
    return Pmax * SFC * t * const.g0


def Wf_climb(H: float, W:float, SFC: float):
    """
    Calculates fuel weight during climb phase in N
    """
    return W * SFC * H * const.g0


def Wf_range(W1, R, SFC: float, eta, CL, CD, S, P_extra, h=500, v_cruise=const.v_cruise):
    """
    Parameters
    ----------
    W1 : array_like, shape (n,)
        Initial weight [N]
    R : array_like, shape (m,)
        Range [m], should be an equally spaced array

    Returns
    -------
    W : numpy.ndarray shape (n, m)
        Weight [N]
    """
    rho = const.m2rho(h)
    CL = CL[..., np.newaxis]
    CD = CD[..., np.newaxis]

    def dw_dr(R, W):
        v = np.sqrt(W * 2 / (rho * S * CL))
        P = W*CD/CL*v + P_extra
        dwdr = (-P * SFC / eta * const.g0) / v
        return np.max(np.where(v > v_cruise, dwdr, -np.inf), axis=0)

    sol = solve_ivp(dw_dr, (0, np.max(R)), W1, rtol=1e-8, atol=1e-8, dense_output=True)
    return sol.sol(R)


def Wf_loiter(W1, t, SFC: float, eta, CL, CD, S, P_extra):
    """
    Parameters
    ----------
    W1 : array_like, shape (n,)
        Initial weight [N]
    t : array_like, shape (m,)
        Time [s], should be an equally spaced array

    Returns
    -------
    W : numpy.ndarray shape (n, m)
        Weight [N]
    """
    rho = const.m2rho(const.h_loiter)
    CL = CL[..., np.newaxis]
    CD = CD[..., np.newaxis]

    def dw_dt(t, W):
        v = np.sqrt(W * 2 / (rho * S * CL))
        P = W * CD / CL * v + P_extra
        dwdt = -P * SFC / eta * const.g0
        return np.max(dwdt, axis=0)

    sol = solve_ivp(dw_dt, (0, np.max(t)), W1, rtol=1e-6, atol=1e-6, dense_output=True)
    return sol.sol(t)


def Wf_payload_mission(CL, CD, S, eta, SFC, P_max, P_aux, t_takeoff=4*60, h_cruise1=500, h_cruise2=500, W0=const.MTOW, drop_payload=0, R=None):
    if R is None:
        R = np.linspace(0, 2*const.R_cruise, 10)

    P_extra = P_aux + const.P_pay_pay
    if type(W0) == float:
        W = np.full_like(t, W0)
    else:
        W = W0.copy()
    W -= Wf_take_off(P_max, SFC, t_takeoff)
    W -= Wf_climb(h_cruise1, W, SFC)
    W = np.diag(Wf_range(W, R, SFC, eta, CL, CD, S, P_extra, h=h_cruise1))
    W = Wf_loiter(W, [const.T_loiter_pay], SFC, eta, CL, CD, S, P_extra)[:,-1]
    W -= drop_payload * const.g0
    W -= Wf_climb(h_cruise2, W, SFC)
    W = np.diag(Wf_range(W, R, SFC, eta, CL, CD, S, P_extra, h=h_cruise2, v_cruise=0))
    
    return R, W0 - (W + drop_payload * const.g0)


def Wf_endurance_mission(CL, CD, S, eta, SFC, P_max, P_aux, t_takeoff=4*60, h_cruise1=500, h_cruise2=500, W0=const.MTOW, R=const.R_cruise, t=None):
    if t is None:
        t = np.linspace(0, 2*const.T_loiter_end, 10)
    
    P_extra = P_aux + const.P_pay_pay
    if type(W0) == float:
        W = np.full_like(t, W0)
    else:
        W = W0.copy()
    W -= Wf_take_off(P_max, SFC, t_takeoff)
    W -= Wf_climb(h_cruise1, W, SFC)
    W = Wf_range(W, [R], SFC, eta, CL, CD, S, P_extra, h=h_cruise1)[:,-1]
    W = np.diag(Wf_loiter(W, t, SFC, eta, CL, CD, S, P_extra)).copy()
    W -= Wf_climb(max(h_cruise2-const.h_loiter,0), W, SFC)
    W = Wf_range(W, [R], SFC, eta, CL, CD, S, P_extra, h=h_cruise2, v_cruise=0)

    return t, W0 - W


def payload_range_diagram(OEW, Wf_max, Payload_max, CL, CD, S, eta, SFC, P_max, P_aux, t_takeoff=4*60, h_cruise1=500, h_cruise2=500, W0=const.MTOW, drop_payload=False):
    if drop_payload == True:
        raise NotImplementedError
    
    R = np.linspace(0, 2*const.R_cruise, 40)
    R, Wf = Wf_payload_mission(CL, CD, S, eta, SFC, P_max, P_aux, t_takeoff, h_cruise1, h_cruise2, W0=W0, R=R)

    i = np.where(Wf < Wf_max)
    R, Wf = R[i], Wf[i]
    Payload = np.minimum([Payload_max], W0-OEW-Wf)

    plt.plot(R/1e3, Wf/const.g0)
    plt.plot(R/1e3, Payload/const.g0)
    plt.show()

    
