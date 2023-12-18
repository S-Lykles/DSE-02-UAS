import numpy as np
from typing import Literal
from numpy.typing import NDArray
from DSE import const

    
def fuel_weight(CL, CD, SFC, S, which='payload', v_cruise=const.v_cruise*1.1, eta=0.75, energy_fraction=1.):
    """
    Calculate fuel weight used for endurance or payload mission
    This is only the weight for cruising and loitering, not including take-off and climb

    Parameters
    ----------
    CL : numpy.ndarray
        Lift coefficient
    CD : numpy.ndarray
        Drag coefficient
    SFC : float
        Specific fuel consumption [kg/W/s]
    S : float
        Wing area [m^2]
    which : str, optional
        'payload' or 'endurance'. The default is 'payload'.
    v_cruise : float, optional
        Cruise velocity [m/s]. The default is minimum cruise speed tines 1.1.
    Energy_fraction : float, optional
        Fraction fuel which is subtracted after each phase. The default is 1. In case
        you have batteries then this would be 0, in the case of hydrogen something like 0.1
    """
    if which == 'payload':
        P_pay = const.P_pay_pay
    elif which == 'endurance':
        P_pay = const.P_pay_end
    
    W1 = const.MTOW
    v1 = np.sqrt(W1 * 2 / (const.rho0 * S * CL))
    D = W1 * CD / CL
    P = D * v1 + const.P_aux + P_pay
    E = const.R_cruise / v1 * P
    idx = np.argmin(E[np.where(v1 > v_cruise)])
    E = E[idx]
    M_fuel1 = E * SFC / eta

    W2 = W1 - M_fuel1 * const.g0 * energy_fraction

    v2 = np.sqrt(W2 * 2 / (const.rho0 * S * CL))
    D = W2 * CD / CL
    P = D * v2 + const.P_aux + P_pay
    if which == 'payload':
        E = const.T_loiter_pay * P
    elif which == 'endurance':
        E = const.T_loiter_end * P
    
    idx = np.argmin(E)
    M_fuel2 = E[idx] * SFC / eta

    W3 = W2 - M_fuel2 * const.g0 * energy_fraction
    
    if which == 'payload':
        W3 -= const.Payload

    v3 = np.sqrt(W3 * 2 / (const.rho0 * S * CL))
    D = W3 * CD / CL
    P = D * v3 + const.P_aux + P_pay
    E = const.R_cruise / v3 * P
    idx = np.argmin(E)
    M_fuel3 = E[idx] * SFC / eta

    return  M_fuel1 + M_fuel2 + M_fuel3

if __name__ == '__main__':
    from DSE.aero.cl_cd import dragpolar_dual
    CL, CD = dragpolar_dual(6, 3.7)
    print(fuel_weight(CL, CD, 320 / (1000 * 1000 * 3600), 6, 'endurance'))

