import numpy as np
from typing import Literal
from numpy.typing import NDArray

# Take of phase
def Mfrac_take_off():
    pass


def Mfrac_climb(H: float, SFC: float):
    """
    Parameters
    ----------
    H : float
        Height required [m]
    SFC : float
        Specific fuel consumption [g/kW/h]
    """
    g0 = 9.81  # [m/s^2]
    SFC = SFC / (1000 * 1000 * 3600)  # Convert SFC from g/kW/h to kg/W/s
    Mfrac = np.exp(-H * g0 * SFC)
    return Mfrac


def Mfrac_cruise(SFC: float, R: float, CL: float, CD: float, eta: float) -> float:
    """
    Caculate fuel mass fraction of cruise phase using Breguet range equation
    R = (CL/CD) * (eta/SFC) * ln(M1/M2)

    Parameters
    ----------
    SFC : float
        Specific fuel consumption [g/kW/h]
    R : float
        Range [km]
    CL : float
        Lift coefficient should correspond to the optimum of CL/CD
    CD : float
        Drag coefficient
    eta : float
        Propulsive efficiency

    Returns
    -------
    Mfrac : float
        Fuel mass fraction of cruise phase
    """
    R = R * 1000  # Convert range from km to m
    SFC = SFC / (1000 * 1000 * 3600)  # Convert SFC from g/kW/h to kg/W/s

    Mfrac = np.exp(-SFC * R / (eta * CL/CD))
    return Mfrac


def Mfrac_loiter(E: float, SFC: float, eta, CD, CL, V1):
    """
    Calculate fuel mass fraction of loiter phase
    From "Elements of airplane performance", second edition equation 15.16
    E = 2 * eta * CL / (SFC * CD * V1) * (sqrt(W1 / W2) - 1)
    For optimum endurance, CL and CD should correspond to maximum of CL^3/CD^2

    Parameters
    ----------
    E : float
        Endurance [hr]
    SFC : float
        Specific fuel consumption [g/kW/h]
    eta : float
        Propulsive efficiency
    CD : float
        Drag coefficient
    CL : float
        Lift coefficient
    V1 : float
        Velocity at starting phase of endurence [m/s]
    
    Returns
    -------
    Mfrac : float
        Fuel mass fraction of loiter phase
    """
    E = E * 3600  # Convert endurance from hr to s
    SFC = SFC / (1000 * 1000 * 3600)  # Convert SFC from g/kWh to kg/Ws
    
    Mfrac = (E * SFC * CD * V1 / (2 * eta * CL) + 1)**2
    return Mfrac
    
def total_mass_frac(CL: NDArray, CD: NDArray, SFC: float, S: float, which: Literal['payload', 'endurance']='payload'):
    rho = 1.225  # [kg/m^3]
    g0 = 9.81  # [m/s^2]
    eta = 0.8
    MTOW = 160 # [kg]
    W1 = MTOW * 9.81  # [N]
    v1 = np.sqrt(W1 * 2 / (rho * S * CL))
    vreq = 42 # [m/s]

    CLIMB_HEIGHT = 500  # [m]
    CRUISE_DIST = 185  # [km]
    if which == 'payload':
       LOITER_TIME = 20/60  # [hr] 
    elif which == 'endurance':
        LOITER_TIME = 10  # [hr]

    CL_CD = CL/CD
    CL_cruise = CL[np.argmax(CL_CD)]
    CD_cruise = CD[np.argmax(CL_CD)]
    v_cruise = v1[np.argmax(CL_CD)]
    Pcruise = W1 * CD_cruise / CL_cruise * v_cruise
    print('At maximum CL/CD:')
    print(f'V: {v_cruise:.3f} m/s')
    print(f'CL/CD: {CL_cruise/CD_cruise:.3f}')
    print(f'CL: {CL_cruise:.3f}')
    print(f'P (at MTOW): {Pcruise:.3f} W')

    if v_cruise < vreq:
        print(f'Warning: Cruise velocity is lower than required velocity of {vreq:.3f} m/s')
        print(f'At maximum CL/CD satisfying vreq:')
        i = np.argmax(CL_CD[np.where(v1 > vreq)])
        CL_cruise = CL[i]
        CD_cruise = CD[i]
        v_cruise = v1[i]
        Pcruise = W1 * CD_cruise / CL_cruise * v_cruise
        print(f'V: {v_cruise:.3f} m/s')
        print(f'CL/CD: {CL_cruise/CD_cruise:.3f}')
        print(f'CL: {CL_cruise:.3f}')
        print(f'P (at MTOW): {Pcruise:.3f} W')

    CL_CD = CL**3/CD**2
    CL_loiter = CL[np.argmax(CL_CD)]
    CD_loiter = CD[np.argmax(CL_CD)]
    v_loiter = v1[np.argmax(CL_CD)]
    Ploiter = W1 * CD_loiter / CL_loiter * v_loiter
    print(f'V loiter: {v_loiter:.3f} m/s')
    print(f'CL/CD loiter: {CL_loiter/CD_loiter:.3f}')
    print(f'CL loiter: {CL_loiter:.3f}')
    print(f'P loiter (at MTOW): {Ploiter:.3f} W')


    Mfrac_take_off = 0.99
    Mfrac_climb_1 = Mfrac_climb(CLIMB_HEIGHT, SFC)
    Mfrac_cruise_1 = Mfrac_cruise(SFC, CRUISE_DIST, CL_cruise, CD_cruise, eta)
    print(f'Mfrac_cruise_1: {Mfrac_cruise_1:.3f}')
    
    W_loiter_1 = MTOW * Mfrac_take_off * Mfrac_climb_1 * Mfrac_cruise_1 * g0

    v_loiter = np.sqrt(W_loiter_1 * 2 / (rho * S * CL_loiter))
    print(f'v_loiter: {v_loiter:.3f} m/s') 
    Mfrac_loiter_1 = Mfrac_loiter(LOITER_TIME, SFC, eta, CD_loiter, CL_loiter, v_loiter)
    print(f'Mfrac_loiter_1: {Mfrac_loiter_1:.3f}')

    # cruise back
    Mfrac_cruise_back = Mfrac_cruise(SFC, CRUISE_DIST, CL_cruise, CD_cruise, eta)
    print(f'Mfrac_cruise_back: {Mfrac_cruise_back:.3f}')

    # land
    Mfrac_land = 0.99

    return Mfrac_take_off * Mfrac_climb_1 * Mfrac_cruise_1 * Mfrac_loiter_1 * Mfrac_cruise_back * Mfrac_land


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    # from inputs import *
    from ..aero.cl_cd import *
    

    b = 4  # [m]
    S = 1.5  # [m^2]

    CL, CD, _ = dragpolar(b, S)
    plt.plot(CD, CL)
    plt.show()
    SFC = 300 # [g/kWh]

    Mfrac_fuel = total_mass_frac(CL, CD, SFC, S, which='payload')
    print(f'Mass fraction for fuel: {Mfrac_fuel:.3f}')

