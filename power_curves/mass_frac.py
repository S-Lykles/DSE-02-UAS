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
    MTOW = 160 # [kg]
    rho = 1.225  # [kg/m^3]
    g0 = 9.81  # [m/s^2]
    eta = 0.8

    CLIMB_HEIGHT = 500  # [m]
    CRUISE_DIST = 185  # [km]
    if which == 'payload':
       LOITER_TIME = 20/60  # [hr] 
    elif which == 'endurance':
        LOITER_TIME = 10  # [hr]

    CL_CD = CL/CD
    Mcruise = MTOW
    CL_cruise = CL[np.argmax(CL_CD)]
    print(f'CL in cruise: {CL_cruise:.3f}')
    CD_cruise = CD[np.argmax(CL_CD)]
    print(f'CL/CD max: {CL_cruise/CD_cruise:.3f}')
    Vcruise = np.sqrt(Mcruise*g0 * 2 / (rho * S * CL_cruise))
    print(f'Vcruise: {Vcruise:.3f} m/s')
    Pcruise = Mcruise * g0 * CD_cruise / CL_cruise * Vcruise
    print(f'Pcruise: {Pcruise:.3f} W')

    CL_CD = CL**3/CD**2
    CL_loiter = CL[np.argmax(CL_CD)]
    CD_loiter = CD[np.argmax(CL_CD)]

    Mfrac_take_off = 0.99
    Mfrac_climb_1 = Mfrac_climb(CLIMB_HEIGHT, SFC)
    Mfrac_cruise_1 = Mfrac_cruise(SFC, CRUISE_DIST, CL_cruise, CD_cruise, eta)
    print(f'Mfrac_cruise_1: {Mfrac_cruise_1:.3f}')
    
    W1 = MTOW * Mfrac_take_off * Mfrac_climb_1 * Mfrac_cruise_1 * g0

    V1 = np.sqrt(W1 * 2 / (rho * S * CL_loiter))
    print(f'V1: {V1:.3f} m/s') 
    Mfrac_loiter_1 = Mfrac_loiter(LOITER_TIME, SFC, eta, CD_loiter, CL_loiter, V1)
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

    b = 6  # [m]
    S = 3.763  # [m^2]

    CL, CD, _ = dragpolar(b, S)
    plt.plot(CD, CL)
    plt.show()
    SFC = 300 # [g/kWh]

    Mfrac_fuel = total_mass_frac(CL, CD, SFC, S, which='payload')
    print(f'Mass fraction for fuel: {Mfrac_fuel:.3f}')

