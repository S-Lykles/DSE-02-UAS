import numpy as np
from typing import Literal
from numpy.typing import NDArray
from . import inputs

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
    CL : flon of loiter phase
    From "Elements of airplane performance", second edition equation 15.16
    E = 2 * eta * CL / (SFC * CD * V1) * (sqrt(W1 / W2) - 1)
    For optimum endurance, CL and CD should correspond to maximum of CL^3/CD^2

    Parametersoat
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
    
def fuel_weight(CL, CD, SFC, S, which='payload'):
    """
    Calculate fuel weight used for endurance or payload mission
    This is only the weight for cruising and loitering, not including take-off and climb
    """



    return 0
    