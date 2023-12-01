import numpy as np

# Take of phase
def Mfrac_take_off():
    pass


def Mfrac_climb():
    pass


def Mfrac_cruise(SFC: float, R: float, CL_CD: float, eta: float) -> float:
    """
    Caculate fuel mass fraction of cruise phase using Breguet range equation
    R = (CL/CD) * (eta/SFC) * ln(M1/M2)

    Parameters
    ----------
    SFC : float
        Specific fuel consumption [g/kW/h]
    R : float
        Range [km]
    CL_CD : float
        (optimum) Lift-to-drag ratio
    eta : float
        Propulsive efficiency

    Returns
    -------
    Mfrac : float
        Fuel mass fraction of cruise phase
    """
    R = R * 1000  # Convert range from km to m

    Mfrac = np.exp(-SFC * R / (eta * CL_CD))
    return Mfrac


def Mfrac_loiter(T: float, SFC: float, ):
    pass


