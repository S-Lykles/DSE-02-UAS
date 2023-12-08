import numpy as np
import matplotlib.pyplot as plt
import const

def rotor_sizing_tool(W, DL, N, V_max, psi_rad=20*const.deg2rad, C_T_sig=0.11):
    """
    Rotor sizing tool based on the ppt from Marilena

    Parameters
    ----------
    W : float
        Weight of the aircraft [N]
    DL : float
        Disk loading [N/m^2]
    N : int
        Number of rotors [-]
    V_max : float
        Maximum speed of the aircraft [m/s]
    psi_rad : float, optional
        ...
    C_T_sig : float, optional
        ...
    """
    R               = np.sqrt(W/(N * DL * np.pi))
    V_tip           = 140*(2*R)**0.171
    D_v             = 0.04*W
    k_dl            = 1 + D_v/W
    omega           = V_tip/R
    Vne             = 1.1*V_max
    mu_Vne          = 1.1*Vne/(omega*R)
    Advance_ratio   = V_max / V_tip

    #level flight
    T_level         = W*k_dl
    C_T_level       = T_level/ (const.rho0*np.pi*R**2*omega**2*R**2)
    sig_level       = C_T_level/C_T_sig

    #turning flight
    n_z             = 1 / np.cos(psi_rad)
    T_turn          = W * k_dl * n_z
    C_T_turn        = T_turn / (const.rho0 * np.pi * R**2 * (omega*R)**2)
    sig_turn        = C_T_turn / C_T_sig

    # T_gust = n_z*k_dl*
    # C_T_gust = T_gust/ (Vne*pi*R**2*omega**2*R**2)

    sig_max = max(sig_level, sig_turn)

    return R, D_v, omega, T_level, sig_max


def P_profile_drag(v, W, N, R, omega, sig_max, Cl_alpha_rot=5.73):
    """
    Calculates the profile drag power using the formula from the ppt.

    Parameters
    ----------
    v : np.ndarray
        Velocity of the aircraft [m/s]
    W : float
        Weight of the aircraft [N]
    N : int
        Number of rotors [-]
    R : float
        Radius of the rotor [m]
    omega : float
        Angular velocity of the rotor [rad/s]
    sig_max : float
        Maximum solidity of the rotor [-]
    Cl_alpha_rot : float, optional
        Lift coefficient of the rotor [-]. The default is 5.73
    """
    advance_ratio = v / (omega*R)
    C_t = (W/N) / (const.rho0 * np.pi * R**2 * (omega*R)**2)
    Cl_bar = 6.6*(C_t / sig_max)
    alpha_m = Cl_bar / Cl_alpha_rot

    CD_p1 = 0.0087 - 0.0216*alpha_m + (0.4 * alpha_m **2)
    CD_p2 = 0.011 + 0.4*(alpha_m**2)
    CD_p3 = 0.009 + 0.73*(alpha_m**2)
    C_D_p = (CD_p1 + CD_p3 + CD_p2) / 3

    P_hov = (1/8)*sig_max*C_D_p*const.rho0*(omega*R)**3*np.pi*(R**2)

    return P_hov * (1 + 4.65*(advance_ratio**2)) * N
    

def P_induced(v, DL, W, k=1.15, k_dl=1.04):
    """
    Calculates the induced power using the formula for
    low speed forward flight from the ppt.

    Parameters
    ----------
    DL : float
        Disk loading [N/m^2]
    W : float
        Weight of the aircraft [N]
    k : float, optional
        Inflow correction [-]. Default is 1.15.
    k_dl : float, optional
        Drag correction [-]. Default is 1.04.
    """
    v_ih = np.sqrt(DL / (2 * const.rho0))
    # v_ibar = 1 / v
    # v_i = v_ibar * 
    a = v_ih**-4
    b = (v**2/v_ih**4)
    c = -1
    x = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
    v_i = np.sqrt(x)
    return k * k_dl * W * v_i


def generate_number_of_blades(R, sigma):
    # We assume that we can integrate between 2-20 blades per rotor
    possible_number_blades = np.arange(2, 20, 1)

    # The chord formula follows from the slides from Marilena
    # This assumes a constant chord!
    chord_array = (sigma * np.pi * R)/ possible_number_blades
    AR = (R**2) / (R * chord_array)  # The aspect ratio of the blades

    # For a certain number of blades, the aspect ratio is found.
    # This is constrained between 14<AR<20 as in de slides
    AR_contrained = AR[(AR > 14) & (AR < 20)]  # The aspect ratio as based on the AR constraint
    # The following connects the possible AR to the number of blades based on index
    number_of_blades = possible_number_blades[np.in1d(AR, AR_contrained).nonzero()[0]]
    print('Number of blades:', number_of_blades)


def generate_power_versus_disk_loading(P_hover_array, DL_array):
    # Assuming units in kg/Watts and kg/m^2 respectively
    print('Hover Power:', P_hover_array)
    PL_array = (M_gross * 2.20462) / (P_hover_array * 0.00135962)
    DL_array = DL_array * (0.204816)

    # Now we have arrays in lb/hp and lb/ft^2
    # We plot this to compare with NASA data (See literature)
    plt.scatter(DL_array, PL_array)
    plt.grid()
    plt.xlabel('Disk Loading in lb/ft^2')
    plt.ylabel('Power efficiency in lb/hp')
    plt.show()

def delta_p_climb(vc, W):
    delta_p = (W*0.224808943)*((vc*3.2808399)/2)
    return delta_p


N_blades = generate_number_of_blades(1.4738579917347323, 0.06659266968898687)
print(N_blades)
