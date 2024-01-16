import numpy as np
from DSE import const
from DSE.power_curves import rotor_tool
from matplotlib import pyplot as plt

R, D_v, omega, T_level, sig_max = rotor_tool.rotor_sizing_tool((180*9.81),500, 4, V_max=10, psi_rad=10*const.deg2rad, C_T_sig=0.16)
print(rotor_tool.generate_number_of_blades(R, sig_max))
print(sig_max)

def find_sigma_for_two_blades(R):
    sigma_values = np.linspace(0.01, 1.0, 1000)

    for sigma in sigma_values:
        number_of_blades = rotor_tool.generate_number_of_blades(R, sigma)
        if 2 in number_of_blades:
            return sigma

    return None

R_value = R
sigma_for_two_blades = find_sigma_for_two_blades(R_value)

print(f"For R={R_value}, sigma={sigma_for_two_blades} results in number_of_blades=2.")

