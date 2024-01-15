from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
# from DSE.plot_setting import report_tex, set_size, report_fast
from DSE import const

# report_tex = report_fast
# Read Xfoil data from the file
# Dual Phase

# Get the directory where the script is located
file_dir = Path(__file__).parent

# Construct the path to the text file
data_23012_airfoil = np.loadtxt(file_dir/r"NACA23012_RE_1.6E6.txt", skiprows=11)
data_23012_wing = np.loadtxt(file_dir/r"T1-42_0 m_s-LLT.txt", skiprows=8)
data_0012_airfoil = np.loadtxt(file_dir/r"NACA0012_RE_2E6.txt", skiprows=11)
data_0012_wing = np.loadtxt(file_dir/r"naca0012.txt", skiprows=11)

# Extract columns
alpha_23012_airfoil = data_23012_airfoil[:, 0]    # Angle of attack
cl_23012_airfoil = data_23012_airfoil[:, 1]       # Lift coefficient
cd_23012_airfoil = data_23012_airfoil[:, 2]       # Drag coefficient
cm_23012_airfoil = data_23012_airfoil[:, 4]       # Moment coefficient

alpha_23012_wing = data_23012_wing[:, 0]    # Angle of attack
cl_23012_wing = data_23012_wing[:, 2]       # Lift coefficient
cd_23012_wing = data_23012_wing[:, 5]       # Drag coefficient
cm_23012_wing = data_23012_wing[:, 8]       # Moment coefficient

alpha_0012_airfoil = data_0012_airfoil[:, 0]    # Angle of attack
cl_0012_airfoil = data_0012_airfoil[:, 1]       # Lift coefficient
cd_0012_airfoil = data_0012_airfoil[:, 2]       # Drag coefficient
cm_0012_airfoil = data_0012_airfoil[:, 4]       # Moment coefficient


# Plot the data Dual Phase
if __name__ == "__main__":
    # plt.rcParams.update(report_tex)
    # size = set_size(fraction=1, subplots=(1,1))
    # plt.figure(figsize=(size[0]*0.7,size[1]))
    plt.figure(figsize=(10, 6))
    plt.subplot(2,2,1)
    plt.plot(alpha_23012_airfoil, cl_23012_airfoil, label='NACA 23012_airfoil')
    plt.plot(alpha_0012_airfoil, cl_0012_airfoil, label='NACA 0012_airfoil')
    plt.title('Xfoil Data - $C_l$ vs alpha')
    plt.xlabel('Angle of Attack (degrees)')
    plt.ylabel('Lift Coefficient ($C_l$)')
    plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
    plt.gca().grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
    plt.minorticks_on()
    plt.tight_layout()
    plt.legend()

    plt.subplot(2,2,2)
    plt.plot(cd_23012_airfoil, cl_23012_airfoil, label='NACA 23012_airfoil')
    plt.plot(cd_0012_airfoil, cl_0012_airfoil, label='NACA 0012_airfoil')
    plt.title('Xfoil Data - $C_l$ vs $C_d$')
    plt.xlabel('Drag Coefficient ($C_d$)')
    plt.ylabel('Lift Coefficient ($C_l$)')
    plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
    plt.gca().grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
    plt.minorticks_on()
    plt.tight_layout()
    plt.legend()

    plt.subplot(2,2,3)
    plt.plot(alpha_23012_airfoil, cm_23012_airfoil, label='NACA 23012_airfoil')
    plt.plot(alpha_0012_airfoil, cm_0012_airfoil, label='NACA 0012_airfoil')
    plt.title('Xfoil Data - $C_m$ vs alpha')
    plt.xlabel('Angle of Attack (degrees)')
    plt.ylabel('Moment Coefficient ($C_m$)')
    plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
    plt.gca().grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
    plt.minorticks_on()
    plt.tight_layout()
    plt.legend()

    plt.subplot(2,2,4)
    plt.plot(alpha_23012_airfoil, cl_23012_airfoil/cd_23012_airfoil, label='NACA 23012_airfoil')
    plt.plot(alpha_0012_airfoil, cl_0012_airfoil/cd_0012_airfoil, label='NACA 0012_airfoil')
    plt.title('Xfoil Data - $Cl/Cd$ vs alpha')
    plt.xlabel('Angle of Attack (degrees)')
    plt.ylabel('Cl/Cd')
    plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
    plt.gca().grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
    plt.minorticks_on()
    plt.tight_layout()
    plt.legend()

    # plt.savefig('aero/dual_airfoil.pdf')
    plt.show()
