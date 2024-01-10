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
data_23012 = np.loadtxt(file_dir/r"T1-42_0 m_s-LLT.txt", skiprows=8)

# Extract columns
alpha_23012 = data_23012[:, 0]    # Angle of attack
cl_23012 = data_23012[:, 2]       # Lift coefficient
cd_23012 = data_23012[:, 5]       # Drag coefficient
cm_23012 = data_23012[:, 8]       # Moment coefficient


# Plot the data Dual Phase
if __name__ == "__main__":
    # plt.rcParams.update(report_tex)
    # size = set_size(fraction=1, subplots=(1,1))
    # plt.figure(figsize=(size[0]*0.7,size[1]))
    plt.figure(figsize=(10, 6))
    plt.subplot(1,3,1)
    plt.plot(alpha_23012, cl_23012, label='NACA 23012')
    plt.title('Xfoil Data - $C_l$ vs alpha')
    plt.xlabel('Angle of Attack (degrees)')
    plt.ylabel('Lift Coefficient ($C_l$)')
    plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
    plt.gca().grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
    plt.minorticks_on()
    plt.tight_layout()
    plt.legend()

    plt.subplot(1,3,2)
    plt.plot(cd_23012, cl_23012, label='NACA 23012')
    plt.title('Xfoil Data - $C_l$ vs $C_d$')
    plt.xlabel('Drag Coefficient ($C_d$)')
    plt.ylabel('Lift Coefficient ($C_l$)')
    plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
    plt.gca().grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
    plt.minorticks_on()
    plt.tight_layout()
    plt.legend()

    plt.subplot(1,3,3)
    plt.plot(alpha_23012, cm_23012, label='NACA 23012')
    plt.title('Xfoil Data - $C_m$ vs alpha')
    plt.xlabel('Angle of Attack (degrees)')
    plt.ylabel('Moment Coefficient ($C_m$)')
    plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
    plt.gca().grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
    plt.minorticks_on()
    plt.tight_layout()
    plt.legend()

    # plt.savefig('aero/dual_airfoil.pdf')
    plt.show()
