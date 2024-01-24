import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from DSE.plot_setting import report_tex, set_size, report_fast

airfoil_selection = False
twist = True
airfoil_variation = False
empennage = False
if airfoil_selection:
    # Read Xfoil data from files
    file_dir = Path(__file__).parent
    airfoils = ['Clark_YH', 'naca_23012', 'naca_43012A', 'RG15']
    data = {name: np.loadtxt(file_dir / f"{name}_data.txt", skiprows=11) for name in airfoils}

    # Plot data
    plt.figure(figsize=(15, 5))

    # First subplot
    plt.subplot(1, 3, 1)
    for name, values in data.items():
        plt.plot(values[:, 0], values[:, 1], label=name)
    plt.title('$C_l$ vs $\\alpha$')
    plt.xlabel('Angle of Attack (deg)')
    plt.ylabel('Lift Coefficient (-)')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.minorticks_on()
    plt.tight_layout()
    plt.legend()

    # Second subplot
    plt.subplot(1, 3, 2)
    for name, values in data.items():
        plt.plot(values[:, 0], values[:, 1] / values[:, 2], label=name)
    plt.title('$C_l/C_D$ vs $\\alpha$')
    plt.xlabel('Angle of attack (deg)')
    plt.ylabel('$C_l/C_D$ (-)')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.minorticks_on()
    plt.tight_layout()
    plt.legend()

    # Third subplot
    plt.subplot(1, 3, 3)
    for name, values in data.items():
        plt.plot(values[:, 0], values[:, 4], label=name)
    plt.title('$C_m$ vs $\\alpha$')
    plt.xlabel('Angle of Attack (deg)')
    plt.ylabel('Moment Coefficient (-)')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.minorticks_on()
    plt.tight_layout()
    plt.legend()

    plt.show()
if twist:
    # Read Xfoil data from files
    file_dir = Path(__file__).parent
    twists = ['0twist', '+1twist', '+2twist', '+3twist', '-1twist', '-2twist', '-3twist']
    data = {name: np.loadtxt(file_dir / f"{name}.txt", skiprows=8) for name in twists}

    # Plot data
    plt.figure(figsize=(15, 5))

    # First subplot
    plt.subplot(1, 3, 1)
    for name, values in data.items():
        plt.plot(values[0:33, 0], values[0:33, 2], label=name)
    plt.title('$C_l$ vs $\\alpha$')
    plt.xlabel('Angle of Attack (deg)')
    plt.ylabel('Lift Coefficient (-)')
    plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
    plt.gca().grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
    plt.minorticks_on()
    plt.tight_layout()
    plt.legend()

    # Second subplot
    plt.subplot(1, 3, 2)
    for name, values in data.items():
        plt.plot(values[0:33:, 0], values[0:33, 2] / values[0:33, 5], label=name)
    plt.title('$C_l/C_D$ vs $\\alpha$')
    plt.xlabel('Angle of attack (deg)')
    plt.ylabel('$C_l/C_D$ (-)')
    plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
    plt.gca().grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
    plt.minorticks_on()
    plt.tight_layout()
    plt.legend()

    # Third subplot
    plt.subplot(1, 3, 3)
    for name, values in data.items():
        plt.plot(values[0:33, 0], values[0:33, 8], label=name)
    plt.title('$C_m$ vs $\\alpha$')
    plt.xlabel('Angle of Attack (deg)')
    plt.ylabel('Moment Coefficient (-)')
    plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
    plt.gca().grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
    plt.minorticks_on()
    plt.tight_layout()
    plt.legend()

    plt.show()
if airfoil_variation:

    # Read Xfoil data from files
    file_dir = Path(__file__).parent
    twists = ['Normal', 'Gradually', 'Suddenly']
    data = {name: np.loadtxt(file_dir / f"{name}.txt", skiprows=8) for name in twists}

    # Plot data
    plt.figure(figsize=(15, 5))

    # First subplot
    plt.subplot(1, 3, 1)
    for name, values in data.items():
        plt.plot(values[:, 0], values[:, 2], label=name)
    plt.title('$C_l$ vs $\\alpha$')
    plt.xlabel('Angle of Attack (deg)')
    plt.ylabel('Lift Coefficient (-)')
    plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
    plt.gca().grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
    plt.minorticks_on()
    plt.tight_layout()
    plt.legend()

    # Second subplot
    plt.subplot(1, 3, 2)
    for name, values in data.items():
        plt.plot(values[:, 0], values[:, 2] / values[:, 5], label=name)
    plt.title('$C_l/C_D$ vs $\\alpha$')
    plt.xlabel('Angle of attack (deg)')
    plt.ylabel('$C_l/C_D$ (-)')
    plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
    plt.gca().grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
    plt.minorticks_on()
    plt.tight_layout()
    plt.legend()

    # Third subplot
    plt.subplot(1, 3, 3)
    for name, values in data.items():
        plt.plot(values[:, 0], values[:, 8], label=name)
    plt.title('$C_m$ vs $\\alpha$')
    plt.xlabel('Angle of Attack (deg)')
    plt.ylabel('Moment Coefficient (-)')
    plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
    plt.gca().grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
    plt.minorticks_on()
    plt.tight_layout()
    plt.legend()

    plt.show()

if empennage:
    # Read Xfoil data from files
    file_dir = Path(__file__).parent
    twists = ['NACA0010', 'NACA0012']
    data = {name: np.loadtxt(file_dir / f"{name}.txt", skiprows=11) for name in twists}

    # Plot data
    plt.figure(figsize=(15, 5))

    # First subplot
    plt.subplot(1, 3, 1)
    for name, values in data.items():
        plt.plot(values[:, 0], values[:, 1], label=name)
    plt.title('$C_l$ vs $\\alpha$')
    plt.xlabel('Angle of Attack (deg)')
    plt.ylabel('Lift Coefficient (-)')
    plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
    plt.gca().grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
    plt.minorticks_on()
    plt.tight_layout()
    plt.legend()

    # Second subplot
    plt.subplot(1, 3, 2)
    for name, values in data.items():
        plt.plot(values[:, 0], values[:, 1] / values[:, 2], label=name)
    plt.title('$C_l/C_D$ vs $\\alpha$')
    plt.xlabel('Angle of attack (deg)')
    plt.ylabel('$C_l/C_D$ (-)')
    plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
    plt.gca().grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
    plt.minorticks_on()
    plt.tight_layout()
    plt.legend()

    # Third subplot
    plt.subplot(1, 3, 3)
    for name, values in data.items():
        plt.plot(values[:, 0], values[:, 4], label=name)
    plt.title('$C_m$ vs $\\alpha$')
    plt.xlabel('Angle of Attack (deg)')
    plt.ylabel('Moment Coefficient (-)')
    plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
    plt.gca().grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
    plt.minorticks_on()
    plt.tight_layout()
    plt.legend()

    plt.show()