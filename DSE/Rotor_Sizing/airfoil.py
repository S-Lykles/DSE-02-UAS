import numpy as np
from pathlib import Path
import re
from DSE import const
from scipy.interpolate import SmoothBivariateSpline
from matplotlib import pyplot as plt
import matplotlib as mpl

x = np.array([])
y = np.array([])
cl_arr = np.array([])
cd_arr = np.array([])
cp_min_arr = np.array([])
airfoil_folder = Path(__file__).parent/'Clark_Y'
for file in airfoil_folder.iterdir():
    match = re.search(r"M(\d+\.\d+)", file.name)
    if not match:
        raise ValueError(f"Could not find Mach number in file name: {file.name}")
    M = float(match.group(1))
    data = np.genfromtxt(file, skip_header=2)
    alpha = data[:, 0]
    cl = data[:, 1]
    cd = data[:, 2]
    cp_min = data[:, 7]

    x = np.append(x, M*np.ones_like(alpha))
    y = np.append(y, alpha)
    cl_arr = np.append(cl_arr, cl)
    cd_arr = np.append(cd_arr, cd)
    cp_min_arr = np.append(cp_min_arr, cp_min)

cl_interp = SmoothBivariateSpline(x, y, cl_arr)# bbox=[0, 0.8, -3, 20])
cd_interp = SmoothBivariateSpline(x, y, cd_arr)# bbox=[0, 0.8, -3, 20])
cp_min_interp = SmoothBivariateSpline(x, y, cp_min_arr)#, bbox=[0, 0.8, -3, 20])

Cl_func_clarky = lambda M, a: cl_interp(M, a, grid=False)
Cd_func_clarky = lambda M, a: cd_interp(M, a, grid=False)
Cp_func_clarky = lambda M, a: cp_min_interp(M, a, grid=False)

if __name__ == "__main__":
    Cpfrac = (1 + (const.gamma-1)/2*M**2) / (1 + (const.gamma-1)/2)
    Cp = 2 / const.gamma/M**2 * (Cpfrac**(const.gamma/(const.gamma-1))-1)

    
    alpha = np.linspace(-2, 20, 100)
    M_grid = np.linspace(0.1, 0.7, 100)
    M_grid, alpha = np.meshgrid(M_grid, alpha)
    cl = Cp_func_clarky(M_grid, alpha)
    fig, ax = plt.subplots()
    # ax.contourf(M_grid, alpha, cl, levels=100)
    CS = ax.contour(M_grid, alpha, cl, colors='k')
    ax.clabel(CS, inline=1, fontsize=10)
    ax.set_xlabel("M")
    ax.set_ylabel(r"$\alpha$")
    ax.set_title("Clark Y")
    plt.show()