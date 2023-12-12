import numpy as np
from plot_setting import *
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import pandas as pd
import const
# import dual_phase.inputs as inputs
from power_curves.mass_frac import fuel_weight
# from aero.cl_cd import dragpolar_dual


def plot_fuel(b, S, polar, b_design=6, S_design=3.763, which='endurance', name=None):
    # plt.style.use('seaborn')
    plt.rcParams.update(tex_fonts)

    fuel_w = np.zeros((len(b), len(S)))
    A = np.zeros((len(b), len(S)))
    unit_conversion = 0.45359237 / 745.7/ 60**2
    SFC = 0.5 * unit_conversion
    for j in range(len(S)):
        for i in range(len(b)):
            CL, CD = polar(b[i], S[j])
            fw = fuel_weight(CL, CD, SFC, S[j], which=which, v_cruise=const.v_cruise*1.1)
            fuel_w[j][i] = fw
            A[j][i] = b[i]**2/S[j]

    original_cmap = plt.get_cmap('viridis')
    new_cmap = mcolors.LinearSegmentedColormap.from_list(
        'truncated_viridis', original_cmap(np.linspace(0.15, 1, 256))
    )
    plt.figure(figsize=set_size(fraction=1, subplots=(2,1)))
    CS = plt.contour(b, S, fuel_w, levels=8, colors='k')
    CSA = plt.contour(b, S, A, levels=8, colors='white', linestyles='dashed', linewidths=0.5)
    Aspect_patch = mpatches.Patch(color='white', linestyle='dashed', label='Aspect ratio', linewidth=0.5)
    A_line = Line2D([0], [0], color='white', linestyle='dashed', linewidth=1)
    plt.clabel(CS, inline=1, colors='k')
    plt.clabel(CSA, inline=1, colors='white')
    cnt = plt.contourf(b, S, fuel_w, levels=100, cmap=new_cmap)
    # Need this rasterization trick because otherwise pdf looks weird: https://stackoverflow.com/questions/15822159/aliasing-when-saving-matplotlib-filled-contour-plot-to-pdf-or-eps
    for c in cnt.collections:
        c.set_rasterized(True)

    cbar = plt.gca().figure.colorbar(cnt, ax = plt.gca())
    cbar.ax.set_ylabel('$M_f$ [kg]', rotation = -90, va = "bottom")
    plt.scatter(b_design, S_design, c='r', marker='x', label='Design point')
    plt.xlabel('$b$ [m]')
    plt.ylabel('$S$ [$\\mathrm{m}^2$]')
    # add Aspect ratio to legend, need some hack because contour doesn't allow a label
    handles, labels = plt.gca().get_legend_handles_labels()
    handles.append(A_line)
    labels.append('Aspect ratio')
    plt.legend(handles=handles, labels=labels)
    plt.tight_layout()
    if name is not None:
        plt.savefig('fuel_plots/'+name)
    plt.show()


def plot_cd0(b, S, CD0func, b_design=6, S_design=3.763):
    cd0 = np.zeros((len(b), len(S)))
    unit_conversion = 0.45359237 / 745.7/ 60**2
    SFC = 0.5 * unit_conversion
    for i in range(len(b)):
        for j in range(len(S)):
            cd0[i][j] = CD0func(b[i], S[j])


    original_cmap = plt.get_cmap('viridis')
    new_cmap = mcolors.LinearSegmentedColormap.from_list(
        'truncated_viridis', original_cmap(np.linspace(0.15, 1, 256))
    )
    plt.figure()
    CS = plt.contour(b, S, cd0, levels=8, colors='k')
    plt.clabel(CS, inline=1, fontsize=10, colors='k')
    plt.contourf(b, S, cd0, levels=100, cmap=new_cmap)
    plt.xlabel('$b$ [m]')
    plt.ylabel('$S$ [$\\mathrm{m}^2$]')
    plt.show()

