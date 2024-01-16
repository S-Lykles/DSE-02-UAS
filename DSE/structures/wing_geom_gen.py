from DSE import const
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd

def wing_geometry(file_name, c_root, taper, semispan, twist, stepsize):
    file_dir = Path(__file__).parent
    data = pd.read_table(file_dir/f"{file_name}",delim_whitespace=True,skiprows=[0],names=['x','y'],index_col=False)
    data.drop(index=data.index[0], axis=0, inplace=True)
    span_range = np.arange(0, semispan, stepsize)
    scaling = (1 - taper)/semispan
    chord_list = []
    geom = []
    ones = np.ones(len(data.x))
    ax = plt.axes(projection='3d')
    for semispan in span_range:
        chord = c_root*(1-(semispan*scaling))
        chord_list.append(chord)
        profloc = [chord * data.x - chord * 0.25, chord * data.y, ones*semispan]
        ax.plot3D(ones*semispan, profloc[0], profloc[1], 'blue')
        geom.append(profloc)
    ax.axes.set_xlim3d(left=0, right=3)
    ax.axes.set_ylim3d(bottom=0, top=3)
    ax.axes.set_zlim3d(bottom=0, top=3)
    plt.show()

    return geom, chord_list, span_range




def wingbox_geom(geom, chords, cstart, cend):
    1


def wingbox_AMOI(geom, t):
    Ixx = []
    Iyy = []
    Ixy = []
    area = []
    denom = []
    Aenc = []
    for i in geom:
        xdif = i[0].diff()
        ydif = i[1].diff()

        segment_length = np.sqrt(xdif ** 2 + ydif ** 2)
        segment_angle = np.arctan(xdif / ydif)
        segment_area = t*segment_length
        segment_dx = i[0] + xdif/2
        segment_dy = i[1] + ydif/2
        segment_enclosed = xdif*i[0]+xdif*ydif/2
        Ixx_section = (t*(segment_length ** 3)/12)*np.sin(segment_angle)*np.sin(segment_angle) + (segment_area*segment_dx*segment_dx)
        Iyy_section = (t*(segment_length ** 3)/12)*np.cos(segment_angle)*np.cos(segment_angle) + (segment_area*segment_dy*segment_dy)
        Ixy_section = (t*(segment_length ** 3)/12)*np.sin(segment_angle)*np.cos(segment_angle)
        Ixxcross = Ixx_section.sum()
        Iyycross = Iyy_section.sum()
        Ixycross = Ixy_section.sum()
        Aenclsec = segment_enclosed.sum()
        denomcross = Ixxcross * Iyycross - (Ixycross ** 2)
        Ixx.append(Ixxcross)
        Iyy.append(Iyycross)
        Ixy.append(Ixycross)
        Aenc.append(Aenclsec)
        denom.append(denomcross)
        area.append(segment_area.sum())
    return Ixx, Iyy, Ixy, area, Aenc, denom