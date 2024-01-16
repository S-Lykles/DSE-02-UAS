print('Structures <3 !!!')
from DSE import const
import numpy as np
import matplotlib.pyplot as plt
import wing_geom_gen as wgg
import Internal_load_model as ilm
from pathlib import Path
import pandas as pd
import itertools

def bending(moment_distributionx, moment_distributionz, Ixx, Iyy, Ixy, denom, geom):
    ybend = moment_distributionx*Iyy - moment_distributionz*Ixy
    xbend = moment_distributionz*Ixx - moment_distributionx*Ixy
    normal_stress = []

    for (yb, xb, d, geo) in itertools.zip_longest(ybend, xbend, denom, geom):
        local_norm = (geo[0]*yb + geo[1]*xb)/d
        #print(local_norm)
        normal_stress.append(local_norm)

    return normal_stress

def geo_to_cart(geom):
    x = []
    y = []
    z = []
    for geo in geom:
        #print(geo)
        x.extend(geo[0].to_numpy())
        y.extend(geo[1].to_numpy())
        z.extend(geo[2])

    #print(len(x), len(y), len(z))
    return x, y, z

step = 0.1
semispan = 3


geom, chords, spans = wgg.wing_geometry('external files/lednicerdatfile.dat', 0.833, 0.4, semispan, 3, step)
loadsx, loadsz, torqueyy, point_range, max_th = ilm.combined_loading(0, semispan, step, False)
moment_distributionx, moment_distributionz, point_range = ilm.moment_distr_from_load_distr(loadsx, loadsz, point_range, 0.01)

Ixx, Iyy, Ixy, area, Aenc, denom = wgg.wingbox_AMOI(geom, 0.001)
normal_stress = bending(np.resize(moment_distributionx, len(Ixx)),np.resize(moment_distributionz, len(Ixx)),Ixx,Iyy,Ixy,denom, geom)

volume = sum(area)*step


torqueyy = np.resize(torqueyy, len(Aenc))
Aenc = np.array(Aenc)
mass = 2700*volume
shear_flow = 0.5*(torqueyy/Aenc)


plt.plot(point_range, shear_flow)
plt.show

print(volume, mass)

x, y, z = geo_to_cart(geom)
plt.rcParams["figure.figsize"] = [7.00, 3.50]
plt.rcParams["figure.autolayout"] = True
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
c = normal_stress
img = ax.scatter(x, y, z, c=c, cmap='YlOrRd', alpha=1)
plt.colorbar(img)
ax.axes.set_xlim3d(left=0, right=3)
ax.axes.set_ylim3d(bottom=0, top=3)
ax.axes.set_zlim3d(bottom=0, top=3)
plt.show()

