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

step = 1/128
semispan = 3


geom, chords, spans = wgg.wing_geometry('external files/lednicerdatfile.dat', 0.833, 0.4, semispan, 0, step)
loadsx, loadsz, torqueyy, point_range, max_th = ilm.combined_loading(0, semispan, step, False)
moment_distributionx, moment_distributionz, point_range = ilm.moment_distr_from_load_distr(loadsx, loadsz, point_range, 1/128)

Ixx, Iyy, Ixy, area, Aenc, denom = wgg.wingbox_AMOI(geom, 0.001)
normal_stress = bending(np.resize(moment_distributionx, len(Ixx)),np.resize(moment_distributionz, len(Ixx)),Ixx,Iyy,Ixy,denom, geom)

volume = sum(area)*step


torqueyy = np.resize(torqueyy, len(Aenc))
Aenc = 0.5*np.array(Aenc)
shear_flow = []
tspar= 0.0010
Emod = 70 * 10 ** 9

for (tyy, aen) in itertools.zip_longest(torqueyy, Aenc):
    shearflow = abs(0.5*tyy/aen)
    shear_flow.append(shearflow)
nrib = 0
seclength = 0
index = 0
x, y, z = geo_to_cart(geom)
h = 0.12*np.array(chords)
for shear in shear_flow:
    seclength += step
    taucrit = Emod*0.9*5*((tspar/seclength)**2)
    print(shear, taucrit, shear / tspar, h[index], seclength)
    index += 1
    if taucrit <= abs(shear/tspar):
        print('rib', seclength)
        nrib += 1
        seclength = 0
        shear_flow[:] -= shear

plt.plot(np.resize(point_range, len(shear_flow)), shear_flow)
plt.xlabel('semi-spanwise location (m)')
plt.ylabel('spanwise shear flow distirbution (N/m)')
plt.show

#print(volume, mass)


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

