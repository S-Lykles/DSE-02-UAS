import numpy as np
import matplotlib.pyplot as plt

def create_boom(N, w, h):
    bottom = np.column_stack((np.linspace(-w/2, w/2, int(N/4)), -h/2*np.ones(int(N/4))))
    top = np.column_stack((np.linspace(-w/2, w/2, int(N/4)), h/2*np.ones(int(N/4))))
    left = np.column_stack((-w/2*np.ones(int(N/4)), np.linspace(-h/2, h/2, int(N/4))))
    right = np.column_stack((w/2*np.ones(int(N/4)), np.linspace(-h/2, h/2, int(N/4))))

    boom = np.vstack((bottom, right, top, left))
    return boom

def sec_prop():

    return Ixx, Iyy

boom = create_boom(400, 6, 2)
plt.scatter(boom[:,0], boom[:,1])
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.grid()
plt.show()
