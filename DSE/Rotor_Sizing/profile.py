import numpy as np
import matplotlib.pyplot as plt
from DSE.Rotor_Sizing.airfoil import M_alpha
from DSE import const

# Function to calculate binomial coefficient
def binomial_coefficient(n, k):
    return np.math.factorial(n) / (np.math.factorial(k) * np.math.factorial(n - k))

# Function to calculate a point on the Bezier curve
def bezier_curve(P, t):
    n = len(P) - 1
    x = np.zeros_like(t)
    for i, point in enumerate(P):
        bernstein = binomial_coefficient(n, i) * (t ** i) * ((1 - t) ** (n - i))
        x += point * bernstein
    return x


# Define your control points
xb = np.linspace(0, 1, 7)
P0_chord = [0.1, 0.1, 0.17, 0.15, 0.09, 0.07, 0.07, 0.07, 0.1]
P_chord_bounds = [(0.05, 0.15), (0.0, 0.6), (0.0, 0.6), (-0.3, 0.5), (-0.3, 0.3), (-0.1, 0.2), (0.0, 0.2), (0.0, 0.2), (0.0, 0.2)]
P0_twist = [13, 13, 7, 0.5]
P_twist_bounds = [(3, 13), (0,15), (2,7), (0, 0.5)]


def chord_dist(r, P=P0_chord):
    R = r[-1]
    c = bezier_curve(P, r/R) * R
    return c

def twist_dist(r, P=P0_twist):
    R = r[-1]
    t = bezier_curve(P, r/R)
    return np.deg2rad(t)

def alpha_dist(r, P=P0_twist):
    R = r[-1]
    t = bezier_curve(P, r/R)
    return np.deg2rad(t)

if __name__ == "__main__":
    omega = 4000 * 2 * np.pi / 60
    R = 0.5
    r = np.linspace(0, R, 100)
    M, alpha = M_alpha.T
    c = chord_dist(r)
    print(np.mean(c))
    a = alpha_dist(r)
    fig, ax = plt.subplots(2, 1)
    ax[0].plot(r, c)
    ax[1].scatter(np.linspace(0, R, len(P_twist_bounds)),list(zip(*P_twist_bounds))[0])
    ax[1].scatter(np.linspace(0, R, len(P_twist_bounds)),list(zip(*P_twist_bounds))[1])
    ax[1].plot(r, np.rad2deg(a))
    ax[1].plot(M*const.a/omega, alpha)
    ax[1].set_xlim(0, 0.5)

    plt.show()
