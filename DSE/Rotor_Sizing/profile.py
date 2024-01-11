import numpy as np
import matplotlib.pyplot as plt

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
P0_chord = [0.1, 0.1, 0.17, 0.15, 0.09, 0.07, 0.02]
P_chord_bounds = [(0.05, 0.15), (0.0, 0.6), (0.0, 0.6), (-0.3, 0.5), (-0.3, 0.3), (0.0, 0.2), (0.01, 0.2)]
P0_twist = [14, 7, 5, 4, 3]
P_twist_bounds = [(10, 20), (0, 20), (0, 20), (0, 20), (0, 3)]


def chord_dist(r, P=P0_chord):
    R = r[-1]
    c = bezier_curve(P, r/R) * R
    return c


def twist_dist(r, P=P0_twist):
    R = r[-1]
    t = bezier_curve(P, r/R)
    return np.deg2rad(t)

if __name__ == "__main__":
    r = np.linspace(0, 1, 100)
    c = chord_dist(r)
    t = twist_dist(r)
    fig, ax = plt.subplots(2, 1)
    ax[0].plot(r, c)
    ax[1].plot(r, t)

    plt.show()
