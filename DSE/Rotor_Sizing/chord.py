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
P0 = [0.14, 0.07, 0.17, 0.15, 0.09, 0.07, 0.02]
P_bounds = [(0.05, 0.15), (0.0, 0.6), (0.0, 0.6), (0.0, 0.5), (0.0, 0.3), (0.0, 0.2), (0.0, 0.2)]

def chord_dist(r, P=P0):
    R = r[-1]
    c = bezier_curve(P, r/R) * R
    return c

