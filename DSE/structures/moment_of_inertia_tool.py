from DSE import const
import numpy as np
import math
import matplotlib.pyplot as plt
from pathlib import Path
from Internal_load_model import momentdistributionz

# This tool can be used to compute the moment of inertia of various cross-sectional shapes
# The formula's were obtained from course AE2135-I or


# Function computing moment of inertia around x-axis with cg on neutral line
def moment_of_inertia_rectangle_x_axis(width, height):
    I_xx = width * height**3 / 12

    return I_xx


# Function computing moment of inertia around y-axis with cg on neutral line
def moment_of_inertia_rectangle_y_axis(width, height):
    I_yy = height * width**3 / 12

    return I_yy

# Function computing the moment of inertia of a solid circular section
def moment_of_inertia_solid_circular_section(diameter):
    I = np.pi * diameter**4 / 64

    return I

# Function computing the moment of inertia of a thin-walled circular section
def moment_of_inertia_thin_walled_circular_section(diameter, thickness):
    I = np.pi * thickness * diameter**3 / 8

    return I

# Parallel axis theorem
def parallel_axis_theorem(inertia, area, distance):
    I = inertia + area*distance**2

    return I

# Function for final moment of inertia for the wing box around the x-axis
# h_front_center = Height of the center section of the front spar
# w_front_center = Width of the center section of the front spar
# h_front_tips = Height of the tip section of the front spar
# w_front_tip = Width of the tip section of the front spar
# h_rear_center = Height of the center section of the rear spar
# w_rear_center = Width of the center section of the rear spar
# h_rear_tips = Height of the tip section of the rear spar
# w_rear_tip = Width of the tip section of the rear spar
# t_top = Thickness of the top sheet
# l_top = Length of the top sheet
# t_bottom = Thickness of the bottom sheet
# l_bottom = Length of the bottom sheet
# d = Distance between the neutral line (y-axis) and the center of the spars

def I_xx_wing_box(h_front_center, w_front_center, h_front_tip, w_front_tip, h_rear_center,
                  w_rear_center, h_rear_tip, w_rear_tip, h_center_center, w_center_center, h_center_tip,
                  w_center_tip, t_top, l_top, t_bottom, l_bottom, number_of_spars):
    if number_of_spars == 1:

        # Moment of inertia of various rectangles in the wing box
        I_spar_center = moment_of_inertia_rectangle_x_axis(w_front_center, h_front_center)
        I_spar_tips = moment_of_inertia_rectangle_x_axis(w_front_tip, h_front_tip)
        I_top_sheet = moment_of_inertia_rectangle_x_axis(l_top, t_top)
        I_bottom_sheet = moment_of_inertia_rectangle_x_axis(l_bottom, t_bottom)

        # Including the parallel axis theorem

        # Moment of inertia for tips of spar
        I_spar_tips = parallel_axis_theorem(I_spar_tips, h_front_tip * w_front_tip,
                                             0.5 * h_front_center + 0.5 * h_front_tip)
        # Moment of inertia top and bottom sheet
        I_top_sheet = parallel_axis_theorem(I_top_sheet, t_top * l_top,
                0.5 * h_front_center + h_front_tip + 0.5 * t_top)
        I_bottom_sheet = parallel_axis_theorem(I_bottom_sheet, t_bottom * l_bottom,
                                               0.5 * h_front_center + h_front_tip + 0.5* t_bottom)

        # Adding various components to get the total moment of inertia of the wing box around the x-axis
        I_xx = I_spar_center + 2 * I_spar_tips + I_top_sheet + I_bottom_sheet
        A = l_top*t_top + l_bottom*t_bottom + 2* w_front_tip * h_front_tip + w_front_center * h_front_center


    elif number_of_spars == 2:
        # Moment of inertia of various rectangles in the wing box
        I_front_center = moment_of_inertia_rectangle_x_axis(w_front_center, h_front_center)
        I_front_tips = moment_of_inertia_rectangle_x_axis(w_front_tip, h_front_tip)
        I_rear_center = moment_of_inertia_rectangle_x_axis(w_rear_center, h_rear_center)
        I_rear_tips = moment_of_inertia_rectangle_x_axis(w_rear_tip, h_rear_tip)
        I_top_sheet = moment_of_inertia_rectangle_x_axis(l_top, t_top)
        I_bottom_sheet = moment_of_inertia_rectangle_x_axis(l_bottom, t_bottom)

        # Including the parallel axis theorem

        # Moment of inertia for tips of spar
        I_front_tips = parallel_axis_theorem(I_front_tips, h_front_tip*w_front_tip, 0.5 * h_front_center + 0.5 * h_front_tip)
        I_rear_tips = parallel_axis_theorem(I_rear_tips, h_rear_tip*w_rear_tip, 0.5 * h_rear_center + 0.5 * h_rear_tip)

        #Moment of inertia top and bottom sheet
        I_top_sheet = parallel_axis_theorem(I_top_sheet, t_top*l_top, (
                    0.5 * h_front_center + h_front_tip + t_top + 0.5 * h_rear_center + h_rear_tip) / 2)
        I_bottom_sheet = parallel_axis_theorem(I_bottom_sheet, t_bottom*l_bottom,
                                               (0.5 * h_front_center + h_front_tip + t_bottom + 0.5 * h_rear_center + h_rear_tip) / 2)

        # Adding various components to get the total moment of inertia of the wing box around the x-axis
        I_xx = I_front_center + 2 * I_front_tips + I_rear_center + 2 * I_rear_tips + I_top_sheet + I_bottom_sheet
        A = (l_top * t_top + l_bottom * t_bottom + 2 * w_front_tip * h_front_tip + w_front_center * h_front_center
             + 2 * w_rear_tip * h_rear_tip + w_rear_center * h_rear_center)

    elif number_of_spars == 3:
        # Moment of inertia of various rectangles in the wing box
        I_front_center = moment_of_inertia_rectangle_x_axis(w_front_center, h_front_center)
        I_front_tips = moment_of_inertia_rectangle_x_axis(w_front_tip, h_front_tip)
        I_center_center = moment_of_inertia_rectangle_x_axis(w_center_center, h_center_center)
        I_center_tips = moment_of_inertia_rectangle_x_axis(w_center_tip, h_center_tip)
        I_rear_center = moment_of_inertia_rectangle_x_axis(w_rear_center, h_rear_center)
        I_rear_tips = moment_of_inertia_rectangle_x_axis(w_rear_tip, h_rear_tip)
        I_top_sheet = moment_of_inertia_rectangle_x_axis(l_top, t_top)
        I_bottom_sheet = moment_of_inertia_rectangle_x_axis(l_bottom, t_bottom)

        # Including the parallel axis theorem

        # Moment of inertia for tips of spar
        I_front_tips = parallel_axis_theorem(I_front_tips, h_front_tip * w_front_tip,
                                             0.5 * h_front_center + 0.5 * h_front_tip)
        I_center_tips = parallel_axis_theorem(I_center_tips, h_center_tip * w_center_tip,
                                             0.5 * h_center_center + 0.5 * h_center_tip)
        I_rear_tips = parallel_axis_theorem(I_rear_tips, h_rear_tip * w_rear_tip,
                                            0.5 * h_rear_center + 0.5 * h_rear_tip)

        # Moment of inertia top and bottom sheet
        I_top_sheet = parallel_axis_theorem(I_top_sheet, t_top * l_top, (
                0.5 * h_front_center + h_front_tip + t_top + 0.5 * h_rear_center + h_rear_tip) / 2)
        I_bottom_sheet = parallel_axis_theorem(I_bottom_sheet, t_bottom * l_bottom,
                                               (0.5 * h_front_center + h_front_tip + 1.5 * t_bottom + 0.5 * h_rear_center
                                                + h_rear_tip + h_center_tip + h_center_center) / 3)

        # Adding various components to get the total moment of inertia of the wing box around the x-axis
        I_xx = (I_front_center + 2 * I_front_tips + I_center_center + 2 * I_center_tips  + I_rear_center + 2 * I_rear_tips
                + I_top_sheet + I_bottom_sheet)
        A = (l_top * t_top + l_bottom * t_bottom + 2 * w_front_tip * h_front_tip + w_front_center * h_front_center
             + 2 * w_rear_tip * h_rear_tip + w_rear_center * h_rear_center + 2 * w_center_tip * h_center_tip + w_center_center * h_center_center)

    return I_xx, A


def I_yy_wing_box(h_front_center, w_front_center, h_front_tip, w_front_tip, h_rear_center,
                  w_rear_center, h_rear_tip, w_rear_tip, h_center_center, w_center_center, h_center_tip,
                  w_center_tip, t_top, l_top, t_bottom, l_bottom, number_of_spars, d):
    if number_of_spars == 1:

        # Moment of inertia of various rectangles in the wing box
        I_spar_center = moment_of_inertia_rectangle_y_axis(w_front_center, h_front_center)
        I_spar_tips = moment_of_inertia_rectangle_y_axis(w_front_tip, h_front_tip)
        I_top_sheet = moment_of_inertia_rectangle_y_axis(l_top, t_top)
        I_bottom_sheet = moment_of_inertia_rectangle_y_axis(l_bottom, t_bottom)

        # Adding various components to get the total moment of inertia of the wing box around the x-axis
        I_yy = I_spar_center + 2 * I_spar_tips + I_top_sheet + I_bottom_sheet
        A = l_top * t_top + l_bottom * t_bottom + 2 * w_front_tip * h_front_tip + w_front_center * h_front_center

    elif number_of_spars == 2:
        # Moment of inertia of various rectangles in the wing box
        I_front_center = moment_of_inertia_rectangle_y_axis(w_front_center, h_front_center)
        I_front_tips = moment_of_inertia_rectangle_y_axis(w_front_tip, h_front_tip)
        I_rear_center = moment_of_inertia_rectangle_y_axis(w_rear_center, h_rear_center)
        I_rear_tips = moment_of_inertia_rectangle_y_axis(w_rear_tip, h_rear_tip)
        I_top_sheet = moment_of_inertia_rectangle_y_axis(l_top, t_top)
        I_bottom_sheet = moment_of_inertia_rectangle_y_axis(l_bottom, t_bottom)

        # Including the parallel axis theorem

        # Moment of inertia for the spars
        I_front_center = parallel_axis_theorem(I_front_center, h_front_center * w_front_center, d)
        I_front_tips = parallel_axis_theorem(I_front_tips, h_front_tip * w_front_tip, d)
        I_rear_center = parallel_axis_theorem(I_rear_center, h_rear_center * w_rear_center, d)
        I_rear_tips = parallel_axis_theorem(I_rear_tips, h_rear_tip * w_rear_tip, d)


        # Adding various components to get the total moment of inertia of the wing box around the x-axis
        I_yy = I_front_center + 2 * I_front_tips + I_rear_center + 2 * I_rear_tips + I_top_sheet + I_bottom_sheet
        A = (l_top * t_top + l_bottom * t_bottom + 2 * w_front_tip * h_front_tip + w_front_center * h_front_center
             + 2 * w_rear_tip * h_rear_tip + w_rear_center * h_rear_center)

    elif number_of_spars == 3:
        # Moment of inertia of various rectangles in the wing box
        I_front_center = moment_of_inertia_rectangle_y_axis(w_front_center, h_front_center)
        I_front_tips = moment_of_inertia_rectangle_y_axis(w_front_tip, h_front_tip)
        I_center_center = moment_of_inertia_rectangle_y_axis(w_center_center, h_center_center)
        I_center_tips = moment_of_inertia_rectangle_y_axis(w_center_tip, h_center_tip)
        I_rear_center = moment_of_inertia_rectangle_y_axis(w_rear_center, h_rear_center)
        I_rear_tips = moment_of_inertia_rectangle_y_axis(w_rear_tip, h_rear_tip)
        I_top_sheet = moment_of_inertia_rectangle_y_axis(l_top, t_top)
        I_bottom_sheet = moment_of_inertia_rectangle_y_axis(l_bottom, t_bottom)

        # Including the parallel axis theorem

        # Moment of inertia for the spars
        I_front_center = parallel_axis_theorem(I_front_center, h_front_center * w_front_center, d)
        I_front_tips = parallel_axis_theorem(I_front_tips, h_front_tip * w_front_tip, d)
        I_rear_center = parallel_axis_theorem(I_rear_center, h_rear_center * w_rear_center, d)
        I_rear_tips = parallel_axis_theorem(I_rear_tips, h_rear_tip * w_rear_tip, d)

        # Adding various components to get the total moment of inertia of the wing box around the x-axis
        I_yy = (I_front_center + 2 * I_front_tips + I_center_center + 2 * I_center_tips + I_rear_center + 2 * I_rear_tips
                + I_top_sheet + I_bottom_sheet)
        A = (l_top * t_top + l_bottom * t_bottom + 2 * w_front_tip * h_front_tip + w_front_center * h_front_center
             + 2 * w_rear_tip * h_rear_tip + w_rear_center * h_rear_center + 2 * w_center_tip * h_center_tip + w_center_center * h_center_center)

    return I_yy, A

# Tool to compute the moment of inertia around the x-axis along the rectangular section of the fuselage
def I_xx_rectangle_section_fuselage(h_top, w_top, h_side, w_side, h_bottom, w_bottom):
    # Moment of inertia for the various rectangles:
    I_sides = moment_of_inertia_rectangle_x_axis(w_side, h_side)
    I_top = moment_of_inertia_rectangle_x_axis(w_top, h_top)
    I_bottom = moment_of_inertia_rectangle_x_axis(w_bottom, h_bottom)

    # Including the parallel axis theorem:
    I_top = parallel_axis_theorem(I_top, w_top * h_top, 0.5 * h_side + 0.5 * h_top)
    I_bottom = parallel_axis_theorem(I_bottom, w_bottom * h_bottom, 0.5 * h_side + 0.5 * h_bottom)

    # Adding various contributions to the total moment of inertia for the fuselage cross-section
    I_xx = 2 * I_sides + I_top + I_bottom

    return I_xx


# Tool to compute the moment of inertia around the y-axis along the rectangular section of the fuselage
def I_yy_rectangle_section_fuselage(h_top, w_top, h_side, w_side, h_bottom, w_bottom):
    # Moment of inertia for the various rectangles:
    I_sides = moment_of_inertia_rectangle_y_axis(w_side, h_side)
    I_top = moment_of_inertia_rectangle_y_axis(w_top, h_top)
    I_bottom = moment_of_inertia_rectangle_y_axis(w_bottom, h_bottom)

    # Including the parallel axis theorem:
    I_sides = parallel_axis_theorem(I_sides, h_side*w_side, 0.5 * w_top - 0.5 * w_side)

    # Adding various contributions to the total moment of inertia for the fuselage cross-section
    I_yy = 2 * I_sides + I_top + I_bottom

    return I_yy

# Torsion calculator for a beam with a certain geometrical cross-section
def compute_torsion(T, rho, J, G, t, A_m, s):
    #Torsion and twist cirucular section
    tau_circ = T * rho / J
    dtheta_dz_circ = T / G / J

    #Torsion and twist thin-walled closed section (assumed constant thickness)
    tau_thin_circ = T / 2 / t / A_m
    dtheta_dz_thin_circ = T / (4 * A_m**2 * G * t)

    # Torsion and twist of a thin plate
    tau_max_thin_plate = 3 * T / s / t**2
    dtheta_dz_thin_plate = 3* T /G /s /t**3

    return tau_circ, dtheta_dz_circ, tau_thin_circ, dtheta_dz_thin_circ, tau_max_thin_plate, dtheta_dz_thin_plate

# Buckling calculator for a beam with a certain geometrical cross-section
def compute_buckling(E, Ixx, buckling, L, v, t, b):
    # Specify buckling type
    if buckling == 'Fixed-Fixed':
        L_e = L /2
        C = 4
    elif buckling == 'Pinned-Fixed':
        L_e = 0.6992* L
        C = 2.046
    elif buckling == 'Fixed-Free':
        L_e = 2*L
        C = 0.25

    # Euler buckling formula
    P_cr = np.pi**2 * E * Ixx / L_e**2

    # Buckling of thin plates
    sigma_cr = C * np.pi**2 * E /12 / (1-v**2) * (t/b)**2

    return P_cr, sigma_cr

# Number of stringers calculator given:
# Source: Airframe structural design, practical design information and data on aircraft structures Michael c. y. Niu
# Spar spacing
# Fixed or free moving
# sheet thickness
# Young modulus for the skin
# Moment at that particular location along the span
# Max sheet height from point in which moment acts


def number_of_stringers_computation(K_c, s, t, E, t_c): #, mu):
    E = E * 0.000145037738
    t = t * 39.3700787
    h = (0.833 - s / 6) * t_c * 39.3700787
    M = momentdistributionz[round(s*100)] * 8.85074579

    N = np.arange(1, 10, 1)
    L = 0.52 * (0.833 - s / 6) * 39.3700787 / N
    print(N, 'N')

    # Compute the maximum compressive force acting on the wing box section
    Fmax = M / h
    sigma_cr = Fmax / (t * L)
    L_check = t / np.sqrt(sigma_cr / (K_c * E))
    #L_check = t / np.sqrt(12 * sigma_cr * (1 - mu**2)/(K_c * E * np.pi**2))
    # Compute the maximum spacing between stringers
    N_req = L/L_check
    print(N_req, "Nreq")
    # return N_stringers

def number_of_stringers_computation_horizontal_tail(K_c, s, t, E, t_c): #, mu):
    E = E * 0.000145037738
    t = t * 39.3700787
    h = 0.118 * t_c * 39.3700787
    M = momentdistributionz[round(s*100)] * 8.85074579

    N = np.arange(1, 5, 1)
    L = 0.5 * 0.118 * 39.3700787 / N
    print(N, 'N')

    # Compute the maximum compressive force acting on the wing box section
    Fmax = M / h
    sigma_cr = Fmax / (t * L)
    L_check = t / np.sqrt(sigma_cr / (K_c * E))
    #L_check = t / np.sqrt(12 * sigma_cr * (1 - mu**2)/(K_c * E * np.pi**2))
    # Compute the maximum spacing between stringers
    N_req = L/L_check
    print(N_req, "Nreq")
    # return N_stringers

def number_of_stringers_computation_vertical_tails(K_c, s, t, E, t_c, A): #, mu):
    E = E * 0.000145037738
    t = t * 39.3700787
    h = (0.43 - 0.2443 * s) * t_c * 39.3700787
    M = momentdistributionz[round(s*100)] * 8.85074579

    N = np.arange(1, 5, 1)
    L = 0.5 * (0.43 - 0.2443 * s) * 39.3700787 / N
    print(N, 'N')

    # Compute the maximum compressive force acting on the wing box section
    Fmax = M / h #+ 109.722 / 2 / A
    sigma_cr = Fmax / (t * L)
    L_check = t / np.sqrt(sigma_cr / (K_c * E))
    #L_check = t / np.sqrt(12 * sigma_cr * (1 - mu**2)/(K_c * E * np.pi**2))
    # Compute the maximum spacing between stringers
    N_req = L/L_check
    print(N_req, "Nreq")
    # return N_stringers
def number_of_ribs_buckling_computation(K_c, s, t, E, T, t_c):
    E = E * 0.000145037738
    t = t * 39.3700787
    h = (0.833 - s / 6) * t_c * 39.3700787
    T = T * 8.85074579

    N = np.arange(1, 10, 1)
    L = 0.52 * (0.833 - s / 6) * 39.3700787 / N
    print(N, 'N')

    # Compute the maximum compressive force acting on the wing box section
    Fmax = T / h
    sigma_cr = Fmax / (t * L)
    L_check = t / np.sqrt(12 * sigma_cr * (1 - mu**2)/(K_c * E * np.pi**2))
    # Compute the maximum spacing between stringers
    N_req = L / L_check
    print(N_req, "Nreq")
    # return N_stringers

# Main wing stringer computations:
# print("Compute number of stringers required at s = 0:", number_of_stringers_computation(3.62, 0, 0.001, 71*10**9,
#                                                                                         0.12))#, 0.3))
# print("Compute number of stringers required at s = 0.273:", number_of_stringers_computation(3.62, 0.273, 0.001, 71*10**9,
#                                                                                         0.12))#, 0.3))
# print("Compute number of stringers required at s = 0.836:", number_of_stringers_computation(3.62, 0.836, 0.001, 71*10**9,
#                                                                                         0.12))#, 0.3))
# print("Compute number of stringers required at s = 1.15:", number_of_stringers_computation(3.62, 1.15, 0.001, 71*10**9,
#                                                                                         0.12))#, 0.3))
# print("Compute number of stringers required at s = 2.25:", number_of_stringers_computation(3.62, 2.25, 0.001, 71*10**9,
#                                                                                         0.12))#, 0.3))
# print("Compute number of stringers required at s = 2.9:", number_of_stringers_computation(3.62, 2.9, 0.001, 71*10**9,
#                                                                                         0.12))#, 0.3))

# Stringer computations for the horizontal tail:
# print("The moment of inertia of the wing box for the horizontal tail is:",
#       I_xx_wing_box(0.008, 0.003, 0.002,0.004,0.008,0.002,
#                     0.004,0.008,0,0,0,0,
#                     0.001,0.5*0.118,0.001,0.5*0.118, 2))
print("Compute the number of stringers required for the horizontal tail at s= 0:",
      number_of_stringers_computation_horizontal_tail(3.62, 0, 0.001, 71*10**9, 0.12))
print("Compute the number of stringers required for the horizontal tail at s= 0.4:",
      number_of_stringers_computation_horizontal_tail(3.62, 0.5625, 0.001, 71*10**9, 0.12))
print("Compute the number of stringers required for the horizontal tail at s= 1.0:",
      number_of_stringers_computation_horizontal_tail(3.62, 1.15, 0.001, 71*10**9, 0.12))
# print("Compute the number of stringers required for the horizontal tail at s= 1.3:",
#       number_of_stringers_computation_horizontal_tail(3.62, 1.3, 0.001, 71*10**9, 0.12))
# print("Compute the number of stringers required for the horizontal tail at s= 1.9:",
#       number_of_stringers_computation_horizontal_tail(3.62, 1.9, 0.001, 71*10**9, 0.12))

# Computations for the vertical tails:
# print("The moment of inertia of the wing box for the vertical tails is:",
#       I_xx_wing_box(0.035*0.75,0.005*0.75,0.007*0.75,0.01*0.75,0.035*0.75,0.005*0.75,0.007*0.75,0.01*0.75,
#                      0,0,0,0,0.001,0.2*0.75,0.001,0.2*0.75,2))
# print("Compute the number of stringers required for the vertical tails at s = 0:",
#       number_of_stringers_computation_vertical_tails(3.62, 0, 0.001, 70* 10**9, 0.12, 0.00103))
# print("Compute the number of stringers required for the vertical tails at s = 0.44:",
#       number_of_stringers_computation_vertical_tails(3.62, 0.44, 0.001, 70* 10**9, 0.12, 0.000654375))


# Check for valid moment of inertia of the wing box for the horizontal tail
# print("The sizing of the moment of inertia for the horizontal tail wing is:"
#       , I_xx_wing_box(0.01, 0.002,0.004,0.008,0.01,0.002,
#                       0.004,0.008,0,0,0,0,0.001,
#                       0.109*0.5, 0.001, 0.109*0.5, 2))
