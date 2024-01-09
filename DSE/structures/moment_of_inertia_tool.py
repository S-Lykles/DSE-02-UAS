from DSE import const
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

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

#Function computing the moment of inertia of a solid circular section
def moment_of_inertia_solid_circular_section(diameter)
    I = np.pi * diameter**4 / 64

    return I

#Function computing the moment of inertia of a thin-walled circular section
def moment_of_inertia_thin_walled_circular_section(diameter, thickness)
    I = np.pi * thickness * diameter**3 / 8

    return I

# Parallel axis theorem
def parallel_axis_theorem(inertia, mass, distance):
    I = inertia + mass*distance**2

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

def I_xx_wing_box_two_spars(h_front_center, w_front_center, h_front_tip, w_front_tip, h_rear_center,
                            w_rear_center, h_rear_tip, w_rear_tip, t_top, l_top, t_bottom, l_bottom, number_of_spars):
    if number_of_spars == 1:

        I_spar_center = moment_of_inertia_rectangle_x_axis(w_front_center, h_front_center)
        I_spar_tips = moment_of_inertia_rectangle_x_axis(w_front_tip, h_front_tip)
        I_top_sheet = moment_of_inertia_rectangle_x_axis(l_top, t_top)
        I_bottom_sheet = moment_of_inertia_rectangle_x_axis(l_bottom, t_bottom)

        # Including the parallel axis theorem

        # Moment of inertia for tips of spar
        I_spar_tips = parallel_axis_theorem(I_spar_tips, h_front_tip * w_front_tip,
                                             0.5 * h_front_center + 0.5 * h_front_tip)
        # Moment of inertia top and bottom sheet
        I_top_sheet = parallel_axis_theorem(I_top_sheet, t_top * l_top, (
                0.5 * h_front_center + h_front_tip + 2 * t_top + 0.5 * h_rear_center + h_rear_tip) / 2)
        I_bottom_sheet = parallel_axis_theorem(I_bottom_sheet, t_bottom * l_bottom,
                                               (
                                                           0.5 * h_front_center + h_front_tip + t_bottom + 0.5 * h_rear_center + h_rear_tip) / 2)

    elif number_of_spars == 2:
        # Moment of inertia of various rectangles in the
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
                    0.5 * h_front_center + h_front_tip + 2 * t_top + 0.5 * h_rear_center + h_rear_tip) / 2)
        I_bottom_sheet = parallel_axis_theorem(I_bottom_sheet, t_bottom*l_bottom,
                                               (0.5 * h_front_center + h_front_tip + t_bottom + 0.5 * h_rear_center + h_rear_tip) / 2)

        # Adding various components to get the total moment of inertia of the wing box around the x-axis
        I_xx = I_front_center + 2 * I_front_tips + I_rear_center + 2 * I_rear_tips + I_top_sheet + I_bottom_sheet

    return I_xx


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

def compute_buckling(E, Ixx, buckling, L, v, t, b)
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

