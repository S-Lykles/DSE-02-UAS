from DSE import const
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

#Dimensions leading rib
h_leading_rib_center =
w_leading_rib_center =
m_leading_rib_center =
h_leading_rib_tips =
w_leading_rib_tips =
m_leading_rib_tips =

#Dimensions trailing rib
h_trailing_rib_center = h_leading_rib_center
w_trailing_rib_center = w_leading_rib_center

h_trailing_rib_tips = h_leading_rib_tips
w_trailing_rib_tips = w_leading_rib_tips

#Dimensions top sheet
t_top_sheet =
l_top_sheet =

#Dimensions bottom sheet
t_bottom_sheet =
l_bottom_sheet =

# Function computing moment of inertia around x-axis with cg on neutral line
def compute_moment_of_inertia_object_x_axis(width, height):
    I_xx = width * height**3 / 12

    return I_xx


# Function computing moment of inertia around y-axis with cg on neutral line
def compute_moment_of_inertia_object_y_axis(width, height):
    I_yy = height * width**3 / 12

    return I_yy


# Parallel axis theorem
def parallel_axis_theorem(inertia, mass, distance):
    I = inertia + mass*distance**2

    return I

# Function for final moment of inertia for the wing box around the x-axis
def I_xx_wing_box(h_leading_rib_center, w_leading_rib_center, h_leading_rib_tips, w_leading_rib_tips, h_trailing_rib_center,
    w_trailing_rib_center, h_trailing_rib_tips, w_trailing_rib_tips, t_top_sheet, l_top_sheet, t_bottom_sheet, l_bottom_sheet):

    # Moment of inertia of various rectangles in the
    I_leading_rip_center = compute_moment_of_inertia_object_x_axis(w_leading_rib_center,h_leading_rib_center)
    I_leading_rip_tips = compute_moment_of_inertia_object_x_axis(w_leading_rib_tips, h_leading_rib_tips)
    I_trailing_rip_center = compute_moment_of_inertia_object_x_axis(w_trailing_rib_center, h_trailing_rib_center)
    I_trailing_rip_tips = compute_moment_of_inertia_object_x_axis(w_trailing_rib_tips, h_trailing_rib_tips)
    I_top_sheet = compute_moment_of_inertia_object_x_axis(l_top_sheet, t_top_sheet)
    I_bottom_sheet = compute_moment_of_inertia_object_x_axis(l_bottom_sheet, t_bottom_sheet)

    # Including the parallel axis theorem
    I_leading_rip_tips = parallel_axis_theorem(I_leading_rip_tips, )

