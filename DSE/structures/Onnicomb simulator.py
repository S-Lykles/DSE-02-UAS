import numpy as np
import matplotlib.pyplot as plt

#INPUTS
#------
import math

def calculate_bending_stress(moment, section_modulus):
    """
    Calculate bending stress.

    Parameters:
    - moment: Bending moment (Nm)
    - section_modulus: Section modulus of the cross-section (m^3)

    Returns:
    - Bending stress (Pa)
    """
    return moment / section_modulus

def calculate_shear_stress(shear_force, shear_area):
    """
    Calculate shear stress.

    Parameters:
    - shear_force: Shear force (N)
    - shear_area: Shear area of the cross-section (m^2)

    Returns:
    - Shear stress (Pa)
    """
    return shear_force / shear_area

def calculate_normal_stress(axial_force, area):
    """
    Calculate normal stress.

    Parameters:
    - axial_force: Axial force (N)
    - area: Cross-sectional area (m^2)

    Returns:
    - Normal stress (Pa)
    """
    return axial_force / area

def main():
    # Example input values (replace with your specific values)
    bending_moment = 1000.0  # Nm
    shear_force = 500.0  # N
    axial_force = 2000.0  # N
    section_modulus = 0.005  # m^3
    shear_area = 0.002  # m^2
    cross_section_area = 0.01  # m^2

    # Calculate stresses
    bending_stress = calculate_bending_stress(bending_moment, section_modulus)
    shear_stress = calculate_shear_stress(shear_force, shear_area)
    normal_stress = calculate_normal_stress(axial_force, cross_section_area)

    # Print results
    print(f"Bending Stress: {bending_stress} Pa")
    print(f"Shear Stress: {shear_stress} Pa")
    print(f"Normal Stress: {normal_stress} Pa")

if __name__ == "__main__":
    main()

