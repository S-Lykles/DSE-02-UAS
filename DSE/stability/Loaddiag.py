import numpy as np
from DSE import const
from DSE.structures.center_of_gravity import class_two_cg_estimation



def load_diagram_plot(empty, fuel, payload_supply):
    """Create potato plot based on Cg calculation
    """
    fuel_frac = np.linspace(0.0, 1.0, 101)
    payload_frac = np.linspace(0.0, 1.0, 101)


load_diagram_plot()
