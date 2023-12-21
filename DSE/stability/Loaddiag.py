import numpy as np
from DSE import const
from DSE.structures.center_of_gravity import class_two_cg_estimation
import matplotlib.pyplot as plt
from DSE import plot_setting



def load_diagram_plot(empty, fuel, payload_supply, plots == False):
    """Create potato plot based on Cg calculation

    needs OEW
    needs fuel weight (fuel location estimation based on class II estimation)
    needs payload weight (weight location based on class II estimation, assumed to be concentrated in c_g of payload)"""

    weight = []
    c_g = []

    #c_g variation for fuel
    for i in range(101):
        weight.append(class_two_cg_estimation(empty, i/100 * fuel)[0])
        c_g.append(class_two_cg_estimation(empty, i/100 * fuel)[1,0])

    #c_g variation for payload, fuel included
    for i in range(101):
        weight.append(class_two_cg_estimation(empty, fuel, i/100 * payload_supply)[0])
        c_g.append(class_two_cg_estimation(empty, fuel, i/100 * payload_supply)[1,0])

    if plots == True:
        plt.plot(c_g, weight)
        plt.show()

    return weight, c_g
